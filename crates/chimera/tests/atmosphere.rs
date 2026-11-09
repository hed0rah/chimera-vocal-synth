//! The atmosphere model against measured gas properties, real-fluid values, and the one acoustic
//! effect of a gas change with a long record in speech science: helium speech.

use chimera::airway::frequency::{find_peaks, transfer_function};
use chimera::anatomy::human::{BODY_MALE, NEUTRAL};
use chimera::atmosphere::{self, Breather, Conditions, DRY_AIR, PRESETS, Species, state};
use chimera::engine::Engine;
use chimera::physics::Gas;
use chimera::voice::human_airway::HumanAirway;

fn close(got: f64, want: f64, tol: f64) -> bool {
    (got / want - 1.0).abs() <= tol
}

#[test]
fn earth_breathed_by_a_human_is_the_tuned_air_exactly() {
    assert_eq!(Conditions::EARTH.gas(), Gas::BODY_AIR);
    let p = atmosphere::preset("earth").unwrap();
    assert_eq!(p.conditions(Breather::Human).gas(), Gas::BODY_AIR);
}

/// Handbook values at 20 °C, 1 atm: dry air c 343.2 m/s, rho 1.204 kg/m³, mu 18.1 µPa s (ISO
/// 9613-1; CRC Handbook); helium c 1007 m/s; SF6 c between 130 and 140 m/s (NIST lists about 134).
#[test]
fn gases_at_room_temperature_match_handbook_values() {
    let air = state(&DRY_AIR, 293.15, 101325.0);
    assert!(close(air.c, 343.2, 0.002), "air c {}", air.c);
    assert!(close(air.rho, 1.204, 0.002), "air rho {}", air.rho);
    assert!(close(air.mu, 1.81e-5, 0.02), "air mu {}", air.mu);
    let mut he = [0.0; atmosphere::N_SPECIES];
    he[Species::He as usize] = 1.0;
    assert!(close(state(&he, 293.15, 101325.0).c, 1007.0, 0.01));
    let mut sf6 = [0.0; atmosphere::N_SPECIES];
    sf6[Species::SF6 as usize] = 1.0;
    let c = state(&sf6, 293.15, 101325.0).c;
    assert!((130.0..140.0).contains(&c), "SF6 c {c}");
}

/// Real-fluid c, rho and mu (CoolProp, GERG-style mixing; Wilke mixing for mu) for the same
/// conditions as the presets, from the databurn planetary-acoustics model (venus_saturn.json).
/// The ideal model's errors are largest where the gas is least ideal: Venus at 92 bar, cold
/// Titan; and its Sutherland viscosity is poor for hydrogen below 150 K (Saturn, Neptune).
#[test]
fn the_ideal_model_tracks_real_fluid_values_at_every_site() {
    let real: &[(&str, f64, f64, f64, f64)] = &[
        // id, c m/s, rho kg/m³, mu Pa s, mu tolerance
        ("venus_surface", 421.05, 64.759, 3.3806e-5, 0.10),
        ("venus_55km", 271.58, 0.86776, 1.5046e-5, 0.10),
        ("mars_jezero", 242.36, 0.013603, 1.1751e-5, 0.10),
        ("mars_hellas", 251.53, 0.027952, 1.2734e-5, 0.10),
        ("mars_olympus", 211.48, 0.0028219, 8.6219e-6, 0.10),
        ("titan", 194.97, 5.3529, 6.3953e-6, 0.10),
        ("saturn", 855.12, 0.20622, 5.8019e-6, 0.25),
        ("neptune", 612.45, 0.43293, 4.0027e-6, 0.25),
        ("pluto", 124.11, 1.0450e-4, 2.1833e-6, 0.10),
    ];
    for &(id, c, rho, mu, mu_tol) in real {
        let p = atmosphere::preset(id).unwrap();
        let s = state(&p.mix, p.temperature, p.pressure);
        assert!(close(s.c, c, 0.025), "{id}: c {} vs {c}", s.c);
        assert!(close(s.rho, rho, 0.04), "{id}: rho {} vs {rho}", s.rho);
        assert!(close(s.mu, mu, mu_tol), "{id}: mu {} vs {mu}", s.mu);
    }
}

#[test]
fn the_airway_boils_below_the_armstrong_limit() {
    let mars = atmosphere::preset("mars_jezero").unwrap();
    let mut c = mars.conditions(Breather::Human);
    assert_eq!(c.airway_humidity, 0.0, "a human preset below 6.3 kPa breathes dry");
    c.airway_humidity = 1.0;
    let a = c.airway();
    assert!(a.boiling && a.water == 1.0);
    assert!(!Conditions::EARTH.airway().boiling);
    let armstrong = atmosphere::water_saturation(atmosphere::BODY_TEMPERATURE);
    assert!((6200.0..6400.0).contains(&armstrong), "Armstrong limit {armstrong} Pa");
}

fn formants(gas: Gas) -> Vec<f64> {
    let mut air = HumanAirway::new_in(96000.0, BODY_MALE, gas);
    air.apply(&NEUTRAL, 0.0);
    let el = air.net.linear_elements();
    let freqs: Vec<f64> = (0..1400).map(|i| 100.0 * (9000.0f64 / 100.0).powf(i as f64 / 1399.0)).collect();
    let h = transfer_function(&el, air.pharynx.first_node(), &freqs, &[air.glottis_port]);
    find_peaks(&freqs, &h.mag, 150.0).iter().map(|p| p.freq).collect()
}

/// Helium speech: the resonances rise with the speed of sound, but F1 rises less, because the
/// yielding walls give the closed tract a resonance of its own that does not scale with c (Fant
/// and Lindqvist 1968, STL-QPSR 9(1); Richards 1982).
#[test]
fn in_helium_the_formants_rise_with_c_and_f1_rises_least() {
    let mut mix = [0.0; atmosphere::N_SPECIES];
    mix[Species::He as usize] = 0.8;
    mix[Species::O2 as usize] = 0.2;
    let heliox = Conditions { mix, ..Conditions::EARTH }.airway();
    let ratio = heliox.gas.c / Gas::BODY_AIR.c;
    assert!(ratio > 1.6, "heliox c ratio {ratio}");
    let air = formants(Gas::BODY_AIR);
    let he = formants(heliox.gas);
    assert!(air.len() >= 3 && he.len() >= 3, "{air:?} {he:?}");
    // match each air formant to the helium peak nearest its c-scaled position
    let shifted: Vec<f64> = air[..3].iter().map(|f| he.iter().copied().min_by(|a, b| (a - f * ratio).abs().total_cmp(&(b - f * ratio).abs())).unwrap()).collect();
    let r: Vec<f64> = air[..3].iter().zip(&shifted).map(|(a, h)| h / a).collect();
    assert!(close(r[2], ratio, 0.06), "F3 ratio {} vs c ratio {ratio}", r[2]);
    assert!(r[0] < r[2] - 0.05, "F1 should rise least: ratios {r:?}, c ratio {ratio}");
}

/// Every preset, breathed by a human and by a native creature, renders a finite, audible vowel.
#[test]
fn every_preset_renders_a_stable_voice() {
    for p in PRESETS {
        for who in [Breather::Human, Breather::Native] {
            let mut e = Engine::new(48000.0);
            e.set_gas(p.conditions(who).gas());
            e.set_gate(1.0);
            let mut out = vec![0.0f32; 48000];
            e.render(&mut out);
            let peak = out.iter().fold(0.0f32, |m, x| m.max(x.abs()));
            assert!(out.iter().all(|x| x.is_finite()), "{} {who:?}: not finite", p.id);
            assert!(peak > 1e-6, "{} {who:?}: silent", p.id);
        }
    }
}

/// Every creature, in every preset as a native of the place, sings a finite, audible call.
#[test]
fn every_creature_is_stable_in_every_preset() {
    use chimera::creatures::{CREATURES, CreatureVoice};
    chimera::creatures::init();
    for def in CREATURES.iter() {
        for p in PRESETS {
            let mut v = CreatureVoice::new_in(def, 48000.0, p.conditions(Breather::Native).gas(), &[]);
            v.set_gate(1.0);
            v.play_call(0, None, 1.0);
            let mut out = vec![0.0f32; 24000];
            v.render(&mut out);
            let peak = out.iter().fold(0.0f32, |m, x| m.max(x.abs()));
            assert!(out.iter().all(|x| x.is_finite()), "{} in {}: not finite", def.id, p.id);
            assert!(peak > 1e-7, "{} in {}: silent", def.id, p.id);
        }
    }
}
