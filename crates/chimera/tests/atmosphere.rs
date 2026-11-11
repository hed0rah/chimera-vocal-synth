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

/// Real-fluid c, rho and mu at each preset's outside conditions, computed with CoolProp: its HEOS
/// mixture backend (GERG-style mixing) for c and rho, Wilke's rule over its pure-fluid
/// viscosities for mu. Sutherland's viscosity is poor for hydrogen below 150 K (Saturn, Neptune),
/// hence the wider tolerance on mu there.
#[test]
fn the_model_tracks_real_fluid_values_at_every_site() {
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

/// Speech in compressed air: c does not change with pressure, but the denser gas stiffens against the
/// yielding walls and raises their closed-tract resonance F_w, which adds in quadrature to every formant,
/// F^2 = F_rigid^2 + F_w^2, with F_w^2 proportional to the density (Fant and Sonesson 1964, STL-QPSR 5(2);
/// Fant and Lindqvist 1968). So F1 moves most, every formant gains about the same F^2, and that gain grows
/// in proportion to the pressure added.
#[test]
fn in_compressed_air_the_walls_add_the_same_f_squared_to_every_formant() {
    let at = |bar: f64| formants(Conditions { pressure: bar * atmosphere::ATM, airway_humidity: 0.0, ..Conditions::EARTH }.gas());
    let (one, two, eleven) = (at(1.0), at(2.0), at(11.0));
    assert!(one.len() >= 3 && two.len() >= 3 && eleven.len() >= 3);
    let gain = |x: &[f64]| -> Vec<f64> { (0..3).map(|k| x[k] * x[k] - one[k] * one[k]).collect() };
    let (g2, g11) = (gain(&two), gain(&eleven));
    let mean = |g: &[f64]| g.iter().sum::<f64>() / g.len() as f64;
    for g in [&g2, &g11] {
        for k in 0..3 {
            assert!(close(g[k], mean(g), 0.15), "F{} gains {} of F^2, the mean is {}: {g:?}", k + 1, g[k], mean(g));
        }
    }
    // 10 bar added against 1: ten times the gain
    assert!(close(mean(&g11) / mean(&g2), 10.0, 0.1), "gain ratio {}", mean(&g11) / mean(&g2));
    assert!(eleven[0] / one[0] > 1.5 && eleven[2] / one[2] < 1.1, "F1 moves most: {one:?} -> {eleven:?}");
}

/// Peng-Robinson against CoolProp's reference equations of state (pure fluids): within 0.5 % in dilute gas,
/// a few % in dense gas, and 15 % next to a critical point, which the airway reports as near_critical.
#[test]
fn the_real_gas_model_tracks_reference_equations_of_state() {
    // species, T K, P Pa, rho kg/m3, c m/s (CoolProp 8.0), tolerance on each
    let refs: &[(Species, f64, f64, f64, f64, f64)] = &[
        (Species::CO2, 737.0, 9.2e6, 65.7008, 417.26, 0.01),     // Venus's surface
        (Species::CO2, 230.0, 600.0, 0.0138, 240.00, 0.005),      // Mars
        (Species::N2, 93.65, 146700.0, 5.4712, 193.36, 0.005),    // Titan
        (Species::He, 303.15, 3.1e6, 4.8525, 1038.15, 0.01),     // a heliox dive
        (Species::H2, 134.0, 1.0e5, 0.1809, 915.05, 0.005),       // Saturn's 1 bar level
        (Species::SF6, 293.15, 101325.0, 6.1453, 133.85, 0.005),
        (Species::CO2, 400.0, 2.0e7, 380.4992, 310.75, 0.05),     // dense
        (Species::CH4, 250.0, 1.0e7, 112.4258, 392.22, 0.05),
        (Species::Ar, 300.0, 2.0e7, 335.8421, 372.04, 0.05),
        (Species::N2, 300.0, 3.0e7, 294.7812, 483.18, 0.08),
        (Species::CO2, 310.15, 9.2e6, 631.2603, 258.59, 0.15),   // 1.02 Tc, 1.25 Pc: a human airway on Venus
    ];
    for &(sp, t, p, rho, c, tol) in refs {
        let mut m = [0.0; atmosphere::N_SPECIES];
        m[sp as usize] = 1.0;
        let st = state(&m, t, p);
        assert!(close(st.rho, rho, tol), "{sp:?} {t} K {p} Pa: rho {} vs {rho}", st.rho);
        assert!(close(st.c, c, tol), "{sp:?} {t} K {p} Pa: c {} vs {c}", st.c);
    }
    // the ideal gas was 4x too thin there; the airway says it is near the critical point
    let venus = atmosphere::preset("venus_surface").unwrap().conditions(Breather::Human).airway();
    assert!(venus.near_critical && venus.state.compressibility < 0.4);
    assert!(!Conditions::EARTH.airway().near_critical);
}

const RATES: [f64; 10] = [8000.0, 16000.0, 22050.0, 24000.0, 32000.0, 44100.0, 48000.0, 88200.0, 96000.0, 192000.0];

/// Pure hydrogen at 800 K, the panel's hottest airway: the fastest gas it can make.
fn hottest_hydrogen() -> Gas {
    let mut mix = [0.0; atmosphere::N_SPECIES];
    mix[Species::H2 as usize] = 1.0;
    Conditions { mix, pressure: atmosphere::ATM, temperature: 800.0, airway_temperature: 800.0, airway_humidity: 0.0 }.gas()
}

/// Every body at every rate, in air, in every preset and in the fastest gas, builds an airway under
/// the Courant limit; in air at 44.1 kHz and above the oversampling is upstream's (the goldens stand).
#[test]
fn every_airway_is_built_under_the_courant_limit() {
    use chimera::anatomy::human::{BODY_CHILD, BODY_FEMALE};
    use chimera::physics::MAX_CFL;
    use chimera::voice::human_voice::{airway_for, oversample_for};
    let mut gases = vec![("air", Gas::BODY_AIR), ("hottest hydrogen", hottest_hydrogen())];
    for p in PRESETS {
        gases.push((p.id, p.conditions(Breather::Human).gas()));
        gases.push((p.id, p.conditions(Breather::Native).gas()));
    }
    for (name, body) in [("male", BODY_MALE), ("female", BODY_FEMALE), ("child", BODY_CHILD)] {
        for fs in RATES {
            for &(id, gas) in &gases {
                let (m, a) = airway_for(fs, body, oversample_for(fs), gas);
                let cfl = a.net.check_cfl();
                assert!(cfl <= MAX_CFL, "{name} at {fs} Hz in {id}: CFL {cfl:.3} at m {m}");
                if id == "air" && fs >= 44100.0 {
                    assert_eq!(m, oversample_for(fs), "{name} at {fs} Hz in air");
                }
            }
        }
    }
}

/// The child voice at a Bluetooth headset's rates sounds as it does at 48 kHz. Its piriform fossa
/// (three sections) was at a Courant number of 1.11 at 8, 16 and 32 kHz: the voice ran at the
/// limiter ceiling, peak 0.75 and three times the RMS.
#[test]
fn the_child_voice_holds_at_low_rates() {
    use chimera::anatomy::human::BodyMix;
    let held = |fs: f64| {
        let mut e = Engine::new(fs);
        e.set_body(BodyMix::new(0.0, 1.0));
        e.set_gate(1.0);
        let mut out = vec![0.0f32; 2 * fs as usize];
        for chunk in out.chunks_mut(128) {
            e.render(chunk);
        }
        let tail = &out[fs as usize..];
        (tail.iter().map(|&x| (x as f64) * (x as f64)).sum::<f64>() / tail.len() as f64).sqrt()
    };
    let reference = held(48000.0);
    for fs in [8000.0, 16000.0, 22050.0, 24000.0, 32000.0] {
        let rms = held(fs);
        assert!(rms > 0.7 * reference && rms < 1.4 * reference, "{fs} Hz: RMS {rms:.3}, 48 kHz {reference:.3}");
    }
}

/// Every creature at an 8 kHz output in the fastest gas builds under the Courant limit and sings.
#[test]
fn every_creature_holds_in_the_fastest_gas_at_a_low_rate() {
    use chimera::creatures::{CREATURES, CreatureVoice};
    use chimera::physics::MAX_CFL;
    chimera::creatures::init();
    let gas = hottest_hydrogen();
    for def in CREATURES.iter() {
        let mut v = CreatureVoice::new_in(def, 8000.0, gas, &[]);
        if let Some(n) = v.sim.net() {
            assert!(n.check_cfl() <= MAX_CFL, "{}: CFL {:.3} at m {}", def.id, n.check_cfl(), v.m);
        }
        v.set_gate(1.0);
        v.play_call(0, None, 1.0);
        let mut out = vec![0.0f32; 8000];
        v.render(&mut out);
        assert!(out.iter().all(|x| x.is_finite()), "{}: not finite", def.id);
    }
}

/// An empty mix has no state; the engine refuses it and keeps sounding in the gas it had.
#[test]
fn the_engine_refuses_a_gas_that_is_not_physical() {
    let empty = Conditions { mix: [0.0; atmosphere::N_SPECIES], pressure: atmosphere::ATM, temperature: 288.15, airway_temperature: 288.15, airway_humidity: 0.0 };
    let g = empty.gas();
    assert!(!g.is_physical(), "an empty mix gave {g:?}");
    let mut e = Engine::new(48000.0);
    let before = e.info();
    let after = e.set_gas(g);
    assert_eq!((before.n_ph, before.n_or), (after.n_ph, after.n_or));
    e.set_gate(1.0);
    let mut out = vec![0.0f32; 24000];
    e.render(&mut out);
    assert!(out.iter().all(|x| x.is_finite()) && out.iter().any(|x| x.abs() > 1e-3));
}
