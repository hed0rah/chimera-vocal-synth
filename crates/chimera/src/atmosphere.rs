//! Atmospheres: the gas an airway breathes, worked out from physical conditions.
//!
//! The controls are the conditions themselves: a mix of gases (mole fractions), the pressure,
//! the temperature of the air, and the temperature and humidity it reaches inside the airway. From
//! those, for an ideal gas:
//!
//! - speed of sound c = sqrt(gamma R T / M), with gamma from each species' heat capacity at T
//!   (translation and rotation, harmonic vibrational modes, and hydrogen's quantised rotation);
//! - density rho = P M / (R T);
//! - viscosity mu from Sutherland's law per species, mixed by Wilke's rule.
//!
//! The airway was tuned in [`Gas::BODY_AIR`]. A condition is applied as a ratio to the same model's
//! own warm, saturated air ([`reference`]), so the model's small errors cancel and the reference
//! conditions give `BODY_AIR` exactly. Against real-fluid values (CoolProp) at ten planetary
//! sites, the ideal model is within 2.5 % on c and 4 % on rho; see `tests/atmosphere.rs`.
//!
//! Units here are SI (K, Pa, kg/m³, m/s, Pa·s); [`Conditions::gas`] converts to the engine's CGS.

use crate::math::{exp, js_max, js_min, sqrt};
use crate::physics::Gas;

/// Molar gas constant, J/(mol K).
pub const R: f64 = 8.314462618;
/// 37 °C, the airway of a warm-blooded body (K).
pub const BODY_TEMPERATURE: f64 = 310.15;
/// One standard atmosphere (Pa).
pub const ATM: f64 = 101325.0;

/// The gases a mix is made of.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Species {
    N2,
    O2,
    Ar,
    CO2,
    He,
    H2,
    CH4,
    H2O,
    SF6,
}

pub const N_SPECIES: usize = 9;
/// Mole fractions, indexed by `Species as usize`; normalised when used.
pub type Mix = [f64; N_SPECIES];

/// What the model needs to know about one gas.
#[derive(Clone, Copy, Debug)]
pub struct SpeciesData {
    pub species: Species,
    /// formula, as the UI shows it
    pub formula: &'static str,
    pub name: &'static str,
    /// kg/mol
    pub molar_mass: f64,
    /// heat capacity cp/R with every vibration frozen: 2.5 monatomic, 3.5 linear, 4 nonlinear
    /// (hydrogen: 2.5, its rotation is added from its quantised levels)
    pub frozen_cp: f64,
    /// harmonic vibrational modes: (characteristic temperature K, degeneracy)
    pub modes: &'static [(f64, f64)],
    /// Sutherland's law: mu0 (Pa s) at T0 (K), constant S (K)
    pub sutherland: (f64, f64, f64),
}

/// Molar masses from IUPAC atomic weights; vibrational temperatures from the fundamental
/// wavenumbers (theta = h c nu / k = 1.4388 K cm x nu); Sutherland constants from White,
/// Viscous Fluid Flow (N2, O2, Ar, CO2, H2, CH4, H2O). He, SF6: fitted to tabulated viscosities
/// near room temperature (MODELLED: rough far from 300 K).
pub const SPECIES: [SpeciesData; N_SPECIES] = [
    SpeciesData { species: Species::N2, formula: "N2", name: "nitrogen", molar_mass: 28.0134e-3, frozen_cp: 3.5, modes: &[(3390.0, 1.0)], sutherland: (1.663e-5, 273.0, 107.0) },
    SpeciesData { species: Species::O2, formula: "O2", name: "oxygen", molar_mass: 31.9988e-3, frozen_cp: 3.5, modes: &[(2256.0, 1.0)], sutherland: (1.919e-5, 273.0, 139.0) },
    SpeciesData { species: Species::Ar, formula: "Ar", name: "argon", molar_mass: 39.948e-3, frozen_cp: 2.5, modes: &[], sutherland: (2.125e-5, 273.0, 144.0) },
    SpeciesData { species: Species::CO2, formula: "CO2", name: "carbon dioxide", molar_mass: 44.0095e-3, frozen_cp: 3.5, modes: &[(960.0, 2.0), (1890.0, 1.0), (3360.0, 1.0)], sutherland: (1.370e-5, 273.0, 222.0) },
    SpeciesData { species: Species::He, formula: "He", name: "helium", molar_mass: 4.002602e-3, frozen_cp: 2.5, modes: &[], sutherland: (1.87e-5, 273.0, 79.4) },
    SpeciesData { species: Species::H2, formula: "H2", name: "hydrogen", molar_mass: 2.01588e-3, frozen_cp: 2.5, modes: &[(6332.0, 1.0)], sutherland: (8.411e-6, 273.0, 97.0) },
    SpeciesData { species: Species::CH4, formula: "CH4", name: "methane", molar_mass: 16.0428e-3, frozen_cp: 4.0, modes: &[(1879.0, 3.0), (2207.0, 2.0), (4197.0, 1.0), (4344.0, 3.0)], sutherland: (1.03e-5, 273.0, 198.0) },
    SpeciesData { species: Species::H2O, formula: "H2O", name: "water vapour", molar_mass: 18.01528e-3, frozen_cp: 4.0, modes: &[(2295.0, 1.0), (5262.0, 1.0), (5404.0, 1.0)], sutherland: (1.12e-5, 350.0, 1064.0) },
    SpeciesData { species: Species::SF6, formula: "SF6", name: "sulfur hexafluoride", molar_mass: 146.0554e-3, frozen_cp: 4.0, modes: &[(1115.1, 1.0), (925.1, 2.0), (1364.0, 3.0), (884.9, 3.0), (753.9, 3.0), (500.7, 3.0)], sutherland: (1.42e-5, 273.0, 240.0) },
];

/// Dry air at the surface: N2, O2, Ar and 420 ppm CO2 (the rest, under 0.003 %, left out).
pub const DRY_AIR: Mix = [0.78084, 0.20946, 0.00934, 0.00042, 0.0, 0.0, 0.0, 0.0, 0.0];

/// Hydrogen's rotational temperature (K).
const THETA_ROT_H2: f64 = 85.35;

/// Rotational heat capacity cv/R of para (even J) or ortho (odd J) hydrogen, from its levels
/// E_J = J(J+1) k theta: the variance of E/kT over the populations.
fn h2_rotation(t: f64, ortho: bool) -> f64 {
    let (mut z, mut m1, mut m2) = (0.0, 0.0, 0.0);
    let mut j = if ortho { 1 } else { 0 };
    while j < 40 {
        let jf = j as f64;
        let e = (jf * (jf + 1.0) * THETA_ROT_H2) / t;
        let w = (2.0 * jf + 1.0) * exp(-e);
        z += w;
        m1 += w * e;
        m2 += w * e * e;
        j += 2;
    }
    let m1 = m1 / z;
    m2 / z - m1 * m1
}

/// Ideal-gas cp/R of one species at T: frozen part, vibrations (Einstein), and hydrogen's
/// rotation as normal hydrogen (ortho:para frozen at 3:1, as a gas cooled without a catalyst).
pub fn cp_over_r(s: &SpeciesData, t: f64) -> f64 {
    let mut cp = s.frozen_cp;
    for &(theta, g) in s.modes {
        let x = theta / t;
        if x < 500.0 {
            let ex = exp(x);
            cp += (g * x * x * ex) / ((ex - 1.0) * (ex - 1.0));
        }
    }
    if s.species == Species::H2 {
        cp += 0.25 * h2_rotation(t, false) + 0.75 * h2_rotation(t, true);
    }
    cp
}

/// Viscosity of one species at T, Pa s (Sutherland).
pub fn viscosity(s: &SpeciesData, t: f64) -> f64 {
    let (mu0, t0, c) = s.sutherland;
    let r = t / t0;
    mu0 * r * sqrt(r) * (t0 + c) / (t + c)
}

/// Saturation vapour pressure of water over liquid water at T (Pa), Buck (1996).
pub fn water_saturation(t: f64) -> f64 {
    let c = t - 273.15;
    611.21 * exp((18.678 - c / 234.5) * (c / (257.14 + c)))
}

/// The state of an ideal gas mix.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasState {
    /// kg/mol
    pub molar_mass: f64,
    pub gamma: f64,
    /// m/s
    pub c: f64,
    /// kg/m³
    pub rho: f64,
    /// Pa s
    pub mu: f64,
    pub temperature: f64,
    pub pressure: f64,
}

impl GasState {
    /// characteristic impedance (rayl, Pa s/m)
    pub fn z(&self) -> f64 {
        self.rho * self.c
    }
    /// kinematic viscosity (m²/s)
    pub fn nu(&self) -> f64 {
        self.mu / self.rho
    }
}

/// The ideal-gas state of `mix` at temperature `t` (K) and pressure `p` (Pa).
pub fn state(mix: &Mix, t: f64, p: f64) -> GasState {
    let total: f64 = mix.iter().map(|x| js_max(0.0, *x)).sum();
    let x = |i: usize| if total > 0.0 { js_max(0.0, mix[i]) / total } else { 0.0 };
    let (mut m, mut cp) = (0.0, 0.0);
    for i in 0..N_SPECIES {
        m += x(i) * SPECIES[i].molar_mass;
        cp += x(i) * cp_over_r(&SPECIES[i], t);
    }
    let gamma = cp / (cp - 1.0);
    // Wilke's mixing rule
    let mut mus = [0.0; N_SPECIES];
    for i in 0..N_SPECIES {
        mus[i] = viscosity(&SPECIES[i], t);
    }
    let mut mu = 0.0;
    for i in 0..N_SPECIES {
        if x(i) == 0.0 {
            continue;
        }
        let mut den = 0.0;
        for j in 0..N_SPECIES {
            if x(j) == 0.0 {
                continue;
            }
            let (mi, mj) = (SPECIES[i].molar_mass, SPECIES[j].molar_mass);
            let a = 1.0 + sqrt(mus[i] / mus[j]) * sqrt(sqrt(mj / mi));
            den += (x(j) * a * a) / sqrt(8.0 * (1.0 + mi / mj));
        }
        mu += (x(i) * mus[i]) / den;
    }
    GasState { molar_mass: m, gamma, c: sqrt((gamma * R * t) / m), rho: (p * m) / (R * t), mu, temperature: t, pressure: p }
}

/// Add water vapour at relative humidity `rh` (0..1) of the saturation pressure at `t`, diluting
/// the rest; the fraction stops at 1 (below the saturation pressure the water would boil).
pub fn humidify(mix: &Mix, t: f64, p: f64, rh: f64) -> (Mix, f64) {
    let total: f64 = mix.iter().map(|x| js_max(0.0, *x)).sum();
    let xw = js_min(1.0, js_max(0.0, rh) * water_saturation(t) / p);
    let mut out = [0.0; N_SPECIES];
    for i in 0..N_SPECIES {
        out[i] = if total > 0.0 { (1.0 - xw) * js_max(0.0, mix[i]) / total } else { 0.0 };
    }
    out[Species::H2O as usize] += xw;
    (out, xw)
}

/// Where a voice is and what it breathes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Conditions {
    /// the air outside (mole fractions)
    pub mix: Mix,
    /// Pa
    pub pressure: f64,
    /// the air's temperature outside (K)
    pub temperature: f64,
    /// the gas's temperature in the airway (K): the body's for a warm-blooded breather
    pub airway_temperature: f64,
    /// relative humidity the airway brings the gas to (0..1), at the airway temperature
    pub airway_humidity: f64,
}

/// The airway's gas as the engine sees it, with the state behind it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Airway {
    pub gas: Gas,
    /// the airway gas, SI
    pub state: GasState,
    /// the air outside, SI
    pub outside: GasState,
    /// mole fraction of water vapour in the airway
    pub water: f64,
    /// the airway is below the water's saturation pressure at its temperature (the Armstrong
    /// limit for a body at 37 °C, about 6.3 kPa): its water would boil
    pub boiling: bool,
}

impl Conditions {
    /// Sea-level air breathed by a human: the conditions the engine was tuned in.
    pub const EARTH: Conditions = Conditions { mix: DRY_AIR, pressure: ATM, temperature: 288.15, airway_temperature: BODY_TEMPERATURE, airway_humidity: 1.0 };

    pub fn airway(&self) -> Airway {
        let p = self.pressure;
        let (mix, water) = humidify(&self.mix, self.airway_temperature, p, self.airway_humidity);
        let st = state(&mix, self.airway_temperature, p);
        let rf = reference();
        let gas = Gas { rho: Gas::BODY_AIR.rho * (st.rho / rf.rho), c: Gas::BODY_AIR.c * (st.c / rf.c), mu: Gas::BODY_AIR.mu * (st.mu / rf.mu) };
        let boiling = self.airway_humidity > 0.0 && water_saturation(self.airway_temperature) >= p;
        Airway { gas, state: st, outside: state(&self.mix, self.temperature, p), water, boiling }
    }

    /// The gas for [`crate::engine::Engine::set_gas`].
    pub fn gas(&self) -> Gas {
        self.airway().gas
    }
}

/// The model's own warm, saturated sea-level air: the state `Gas::BODY_AIR` stands for.
pub fn reference() -> GasState {
    let c = Conditions::EARTH;
    let (mix, _) = humidify(&c.mix, c.airway_temperature, c.pressure, c.airway_humidity);
    state(&mix, c.airway_temperature, c.pressure)
}

/// Who is breathing: sets the airway's temperature and humidity for a preset.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Breather {
    /// a warm-blooded body (37 °C) that humidifies its breath, where the pressure lets it
    Human,
    /// a creature of the place: its airway at the outside temperature, dry
    Native,
}

/// A place, with its conditions and where they come from.
#[derive(Clone, Copy, Debug)]
pub struct Preset {
    pub id: &'static str,
    pub name: &'static str,
    pub mix: Mix,
    pub pressure: f64,
    pub temperature: f64,
    pub source: &'static str,
}

impl Preset {
    pub fn conditions(&self, who: Breather) -> Conditions {
        match who {
            Breather::Human => Conditions {
                mix: self.mix,
                pressure: self.pressure,
                temperature: self.temperature,
                airway_temperature: BODY_TEMPERATURE,
                // humid breath where water can stay liquid at 37 °C, dry below the Armstrong limit
                airway_humidity: if self.pressure > water_saturation(BODY_TEMPERATURE) { 1.0 } else { 0.0 },
            },
            Breather::Native => Conditions { mix: self.mix, pressure: self.pressure, temperature: self.temperature, airway_temperature: self.temperature, airway_humidity: 0.0 },
        }
    }
}

const fn mix(n2: f64, o2: f64, ar: f64, co2: f64, he: f64, h2: f64, ch4: f64, h2o: f64, sf6: f64) -> Mix {
    [n2, o2, ar, co2, he, h2, ch4, h2o, sf6]
}

/// Places to speak in. Surface (or 1 bar for the giants) conditions.
pub const PRESETS: &[Preset] = &[
    Preset { id: "earth", name: "Earth, sea level", mix: DRY_AIR, pressure: ATM, temperature: 288.15, source: "ISO 2533 standard atmosphere; dry-air composition (CO2 420 ppm)" },
    Preset { id: "mars_jezero", name: "Mars, Jezero crater", mix: mix(0.0259, 0.0, 0.0194, 0.951, 0.0, 0.0, 0.0, 0.0, 0.0), pressure: 600.0, temperature: 230.0, source: "Maurice et al. 2022, Nature 605 (Perseverance SuperCam); Franz et al. 2017 (SAM composition)" },
    Preset { id: "mars_hellas", name: "Mars, Hellas basin floor", mix: mix(0.0259, 0.0, 0.0194, 0.951, 0.0, 0.0, 0.0, 0.0, 0.0), pressure: 1340.0, temperature: 250.0, source: "Mars Climate Database; Franz et al. 2017 (composition)" },
    Preset { id: "mars_olympus", name: "Mars, Olympus Mons summit", mix: mix(0.0259, 0.0, 0.0194, 0.951, 0.0, 0.0, 0.0, 0.0, 0.0), pressure: 92.0, temperature: 170.0, source: "Mars Climate Database; Franz et al. 2017 (composition)" },
    Preset { id: "venus_surface", name: "Venus, surface", mix: mix(0.035, 0.0, 0.0, 0.965, 0.0, 0.0, 0.0, 0.0, 0.0), pressure: 9.2e6, temperature: 737.0, source: "Seiff et al. 1985, VIRA (Venus International Reference Atmosphere)" },
    Preset { id: "venus_55km", name: "Venus, 55 km cloud deck", mix: mix(0.035, 0.0, 0.0, 0.965, 0.0, 0.0, 0.0, 0.0, 0.0), pressure: 49600.0, temperature: 299.4, source: "Seiff et al. 1985, VIRA" },
    Preset { id: "titan", name: "Titan, surface", mix: mix(0.9435, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0565, 0.0, 0.0), pressure: 146700.0, temperature: 93.65, source: "Fulchignoni et al. 2005 (Huygens HASI); Niemann et al. 2010, JGR 115 (Huygens GCMS, 5.65 % methane at the surface)" },
    Preset { id: "saturn", name: "Saturn, 1 bar level", mix: mix(0.0, 0.0, 0.0, 0.0, 0.11, 0.8855, 0.0045, 0.0, 0.0), pressure: 1.0e5, temperature: 134.0, source: "Lindal et al. 1985 (Voyager radio occultation); Koskinen and Guerlet 2018 (helium)" },
    Preset { id: "neptune", name: "Neptune, 1 bar level", mix: mix(0.0, 0.0, 0.0, 0.0, 0.188, 0.797, 0.015, 0.0, 0.0), pressure: 1.0e5, temperature: 72.0, source: "Lindal 1992 (Voyager 2 radio occultation)" },
    Preset { id: "pluto", name: "Pluto, surface", mix: mix(0.995, 0.0, 0.0, 0.0, 0.0, 0.0, 0.005, 0.0, 0.0), pressure: 1.15, temperature: 37.0, source: "Gladstone et al. 2016, Science 351 (New Horizons)" },
    Preset { id: "heliox_dive", name: "Earth, heliox at 300 m (saturation dive)", mix: mix(0.0, 0.013, 0.0, 0.0, 0.987, 0.0, 0.0, 0.0, 0.0), pressure: 3.1e6, temperature: 303.15, source: "saturation diving practice: O2 held near 0.4 bar partial pressure in helium" },
    Preset { id: "helium_breath", name: "Earth, a breath of balloon helium", mix: mix(0.156, 0.042, 0.002, 0.0, 0.8, 0.0, 0.0, 0.0, 0.0), pressure: ATM, temperature: 293.15, source: "MODELLED: about 80 % helium left in the airway after one breath" },
    Preset { id: "sf6_breath", name: "Earth, a breath of SF6", mix: mix(0.156, 0.042, 0.002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.8), pressure: ATM, temperature: 293.15, source: "MODELLED: about 80 % SF6 left in the airway after one breath" },
];

pub fn preset(id: &str) -> Option<&'static Preset> {
    PRESETS.iter().find(|p| p.id == id)
}
