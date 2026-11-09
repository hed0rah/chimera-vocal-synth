//! Physical constants for airway acoustics, CGS units throughout
//! (cm, g, s; pressure in dyn/cm² = 0.1 Pa; volume velocity in cm³/s).

/// Air density at body temperature, saturated (g/cm³).
pub const RHO: f64 = 1.14e-3;
/// Speed of sound in warm humid air (cm/s).
pub const C_SOUND: f64 = 35000.0;
/// Dynamic viscosity of air (poise = g/(cm·s)).
pub const MU: f64 = 1.86e-4;
/// Kinematic viscosity (cm²/s).
pub const NU: f64 = MU / RHO;
/// ρc² — bulk modulus of air (dyn/cm²).
pub const RHO_C2: f64 = RHO * C_SOUND * C_SOUND;
/// 1 cm H2O in dyn/cm².
pub const CM_H2O: f64 = 980.665;

/// The gas filling an airway: everything the acoustics needs from it. An airway takes its gas
/// when it is built (there is no global gas); [`Gas::BODY_AIR`] is the constants above, so an
/// airway built in it renders exactly as before. Other gases come from [`crate::atmosphere`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gas {
    /// density (g/cm³)
    pub rho: f64,
    /// speed of sound (cm/s)
    pub c: f64,
    /// dynamic viscosity (poise)
    pub mu: f64,
}

impl Gas {
    /// Warm, saturated air in a human airway: the constants the engine was tuned in.
    pub const BODY_AIR: Gas = Gas { rho: RHO, c: C_SOUND, mu: MU };

    /// kinematic viscosity (cm²/s)
    #[inline]
    pub fn nu(&self) -> f64 {
        self.mu / self.rho
    }

    /// ρc², the bulk modulus (dyn/cm²), multiplied in the same order as [`RHO_C2`]
    #[inline]
    pub fn rho_c2(&self) -> f64 {
        self.rho * self.c * self.c
    }

    /// characteristic impedance ρc (g/(cm²·s), 10 rayl)
    #[inline]
    pub fn z(&self) -> f64 {
        self.rho * self.c
    }
}

impl Default for Gas {
    fn default() -> Gas {
        Gas::BODY_AIR
    }
}

/// Soft-tissue wall parameters per unit area (Ishizaka, French & Flanagan 1975, cheek/neck).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallParams {
    /// mass per unit area (g/cm²)
    pub mass: f64,
    /// resistance per unit area (dyn·s/cm³)
    pub resistance: f64,
    /// stiffness per unit area (dyn/cm³)
    pub stiffness: f64,
}

pub const WALL_SOFT: WallParams = WallParams { mass: 1.5, resistance: 1600.0, stiffness: 3.0e5 };
/// Cartilaginous trachea: stiffer and lighter.
pub const WALL_TRACHEA: WallParams = WallParams { mass: 0.6, resistance: 1000.0, stiffness: 1.0e6 };
/// Nasal mucosa: lossy.
pub const WALL_NASAL: WallParams = WallParams { mass: 1.0, resistance: 3000.0, stiffness: 1.0e6 };

/// Critical Reynolds number for onset of turbulence noise.
pub const RE_CRIT: f64 = 1800.0;
