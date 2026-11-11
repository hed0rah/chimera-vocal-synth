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
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
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

    /// Every property finite and positive. An empty gas mix has no state (0/0 in the equation of
    /// state), and a sound speed of 0 sizes a tube with infinitely many sections.
    pub fn is_physical(&self) -> bool {
        [self.rho, self.c, self.mu].iter().all(|x| x.is_finite() && *x > 0.0)
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

/// Highest Courant number an airway is built with. The leapfrog is stable to 1; tubes are cut at
/// 0.93 (`human_airway.rs`, `creatures/build.rs`), so this is a margin over the cut. A tube with a
/// fixed section count can exceed it: in air, the child body's piriform fossa (three sections) is
/// at 1.11 at 8, 16 and 32 kHz (and the voice ran at the limiter ceiling there); in a fast gas
/// (helium, hydrogen) every short tube does. [`stable_oversample`] raises the internal rate until
/// the airway is under it.
pub const MAX_CFL: f64 = 0.95;

/// Highest internal rate the stability search goes to, Hz: a bound on the loop, not a target.
/// The fastest gas the atmosphere panel can make, pure hydrogen at 800 K (about 6 times air's
/// sound speed), needs about 450 kHz for the child body at an 8 kHz output.
pub const MAX_SIM_RATE: f64 = 1.0e6;

/// The smallest oversampling factor from `m` up at which the airway `build(m)` makes is under
/// [`MAX_CFL`], and that airway. `cfl` reads an airway's Courant number (None: nothing to check).
/// A tube with a fixed section count has a Courant number in proportion to 1/m, so the first step
/// jumps straight to the factor that satisfies it; the loop then checks.
pub fn stable_oversample<T>(fs_out: f64, m: usize, build: impl Fn(usize) -> T, cfl: impl Fn(&T) -> Option<f64>) -> (usize, T) {
    let mut m = m.max(1);
    let mut x = build(m);
    loop {
        let Some(c) = cfl(&x) else { break };
        if c <= MAX_CFL || fs_out * (m + 1) as f64 > MAX_SIM_RATE {
            break;
        }
        let jump = ((m as f64) * c / MAX_CFL).ceil() as usize;
        m = jump.max(m + 1).min(((MAX_SIM_RATE / fs_out).floor() as usize).max(m + 1));
        x = build(m);
    }
    (m, x)
}

/// Critical Reynolds number for onset of turbulence noise.
pub const RE_CRIT: f64 = 1800.0;
