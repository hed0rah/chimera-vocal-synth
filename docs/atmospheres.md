# Atmospheres

The voice can breathe something other than warm air: a mix of gases at a pressure and a temperature, in a
body that warms and humidifies it (a human) or one at the outside temperature (a creature of the place).
Presets put it on Mars, Venus, Titan, the giants, Pluto, in a heliox dive or after a breath of helium or SF6.

## What the gas changes

An airway's acoustics need three numbers from its gas, carried by `physics::Gas`:

| quantity | where it acts | what you hear |
|---|---|---|
| speed of sound c | node compliance (V / rho c^2), radiation load, the section lengths | resonances move by c / c_air (helium speech) |
| density rho | branch inertance, kinetic (Bernoulli) losses, radiation, lung port | how much air the folds move, how the tract loads them |
| viscosity mu | boundary-layer losses, the Reynolds number of turbulence | formant bandwidths, when a fricative starts to hiss |

An airway takes its gas when it is built (`AirwayNetwork::new_in`, `HumanAirway::new_in`,
`HumanVoice::new_in`, `CreatureVoice::new_in`, the creatures' `create(fs, gas)`); there is no global gas.
`Gas::BODY_AIR` is the engine's old constants, multiplied in the same order, so a voice built in it renders
bit-identically and no golden test moves. `Engine::set_gas` rebuilds the voice as `set_body` does.

**The trap: stability.** Section lengths come from the leapfrog's Courant limit, c T / dx <= 0.93, but each
tube keeps a minimum number of sections. In a fast gas (hydrogen at 37 °C is about 1050 m/s, heliox about
990) those sections are shorter than the limit allows at the usual rate, and the voice blows up (a human on
Saturn: CFL 1.62 at 96 kHz). `human_voice::airway_for` and `CreatureVoice::new_in` raise the oversampling
until the airway's Courant number is under `physics::MAX_CFL` (0.95). Air's worst case is 0.91 (the child
body), so air never raises it.

## From conditions to a gas (`atmosphere.rs`)

The controls are the physical conditions: mole fractions of N2, O2, Ar, CO2, He, H2, CH4, H2O and SF6, the
pressure, the outside temperature, and the airway's temperature and relative humidity. For an ideal gas:

- c = sqrt(gamma R T / M), with gamma from each species' cp(T): translation and rotation, harmonic
  vibrational modes (Einstein), and hydrogen's quantised rotation as normal hydrogen (ortho:para 3:1);
- rho = P M / (R T);
- mu from Sutherland's law per species, mixed by Wilke's rule;
- water vapour from Buck's saturation pressure; below it the airway's water would boil (the Armstrong
  limit, 6.3 kPa at 37 °C), so a human preset below it breathes dry.

A condition is applied as a ratio to the same model's warm saturated sea-level air (`reference()`), so the
model's small errors cancel and `Conditions::EARTH` gives `Gas::BODY_AIR` exactly.

## How good it is (`crates/chimera/tests/atmosphere.rs`)

- Handbook values at 20 °C: air 343.2 m/s, 1.204 kg/m³, 18.1 µPa s; helium 1007 m/s; SF6 130 to 140 m/s.
- Against real-fluid values (CoolProp) at nine sites: c within 2.5 %, rho within 4 % (Venus at 92 bar and
  cold Titan are the least ideal); mu within 10 %, 25 % for hydrogen below 150 K.
- Helium speech: in 80/20 heliox F3 rises with c (within 6 %) and F1 rises least, as measured in divers
  (Fant and Lindqvist 1968): the walls give the closed tract a resonance of its own.
- Every preset renders a finite, audible human voice, and every creature a finite call, in every preset.

## Loudness

The engine's output is the radiated flow's derivative. A listener hears pressure, rho_outside dU/dt, so
`atmosphere::loudness` renders a reference vowel in the gas and in air and reports both: the output (what the
source and tract do) and the pressure at a listener. The page plays every place as loud as Earth by default;
"Real loudness" plays the pressure a listener would hear, capped at +12 dB.

| preset, human | output | at a listener |
|---|---|---|
| Mars, Jezero | +10.6 dB | -28.4 dB |
| Venus, surface | -47.6 dB | -13.0 dB |
| Titan | -3.5 dB | +9.0 dB |
| Saturn, 1 bar | +6.1 dB | -9.4 dB |
| Pluto | +10.2 dB | -71.2 dB |
| heliox at 300 m | -20.1 dB | -7.3 dB |

## Not modelled

- Fluid loading of the folds: the glottis is parametric, so its pitch is the performer's in every gas. Dense
  gas would lower a self-oscillating fold.
- Propagation between the mouth and the listener: absorption, dispersion and delay over distance.
- Real-gas behaviour at high pressure, beyond the ideal-gas ratios above.
