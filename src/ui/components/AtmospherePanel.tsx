import { ATMOSPHERES } from "../../wasm/engine.ts";
import { LEVEL_CAP, useStore, type Breather } from "../store.ts";
import { Segmented, Slider, Toggle } from "./Controls.tsx";
import { Chip } from "./Screen.tsx";

const fmt = (x: number, digits = 3) => (Math.abs(x) >= 1000 ? Math.round(x).toLocaleString("en-US") : x.toPrecision(digits));
const pressure = (pa: number) => (pa >= 1e5 ? `${fmt(pa / 1e5)} bar` : pa >= 1000 ? `${fmt(pa / 1000)} kPa` : `${fmt(pa)} Pa`);
const kelvin = (k: number) => `${Math.round(k)} K (${Math.round(k - 273.15)} °C)`;
const db = (r: number) => `${r >= 1 ? "+" : ""}${(20 * Math.log10(r)).toFixed(1)} dB`;

/** The top bar's air button: where the voice is, and the panel with its conditions. */
export function AtmosphereButton() {
  const a = useStore((s) => s.atmosphere);
  const setOpen = useStore((s) => s.setAtmosphereOpen);
  const name = a.preset ? ATMOSPHERES.presets.find((p) => p.id === a.preset)?.name : "Custom air";
  const load = useStore((s) => s.loadAtmosphere);
  const breathe: [string, string][] = [["earth", "Air"], ["helium_breath", "Helium"], ["sf6_breath", "SF6"]];
  return (
    <div className="atmo">
      <div className="atmo-breathe" role="group" aria-label="Breathe">
        <span className="atmo-label">Breathe</span>
        {breathe.map(([id, label]) => (
          <button key={id} className={`atmo-quick${a.preset === id && a.breather === "human" ? " on" : ""}`} onClick={() => load(id, "human")}>
            {label}
          </button>
        ))}
      </div>
      <button className={`atmo-button${a.open ? " on" : ""}`} aria-expanded={a.open} onClick={() => setOpen(!a.open)} title="The air the voice breathes and speaks into">
        <span className="atmo-label">Air</span>
        <span className="atmo-name">{name}</span>
      </button>
      {a.open && <AtmospherePanel />}
    </div>
  );
}

function AtmospherePanel() {
  const a = useStore((s) => s.atmosphere);
  const load = useStore((s) => s.loadAtmosphere);
  const set = useStore((s) => s.setConditions);
  const setReal = useStore((s) => s.setRealLevel);
  const setOpen = useStore((s) => s.setAtmosphereOpen);
  const c = a.conditions;
  const r = a.report;
  const preset = a.preset ? ATMOSPHERES.presets.find((p) => p.id === a.preset) : null;
  const total = c.mix.reduce((x, y) => x + Math.max(0, y), 0) || 1;
  const l = a.loudness;
  const setMix = (i: number, v: number) => set({ mix: c.mix.map((x, j) => (j === i ? v : x)) });
  return (
    <div className="atmo-panel" data-layout="overlay" role="dialog" aria-label="Atmosphere">
      <div className="atmo-head">
        <h2>Atmosphere</h2>
        <Segmented<Breather>
          label="Who breathes it"
          value={a.breather}
          onChange={(b) => (a.preset ? load(a.preset, b) : useStore.setState((s) => ({ atmosphere: { ...s.atmosphere, breather: b } })))}
          options={[
            { value: "human", label: "Human body" },
            { value: "native", label: "Native" },
          ]}
        />
        <button className="atmo-close" onClick={() => setOpen(false)} aria-label="Close">
          ×
        </button>
      </div>
      <p className="atmo-hint">
        {a.breather === "human"
          ? "A body at 37 °C warms and humidifies what it breathes: the formants follow the airway's temperature, not the outside air's."
          : "A creature of the place: its airway is at the outside temperature, and dry."}
      </p>
      <div className="atmo-presets">
        {ATMOSPHERES.presets.map((p) => (
          <Chip key={p.id} on={a.preset === p.id} onClick={() => load(p.id)}>
            {p.name}
          </Chip>
        ))}
      </div>
      {preset && <p className="atmo-source">Source: {preset.source}</p>}

      <div className="atmo-cols">
        <div className="atmo-col">
          <h3>Conditions</h3>
          <Slider label="Pressure" value={c.pressure} min={1} max={1e7} log onChange={(v) => set({ pressure: v })} format={pressure} />
          <Slider label="Air temperature" value={c.temperature} min={30} max={800} onChange={(v) => set({ temperature: v })} format={kelvin} />
          <Slider label="Airway temperature" value={c.airwayTemperature} min={30} max={800} onChange={(v) => set({ airwayTemperature: v })} format={kelvin} />
          <Slider label="Airway humidity" value={c.airwayHumidity} min={0} max={1} onChange={(v) => set({ airwayHumidity: v })} format={(v) => `${Math.round(v * 100)} %`} />
          <h3>Gas mix</h3>
          {ATMOSPHERES.species.map((s, i) => (
            <Slider key={s.formula} label={`${s.formula} ${s.name}`} value={c.mix[i]} min={0} max={1} onChange={(v) => setMix(i, v)} format={() => `${((100 * Math.max(0, c.mix[i])) / total).toFixed(c.mix[i] / total < 0.01 && c.mix[i] > 0 ? 2 : 1)} %`} />
          ))}
        </div>
        <div className="atmo-col">
          <h3>In the airway</h3>
          <dl className="atmo-read">
            <dt>speed of sound</dt>
            <dd>{fmt(r.airway.c)} m/s</dd>
            <dt>resonances move</dt>
            <dd>
              ×{r.cRatio.toFixed(3)} ({r.semitones >= 0 ? "+" : ""}
              {r.semitones.toFixed(1)} semitones)
            </dd>
            <dt>density</dt>
            <dd>{fmt(r.airway.rho)} kg/m³</dd>
            <dt>viscosity, kinematic</dt>
            <dd>{fmt(r.airway.nu * 1e6)} mm²/s</dd>
            <dt>gamma, molar mass</dt>
            <dd>
              {r.airway.gamma.toFixed(3)}, {fmt(r.airway.molarMass * 1000)} g/mol
            </dd>
            <dt>water vapour</dt>
            <dd>{(r.water * 100).toFixed(1)} %</dd>
            <dt>compressibility Z</dt>
            <dd>{r.airway.compressibility.toFixed(3)}{r.airway.compressibility < 0.9 ? " (a dense fluid)" : ""}</dd>
          </dl>
          {r.nearCritical && <p className="atmo-warn">Near a critical point: the airway gas is between a gas and a liquid, and the real-gas model (Peng-Robinson) is good to about 10 to 15 % here.</p>}
          {r.boiling && <p className="atmo-warn">Below the water's saturation pressure at this airway temperature (the Armstrong limit at 37 °C is 6.3 kPa): a body's water would boil.</p>}
          <h3>Outside</h3>
          <dl className="atmo-read">
            <dt>speed of sound</dt>
            <dd>{fmt(r.outside.c)} m/s</dd>
            <dt>density</dt>
            <dd>{fmt(r.outside.rho)} kg/m³</dd>
            <dt>impedance</dt>
            <dd>{fmt(r.outside.z)} rayl</dd>
            <dt>voice's output vs Earth</dt>
            <dd>{l ? db(l.output) : "measuring"}</dd>
            <dt>at a listener vs Earth</dt>
            <dd>{l ? db(l.real) : "measuring"}</dd>
          </dl>
          <Toggle label={`Real loudness${l && l.real > LEVEL_CAP ? " (capped at +12 dB)" : ""}`} checked={a.realLevel} onChange={setReal} />
          <p className="atmo-hint">
            Off, every place plays as loud as Earth and only the sound changes. On, you hear the pressure a listener would: the voice's airflow, measured on a reference vowel, times the density of the air it travels in. Real gas (Peng-Robinson): within 1 % of reference values at every preset's outside conditions.
          </p>
        </div>
      </div>
    </div>
  );
}
