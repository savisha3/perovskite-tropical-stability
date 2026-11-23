// Compare the browser model (calculator/model.js) with Table 1 of the review.
// Usage (from the repository root): node tests/check_js_model.js
const fs = require('fs'), path = require('path');
const SB = require('../calculator/model.js');
const root = path.join(__dirname, '..');
const cities = fs.readFileSync(path.join(root, 'data/cities.csv'), 'utf8').trim().split('\n').slice(1).map(l => l.split(','));
const ref = {}; const rl = fs.readFileSync(path.join(root, 'tests/table1_reference.csv'), 'utf8').trim().split('\n');
const rh = rl[0].split(','); rl.slice(1).forEach(l => { const f = l.split(','); ref[f[0]] = {}; rh.forEach((h, k) => { if (k) ref[f[0]][h] = +f[k]; }); });
let worst = 0, bad = 0;
for (const c of cities) {
  const [city, , file, , , tilt] = c;
  const rec = SB.parsePower(fs.readFileSync(path.join(root, file), 'utf8'));
  const p = SB.prepare(rec, { tilt: +tilt });
  const o = SB.moduleConditions(p, 'open'), r = SB.moduleConditions(p, 'roof');
  const got = {
    'damp_heat_open_Ea0.6': SB.fieldYears(o, 'DH', 1000, { Ea: 0.6 }),
    'damp_heat_open_Ea0.3': SB.fieldYears(o, 'DH', 1000, { Ea: 0.3 }),
    'damp_heat_open_Ea0.95': SB.fieldYears(o, 'DH', 1000, { Ea: 0.95 }),
    'light_heat_open_Ea0.61': SB.fieldYears(o, 'L', 1000, { T: 85, Ea: 0.61 }),
    'light_heat_open_Ea0.43': SB.fieldYears(o, 'L', 1000, { T: 85, Ea: 0.43 }),
    'dark_heat_open': SB.fieldYears(o, 'D', 1000, { T: 85, Ea: 0.6 }),
    'thermal_cycling_open': SB.fieldYears(o, 'TC', 200),
    'thermal_cycling_roof': SB.fieldYears(r, 'TC', 200),
  };
  const line = [];
  for (const k in got) {
    const rel = Math.abs(got[k] - ref[city][k]) / ref[city][k];
    worst = Math.max(worst, rel);
    if (Math.abs(got[k] - ref[city][k]) > 0.051) { bad++; line.push(`${k}: ${ref[city][k]} vs ${got[k].toFixed(3)}`); }
  }
  console.log(city.padEnd(17), line.length ? 'DIFF ' + line.join('; ') : 'ok', p.issues.join('; '));
}
console.log(`Largest relative difference from Table 1: ${(worst * 100).toFixed(2)}%`);
console.log(bad ? `${bad} values differ by more than 0.05 yr` : 'ALL WITHIN 0.05 yr (Table 1 is reported to 0.1 yr)');
