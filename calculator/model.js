/* Tropical stress budget: browser (and Node) version of stressbudget.py.
 * Same equations as the Python model and ESI Section S2:
 *   solar position (NOAA/Meeus) + refraction, Erbs decomposition, isotropic sky,
 *   Faiman module temperature, constant-vapour-pressure surface RH,
 *   Peck (damp heat), Arrhenius x (G/1000)^0.72 (light + heat), Arrhenius (dark heat),
 *   Coffin-Manson (thermal cycling).
 */
(function (root) {
  'use strict';
  var KB = 8.617e-5;
  var MOUNTS = { open: [25.0, 6.84], roof: [20.0, 0.0] };
  var WIND_SCALE = 0.72;
  var REQUIRED = ['T2M', 'RH2M', 'ALLSKY_SFC_SW_DWN', 'WS10M'];
  var D2R = Math.PI / 180, R2D = 180 / Math.PI;

  function psat(T) { return 6.112 * Math.exp(17.62 * T / (243.12 + T)); }

  /* ---------- NASA POWER CSV parser ---------- */
  function parsePower(text) {
    var lines = text.split(/\r?\n/);
    var end = -1;
    for (var i = 0; i < lines.length; i++) { if (lines[i].trim() === '-END HEADER-') { end = i; break; } }
    if (end < 0) throw new Error('Not a NASA POWER CSV file (no "-END HEADER-" line found).');
    var header = lines.slice(0, end).join('\n');
    var m = header.match(/Latitude\s+(-?[\d.]+)\s+Longitude\s+(-?[\d.]+)/);
    var lat = m ? parseFloat(m[1]) : null, lon = m ? parseFloat(m[2]) : null;
    var s = header.match(/in (LST|UTC)/);
    var std = s ? s[1] : null;
    var cols = lines[end + 1].split(',').map(function (c) { return c.trim(); });
    var idx = {};
    cols.forEach(function (c, k) { idx[c] = k; });
    ['YEAR', 'MO', 'DY', 'HR'].concat(REQUIRED).forEach(function (c) {
      if (!(c in idx)) throw new Error('Missing column ' + c + '. The file needs T2M, RH2M, ALLSKY_SFC_SW_DWN and WS10M.');
    });
    var data = { YEAR: [], MO: [], DY: [], HR: [] };
    REQUIRED.forEach(function (c) { data[c] = []; });
    for (var j = end + 2; j < lines.length; j++) {
      var L = lines[j];
      if (!L.trim()) continue;
      var f = L.split(',');
      ['YEAR', 'MO', 'DY', 'HR'].forEach(function (c) { data[c].push(parseInt(f[idx[c]], 10)); });
      REQUIRED.forEach(function (c) {
        var v = parseFloat(f[idx[c]]);
        data[c].push(v === -999 || isNaN(v) ? NaN : v);
      });
    }
    return { lat: lat, lon: lon, std: std, data: data, n: data.YEAR.length };
  }

  function interpolate(a) {          // linear, edges filled with nearest value
    var n = a.length, missing = 0, last = -1;
    for (var i = 0; i < n; i++) {
      if (isNaN(a[i])) { missing++; continue; }
      if (last < 0) { for (var k = 0; k < i; k++) a[k] = a[i]; }
      else if (i - last > 1) { for (var k2 = last + 1; k2 < i; k2++) a[k2] = a[last] + (a[i] - a[last]) * (k2 - last) / (i - last); }
      last = i;
    }
    if (last < 0) throw new Error('A required parameter has no valid values.');
    for (var k3 = last + 1; k3 < n; k3++) a[k3] = a[last];
    return missing;
  }

  /* ---------- solar position (NOAA / Meeus) ---------- */
  function solarPosition(jd, lat, lon) {
    var T = (jd - 2451545.0) / 36525.0;
    var L0 = (280.46646 + T * (36000.76983 + T * 0.0003032)) % 360;
    var M = 357.52911 + T * (35999.05029 - 0.0001537 * T);
    var e = 0.016708634 - T * (0.000042037 + 0.0000001267 * T);
    var Mr = M * D2R;
    var C = Math.sin(Mr) * (1.914602 - T * (0.004817 + 0.000014 * T)) + Math.sin(2 * Mr) * (0.019993 - 0.000101 * T) + Math.sin(3 * Mr) * 0.000289;
    var trueLong = L0 + C;
    var omega = 125.04 - 1934.136 * T;
    var lambda = trueLong - 0.00569 - 0.00478 * Math.sin(omega * D2R);
    var eps0 = 23 + (26 + (21.448 - T * (46.815 + T * (0.00059 - T * 0.001813))) / 60) / 60;
    var eps = eps0 + 0.00256 * Math.cos(omega * D2R);
    var decl = Math.asin(Math.sin(eps * D2R) * Math.sin(lambda * D2R));
    var y = Math.pow(Math.tan(eps * D2R / 2), 2);
    var L0r = L0 * D2R;
    var eqTime = 4 * R2D * (y * Math.sin(2 * L0r) - 2 * e * Math.sin(Mr) + 4 * e * y * Math.sin(Mr) * Math.cos(2 * L0r)
      - 0.5 * y * y * Math.sin(4 * L0r) - 1.25 * e * e * Math.sin(2 * Mr));       // minutes
    var utcMin = ((jd + 0.5) % 1) * 1440;
    var tst = utcMin + eqTime + 4 * lon;                                            // true solar time, minutes
    var ha = tst / 4 - 180;
    ha = ((ha + 180) % 360 + 360) % 360 - 180;
    var latr = lat * D2R, har = ha * D2R;
    var cosz = Math.sin(latr) * Math.sin(decl) + Math.cos(latr) * Math.cos(decl) * Math.cos(har);
    cosz = Math.max(-1, Math.min(1, cosz));
    var zen = Math.acos(cosz) * R2D;
    var az = Math.atan2(Math.sin(har), Math.cos(har) * Math.sin(latr) - Math.tan(decl) * Math.cos(latr)) * R2D + 180; // from north, eastward
    var elev = 90 - zen, refr = 0;
    if (elev >= -(0.26667 + 0.5667)) {      // SPA refraction, 1010 mbar / 12 degC
      refr = (1013.25 / 1010) * (283 / (273 + 12)) * 1.02 / (60 * Math.tan((elev + 10.3 / (elev + 5.11)) * D2R));
    }
    return { zenith: zen, apparentZenith: zen - refr, azimuth: (az + 360) % 360 };
  }

  function julianDay(ms) { return ms / 86400000 + 2440587.5; }

  function dayOfYear(ms) {
    var d = new Date(ms);
    var start = Date.UTC(d.getUTCFullYear(), 0, 1);
    return Math.floor((ms - start) / 86400000) + 1;
  }

  /* ---------- module conditions ---------- */
  function prepare(rec, opts) {
    opts = opts || {};
    var lat = opts.lat != null ? opts.lat : rec.lat, lon = opts.lon != null ? opts.lon : rec.lon;
    if (lat == null || lon == null || isNaN(lat) || isNaN(lon)) throw new Error('Latitude and longitude are missing from the file; enter them manually.');
    var issues = [];
    var std = rec.std;
    if (!std) { issues.push('time standard not stated (LST assumed)'); std = 'LST'; }
    var d = rec.data, n = rec.n;
    REQUIRED.forEach(function (c) { var k = interpolate(d[c]); if (k) issues.push(c + ': ' + k + ' missing hours (interpolated)'); });
    var shiftMs = lon / 15 * 3600000;
    var utc = new Float64Array(n), dayKey = new Int32Array(n);
    for (var i = 0; i < n; i++) {
      var t0 = Date.UTC(d.YEAR[i], d.MO[i] - 1, d.DY[i], d.HR[i], 30);
      var u = std === 'UTC' ? t0 : t0 - shiftMs;
      utc[i] = u;
      dayKey[i] = Math.floor((u + shiftMs) / 86400000);          // local-solar-time calendar day
    }
    var tilt = opts.tilt != null ? opts.tilt : Math.max(10, Math.round(Math.abs(lat)));
    var surfAz = lat > 0 ? 180 : 0, albedo = 0.25;
    var ct = Math.cos(tilt * D2R), st = Math.sin(tilt * D2R);
    var G = new Float64Array(n);
    for (var j = 0; j < n; j++) {
      var sp = solarPosition(julianDay(utc[j]), lat, lon);
      var ghi = Math.max(d.ALLSKY_SFC_SW_DWN[j], 0);
      var B = 2 * Math.PI * (dayOfYear(utc[j]) - 1) / 365;
      var dniExtra = 1366.1 * (1.00011 + 0.034221 * Math.cos(B) + 0.00128 * Math.sin(B) + 0.000719 * Math.cos(2 * B) + 7.7e-5 * Math.sin(2 * B));
      var cz = Math.cos(sp.zenith * D2R);
      var kt = ghi / (dniExtra * Math.max(cz, 0.065));
      kt = Math.min(Math.max(kt, 0), 1);
      var df = kt <= 0.22 ? 1 - 0.09 * kt : kt <= 0.8 ? 0.9511 - 0.1604 * kt + 4.388 * kt * kt - 16.638 * Math.pow(kt, 3) + 12.336 * Math.pow(kt, 4) : 0.165;
      var dhi = df * ghi, dni = (ghi - dhi) / cz;
      if (sp.zenith > 87 || ghi < 0 || dni < 0 || !isFinite(dni)) { dni = 0; dhi = ghi; }
      var za = sp.apparentZenith * D2R;
      var proj = ct * Math.cos(za) + st * Math.sin(za) * Math.cos((sp.azimuth - surfAz) * D2R);
      proj = Math.max(-1, Math.min(1, proj));
      var poa = Math.max(dni * proj, 0) + dhi * (1 + ct) * 0.5 + ghi * albedo * (1 - ct) * 0.5;
      G[j] = isFinite(poa) ? Math.min(Math.max(poa, 0), 1400) : 0;
    }
    return { G: G, dayKey: dayKey, Ta: d.T2M, RH: d.RH2M, WS: d.WS10M, n: n, years: n / 8766, lat: lat, lon: lon, tilt: tilt, std: std, issues: issues };
  }

  function moduleConditions(p, mount) {
    var U = MOUNTS[mount], n = p.n;
    var Tm = new Float64Array(n), RHm = new Float64Array(n);
    for (var i = 0; i < n; i++) {
      Tm[i] = p.Ta[i] + p.G[i] / (U[0] + U[1] * p.WS[i] * WIND_SCALE);
      RHm[i] = Math.min(100, Math.max(0, p.RH[i] * psat(p.Ta[i]) / psat(Tm[i])));
    }
    return { G: p.G, Tm: Tm, RH: RHm, dayKey: p.dayKey, years: p.years, n: n };
  }

  /* Laboratory hours (cycles for TC) equivalent to one field year */
  function labHoursPerYear(m, test, o) {
    o = o || {};
    var T = o.T != null ? o.T : 85, RHt = o.RH != null ? o.RH : 85, Ea = o.Ea != null ? o.Ea : 0.6;
    var nexp = o.n != null ? o.n : 2.66, order = o.order != null ? o.order : 0.72, dT = o.dT != null ? o.dT : 125;
    var s = 0, i, inv = 1 / (T + 273.15);
    if (test === 'TC') {
      var cur = null, mx = -Infinity, mn = Infinity;
      for (i = 0; i <= m.n; i++) {
        if (i === m.n || m.dayKey[i] !== cur) {
          if (cur !== null) s += Math.pow((mx - mn) / dT, 2);
          if (i === m.n) break;
          cur = m.dayKey[i]; mx = -Infinity; mn = Infinity;
        }
        if (m.Tm[i] > mx) mx = m.Tm[i];
        if (m.Tm[i] < mn) mn = m.Tm[i];
      }
      return s / m.years;
    }
    for (i = 0; i < m.n; i++) {
      var a = Math.exp(-Ea / KB * (1 / (m.Tm[i] + 273.15) - inv));
      if (test === 'DH') s += Math.pow(m.RH[i] / RHt, nexp) * a;
      else if (test === 'L') s += Math.pow(m.G[i] / 1000, order) * a;
      else if (test === 'D') s += a;
      else throw new Error('Unknown test ' + test);
    }
    return s / m.years;
  }

  function fieldYears(m, test, duration, o) { return duration / labHoursPerYear(m, test, o); }

  var api = { parsePower: parsePower, prepare: prepare, moduleConditions: moduleConditions,
    labHoursPerYear: labHoursPerYear, fieldYears: fieldYears, solarPosition: solarPosition, psat: psat };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.StressBudget = api;
})(this);
