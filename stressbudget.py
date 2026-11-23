"""Tropical stress budget for perovskite solar cells.

Converts hourly NASA POWER weather data into the number of field years that
one accelerated laboratory test (damp heat, light and heat, dark heat or
thermal cycling) is equivalent to.

Model (see Section 3 and ESI Section S2 of the accompanying review):
  * plane-of-array irradiance: Erbs decomposition + isotropic sky (pvlib)
  * module temperature: Faiman model, open rack U0 = 25, U1 = 6.84;
    rooftop U0 = 20, U1 = 0; 10 m wind scaled to 2 m by 0.72
  * module-surface relative humidity: constant vapour pressure, Magnus (WMO)
  * damp heat: Peck model, RH exponent n = 2.66, Ea = 0.6 eV
  * light and heat: (G/1000)^0.72 x Arrhenius, Ea = 0.61 eV
  * dark heat: Arrhenius, Ea = 0.6 eV
  * thermal cycling: Coffin-Manson, exponent 2, reference swing 125 K

Usage:
  python stressbudget.py data/weather/kuala_lumpur.csv --tilt 10
"""
import argparse
import re

import numpy as np
import pandas as pd
import pvlib

KB = 8.617e-5                                   # Boltzmann constant, eV/K
REQUIRED = ['T2M', 'RH2M', 'ALLSKY_SFC_SW_DWN', 'WS10M']
MOUNTS = {'open': (25.0, 6.84), 'roof': (20.0, 0.0)}   # Faiman U0, U1
WIND_SCALE = 0.72                               # (2/10)^0.2: 10 m -> 2 m


def psat(T):
    """Saturation vapour pressure over water (hPa), Magnus form (WMO)."""
    return 6.112 * np.exp(17.62 * T / (243.12 + T))


def arrhenius(Ea, Tm, Ttest):
    """Arrhenius rate at module temperature Tm relative to test temperature Ttest (deg C)."""
    return np.exp(-Ea / KB * (1 / (Tm + 273.15) - 1 / (Ttest + 273.15)))


def load(path):
    """Read a NASA POWER hourly CSV (LST or UTC) and return a weather record."""
    lines = open(path).read().splitlines()
    end = lines.index('-END HEADER-')
    header = '\n'.join(lines[:end])
    m = re.search(r'Latitude\s+(-?[\d.]+)\s+Longitude\s+(-?[\d.]+)', header)
    lat, lon = map(float, m.groups())
    std = re.search(r'in (LST|UTC)', header)
    std = std.group(1) if std else None
    d = pd.read_csv(path, skiprows=end + 1).replace(-999, np.nan)
    issues = []
    if std is None:
        issues.append('time standard not stated (LST assumed)')
    missing = [c for c in REQUIRED if c not in d]
    if missing:
        raise ValueError(f'missing parameters {missing}')
    for c in REQUIRED:
        n = d[c].isna().sum()
        if n:
            issues.append(f'{c}: {n} missing hours (interpolated)')
            d[c] = d[c].interpolate(limit_direction='both')
    years = len(d) / 8766
    t0 = pd.to_datetime(dict(year=d.YEAR, month=d.MO, day=d.DY, hour=d.HR)) + pd.Timedelta(minutes=30)
    shift = pd.to_timedelta(lon / 15, unit='h')
    utc = pd.DatetimeIndex(t0 if std == 'UTC' else t0 - shift).tz_localize('UTC')
    lst = utc.tz_localize(None) + shift
    return dict(d=d, lat=lat, lon=lon, std=std, years=years, utc=utc, lst=lst, issues=issues)


def module(w, tilt, mount):
    """Hourly plane-of-array irradiance G, module temperature Tm and surface RH."""
    d, lat, lon, utc = w['d'], w['lat'], w['lon'], w['utc']
    sp = pvlib.solarposition.get_solarposition(utc, lat, lon)
    ghi = d.ALLSKY_SFC_SW_DWN.clip(lower=0).values
    dec = pvlib.irradiance.erbs(ghi, sp.zenith.values, utc)
    G = pvlib.irradiance.get_total_irradiance(
        tilt, 180 if lat > 0 else 0, sp.apparent_zenith.values, sp.azimuth.values,
        dec['dni'], ghi, dec['dhi'], model='isotropic')['poa_global']
    G = np.nan_to_num(np.asarray(G, dtype=float)).clip(0, 1400)
    U0, U1 = MOUNTS[mount]
    Ta = d.T2M.values
    Tm = Ta + G / (U0 + U1 * d.WS10M.values * WIND_SCALE)
    RHm = np.clip(d.RH2M.values * psat(Ta) / psat(Tm), 0, 100)
    return dict(G=G, Tm=Tm, RH=RHm, day=w['lst'].date, years=w['years'])


def lab_hours_per_year(m, test, T=85, RH=85, Ea=0.6, n=2.66, order=0.72, dT=125):
    """Laboratory hours (cycles for TC) equivalent to one field year."""
    if test == 'DH':
        return ((m['RH'] / RH) ** n * arrhenius(Ea, m['Tm'], T)).sum() / m['years']
    if test == 'L':
        return ((m['G'] / 1000) ** order * arrhenius(Ea, m['Tm'], T)).sum() / m['years']
    if test == 'D':
        return arrhenius(Ea, m['Tm'], T).sum() / m['years']
    if test == 'TC':
        swing = pd.Series(m['Tm']).groupby(m['day']).agg(lambda s: s.max() - s.min())
        return ((swing / dT) ** 2).sum() / m['years']
    raise ValueError(test)


def field_years(m, test, duration, **kw):
    """Field years equivalent to a test of the given duration (h, or cycles for TC)."""
    return duration / lab_hours_per_year(m, test, **kw)


def summary(path, tilt=None):
    """Field years for the standard tests (as in Table 1 of the review)."""
    w = load(path)
    tilt = tilt if tilt is not None else max(10, round(abs(w['lat'])))
    out = {'lat': w['lat'], 'lon': w['lon'], 'time_standard': w['std'], 'years': round(w['years'], 2),
           'tilt': tilt, 'issues': '; '.join(w['issues']) or 'none'}
    for mt in MOUNTS:
        m = module(w, tilt, mt)
        for Ea in (0.6, 0.3, 0.95):
            out[f'{mt} damp heat 1000 h, Ea {Ea}'] = field_years(m, 'DH', 1000, Ea=Ea)
        for Ea in (0.61, 0.43):
            out[f'{mt} light+heat 85C 1000 h, Ea {Ea}'] = field_years(m, 'L', 1000, T=85, Ea=Ea)
        out[f'{mt} dark heat 85C 1000 h'] = field_years(m, 'D', 1000, T=85, Ea=0.6)
        out[f'{mt} thermal cycling 200 cycles'] = field_years(m, 'TC', 200)
    return out


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('file', help='NASA POWER hourly CSV with T2M, RH2M, ALLSKY_SFC_SW_DWN and WS10M')
    p.add_argument('--tilt', type=float, default=None, help='module tilt (deg); default: rounded latitude, minimum 10')
    a = p.parse_args()
    for k, v in summary(a.file, a.tilt).items():
        print(f'{k:40s} {v:.2f}' if isinstance(v, float) else f'{k:40s} {v}')
