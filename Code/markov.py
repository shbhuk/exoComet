#!/usr/bin/env python3
"""
markov.py — 1D Markov-chain model of cometary orbital evolution.

Simulates the orbital evolution of a population of comets originating
in a Kuiper belt, perturbed by a chain of planets. Comets start just
beyond the outermost planet and random-walk inward/outward in a
logarithmic distance coordinate until they are absorbed by:

  * colliding with the target planet,
  * colliding with some other (small or giant) planet,
  * being dynamically ejected from the system, or
  * reaching the outer edge of the grid (effectively unbound).

The model also returns the mean cometary lifetime.

Can be used either as a library (call `run(...)` with parameters
directly) or as a command-line tool (`python markov.py markov.in`).

Python 3.9+ port of the original Fortran `markov.f90`.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.constants import G, GM_sun, M_earth, au
from scipy.linalg import solve_banded

# --------------------------------------------------------------------
# Physical constants (cgs), sourced from astropy rather than hard-coded
# --------------------------------------------------------------------
GRAV_CONST = G.cgs.value          # cm^3 g^-1 s^-2
M_SUN = (GM_sun / G).cgs.value    # g  (derived the same way the original did)
M_EARTH = M_earth.cgs.value       # g
AU = au.cgs.value                 # cm
YEAR = u.yr.to(u.s)               # s  (Julian year)

# --------------------------------------------------------------------
# Model-specific constants
# --------------------------------------------------------------------
ECC = 0.5            # typical comet eccentricity when apoapse is inside the system
A_BOUND = 5e4 * AU    # distance of the outermost bound orbit
N_RH = 5.0            # target planet controls dynamics within N_RH Hill radii


# ======================================================================
# Linear algebra
# ======================================================================
def tridag(sub: np.ndarray, diag: np.ndarray, sup: np.ndarray,
           rhs: np.ndarray) -> np.ndarray:
    """
    Solve the tridiagonal system T u = rhs for u.

    `sub`, `diag`, `sup` and `rhs` are all arrays of length n, using the
    Numerical-Recipes convention: `sub[0]` and `sup[-1]` are unused
    (there is no sub-diagonal entry for row 0, nor a super-diagonal
    entry for the last row).

    Thin, vectorised wrapper around `scipy.linalg.solve_banded`
    (LAPACK ``?gtsv``) rather than a hand-rolled Thomas-algorithm loop.
    """
    n = len(diag)
    ab = np.zeros((3, n))
    ab[0, 1:] = sup[:-1]   # super-diagonal
    ab[1, :] = diag        # main diagonal
    ab[2, :-1] = sub[1:]   # sub-diagonal
    return solve_banded((1, 1), ab, rhs)


# ======================================================================
# Physics
# ======================================================================
def collision_ejection_probabilities(a, mplan, mstar, aout, density):
    """
    Per-state probability of a comet colliding with, or being ejected
    by, its controlling planet.

    Returns (pcol, pejc, prip), where `prip` = max(pcol, pejc) is the
    combined "removed from the chain" probability (collision and
    ejection cross sections overlap).
    """
    a_eff = np.minimum(a, aout)
    r_plan = (3 * mplan / (4 * np.pi * density)) ** (1 / 3)
    mscl = mplan / mstar

    pcol = (r_plan ** 2 / (a_eff ** 2 * ECC)
            * (1 + 2 * mscl * a_eff / (r_plan * ECC ** 2)))
    pejc = 4 * mscl ** 2 * a ** 2 / (a_eff ** 2 * ECC ** 3)

    pcol = np.minimum(pcol, 1.0)
    pejc = np.minimum(pejc, 1.0)

    prip = np.maximum(pcol, pejc)
    pejc = prip - pcol

    return pcol, pejc, prip


@dataclass
class RadialGrid:
    """The non-uniform radial grid of Markov-chain states."""
    n: int
    k: int
    a: np.ndarray
    y: np.ndarray
    torb: np.ndarray
    mplan: np.ndarray
    atc_lo: float
    atc_hi: float
    atran: float
    mout: float
    itc_lo: int
    itc_hi: int
    itran: int


def build_radial_grid(mstar, ain, atarg, atran, aout, mtarg, msmall,
                       mtran, mout, opt_target, opt_system,
                       verbose=True) -> RadialGrid:
    """
    Builds the logarithmic radial grid of states, and computes the
    orbital period and controlling-planet mass at each state.

    All arguments are in cgs. `atran` and `mout` may be adjusted
    internally (mirroring the original Fortran `intent(inout)`
    behaviour) and the adjusted values are returned via the result.
    """
    # Target planet controls the region within N_RH Hill radii
    h = (mtarg / 3 / mstar) ** (1 / 3)
    atc_lo = max(atarg * (1 - N_RH * h), ain)
    atc_hi = atarg * (1 + N_RH * h)
    atran = max(atran, atc_hi)

    if opt_target == 1:
        atc_hi = atarg

    if opt_system == 1:                 # no giant planets
        atran = aout
        expo = 0.0
        mout = msmall
    else:
        expo = np.log(mtran / mout) / np.log(atran / aout)

    if opt_system == 2:                 # no small planets
        atran = atc_hi

    # Y-coordinate of a few important states
    if expo == 0.0:
        ytran = np.log(atran / aout)
    else:
        ytran = (1 - (aout / atran) ** expo) / expo

    ytc_hi = ytran + mout / msmall * np.log(atc_hi / atran)
    ytarg = ytc_hi + mout / mtarg * np.log(atarg / atc_hi)
    ytc_lo = ytarg + mout / mtarg * np.log(atc_lo / atarg)
    yin = ytc_lo + mout / msmall * np.log(ain / atc_lo)
    ybound = 1 - aout / A_BOUND

    dy = 10 * mout / mstar          # width of each non-absorbing state in Y
    n = int((ybound - yin) / dy)    # number of non-absorbing states
    k = int(-yin / dy)              # starting state for comets

    # Y values and distances of each non-absorbing state
    y = yin + dy * np.arange(1, n + 1)

    if expo == 0.0:
        a_giant = aout * np.exp(y)
    else:
        a_giant = aout / (1 - expo * y) ** (1 / expo)

    a = np.select(
        condlist=[y >= 0, y < ytc_lo, y < ytc_hi, y < ytran],
        choicelist=[
            aout / (1 - y),
            atc_lo * np.exp(msmall * (y - ytc_lo) / mout),
            atc_hi * np.exp(mtarg * (y - ytc_hi) / mout),
            atran * np.exp(msmall * (y - ytran) / mout),
        ],
        default=a_giant,
    )

    # a(i) is monotonically increasing with i, so the smallest index
    # whose distance exceeds a threshold is just the first True entry.
    def first_index_above(threshold):
        return int(np.argmax(a > threshold)) + 1   # 1-based, for display

    itran = first_index_above(atran)
    itc_hi = first_index_above(atc_hi)
    itc_lo = first_index_above(atc_lo)

    if verbose:
        print(f' i, ytran:  {itran}   {ytran}')
        print(f' i, ytc_hi  {itc_hi}   {ytc_hi}')
        print(f' i, ytc_lo  {itc_lo}   {ytc_lo}')
        print(f' i, yin     1   {yin}')

    tn = 2 * np.pi * np.sqrt(aout ** 3 / GRAV_CONST / mstar)  # period of outermost planet
    torb = tn * (a / aout) ** 1.5

    mplan = np.select(
        condlist=[a > aout, (atc_lo < a) & (a < atc_hi), a < atran],
        choicelist=[mout, mtarg, msmall],
        default=mout * (a / aout) ** expo,
    )

    return RadialGrid(n=n, k=k, a=a, y=y, torb=torb, mplan=mplan,
                       atc_lo=atc_lo, atc_hi=atc_hi, atran=atran,
                       mout=mout, itc_lo=itc_lo, itc_hi=itc_hi,
                       itran=itran)


# ======================================================================
# Input parsing
# ======================================================================
def read_input(path) -> dict:
    """
    Parse a `markov.in`-style input file.

    Each line supplies one value, in the order:
        mstar, ain, atarg, atran, aout,
        mtarg, msmall, mtran, mout,
        density, opt_target, opt_system
    Anything after the first token on a line (Fortran-style '!'
    comments, or stray extra tokens) is ignored, mirroring Fortran's
    list-directed READ.
    """
    keys = ['mstar', 'ain', 'atarg', 'atran', 'aout',
            'mtarg', 'msmall', 'mtran', 'mout',
            'density', 'opt_target', 'opt_system']

    values = []
    with open(path) as f:
        for line in f:
            token = line.split('!', 1)[0].split()
            if not token:
                continue
            values.append(token[0].replace('d', 'e').replace('D', 'e'))
            if len(values) == len(keys):
                break

    if len(values) != len(keys):
        raise ValueError(f'{path}: expected {len(keys)} values, found {len(values)}')

    params = dict(zip(keys, (float(v) for v in values)))
    params['opt_target'] = int(params['opt_target'])
    params['opt_system'] = int(params['opt_system'])
    return params


# ======================================================================
# Top-level driver
# ======================================================================
@dataclass
class MarkovResult:
    """Outcome probabilities, mean lifetime, and the underlying grid."""
    grid: RadialGrid
    pcol: np.ndarray
    pejc: np.ndarray
    prip: np.ndarray
    ftarg: float
    fcol1: float
    fcol2: float
    fejc1: float
    fejc2: float
    fejc3: float
    fout: float
    tbar: float
    params_cgs: dict = field(repr=False)

    @property
    def total_probability(self) -> float:
        return (self.ftarg + self.fcol1 + self.fcol2
                + self.fejc1 + self.fejc2 + self.fejc3 + self.fout)


def run(mstar, ain, atarg, atran, aout, mtarg, msmall, mtran, mout,
        density, opt_target, opt_system,
        write_files: bool = True, outdir='.',
        verbose: bool = True) -> MarkovResult:
    """
    Run the full Markov-chain calculation.

    Parameters are given in the same "natural" units as `markov.in`:
    masses in solar/Earth masses, distances in AU, density in g/cm^3.
    Unit conversion to cgs is handled internally.

    Set `write_files=False` to skip writing `zones.out`/`log.out`
    (e.g. when calling `run()` programmatically many times).
    """
    outdir = Path(outdir)

    # --- convert to cgs -------------------------------------------------
    mstar_cgs = mstar * M_SUN
    mtarg_cgs = mtarg * M_EARTH
    msmall_cgs = msmall * M_EARTH
    mtran_cgs = mtran * M_EARTH
    mout_cgs = mout * M_EARTH
    ain_cgs = ain * AU
    atarg_cgs = atarg * AU
    atran_cgs = atran * AU
    aout_cgs = aout * AU

    grid = build_radial_grid(
        mstar_cgs, ain_cgs, atarg_cgs, atran_cgs, aout_cgs,
        mtarg_cgs, msmall_cgs, mtran_cgs, mout_cgs,
        opt_target, opt_system, verbose=verbose)

    pcol, pejc, prip = collision_ejection_probabilities(
        grid.a, grid.mplan, mstar_cgs, aout_cgs, density)

    if write_files:
        write_zones_file(outdir / 'zones.out', grid, pcol, pejc)

    n, k, a = grid.n, grid.k, grid.a
    atc_lo, atc_hi, atran_cgs = grid.atc_lo, grid.atc_hi, grid.atran
    ki = k - 1   # 0-based index of the comets' starting state

    # --- tridiagonal system describing one Markov step ------------------
    sub = -0.5 * (1 - prip)
    diag = np.ones(n)
    sup = -0.5 * (1 - prip)
    diag[0] = 1.0 if opt_target == 1 else 1 - 0.5 * (1 - prip[0])

    def absorbed_fraction(rhs: np.ndarray) -> float:
        return tridag(sub, diag, sup, rhs)[ki]

    # fraction reaching infinite distance
    rhs = np.zeros(n)
    rhs[-1] = 0.5 * (1 - prip[-1])
    fout = absorbed_fraction(rhs)

    # fraction hitting the target
    rhs = np.zeros(n)
    mask = (a < atc_hi) & (a > atc_lo)
    rhs[mask] = pcol[mask]
    if opt_target == 1:
        rhs[0] = 0.5 * (1 - prip[0])
    ftarg = absorbed_fraction(rhs)

    # fraction hitting another small planet
    rhs = np.zeros(n)
    if opt_system == 1:
        mask = (a >= atc_hi) | (a < atc_lo)
    else:
        mask = (a < atran_cgs) & ((a >= atc_hi) | (a < atc_lo))
    rhs[mask] = pcol[mask]
    fcol1 = absorbed_fraction(rhs)

    # fraction hitting another giant planet
    rhs = np.zeros(n)
    if opt_system >= 2:
        mask = a >= atran_cgs
        rhs[mask] = pcol[mask]
    fcol2 = absorbed_fraction(rhs)

    # fraction ejected by the target
    rhs = np.zeros(n)
    mask = (a < atc_hi) & (a > atc_lo)
    rhs[mask] = pejc[mask]
    fejc1 = absorbed_fraction(rhs)

    # fraction with a < a_out ejected by an intermediate planet
    rhs = np.zeros(n)
    mask = (a < aout_cgs) & ((a >= atc_hi) | (a < atc_lo))
    rhs[mask] = pejc[mask]
    fejc2 = absorbed_fraction(rhs)

    # fraction with a > a_out ejected by an intermediate planet
    rhs = np.zeros(n)
    mask = a >= aout_cgs
    rhs[mask] = pejc[mask]
    fejc3 = absorbed_fraction(rhs)

    # --- mean cometary lifetime, via the transposed system --------------
    sub_t = np.zeros(n)
    sub_t[1:] = -0.5 * (1 - prip[:-1])
    sup_t = np.zeros(n)
    sup_t[:-1] = -0.5 * (1 - prip[1:])
    diag_t = diag  # same main diagonal

    rhs = np.zeros(n)
    rhs[ki] = 1.0
    uvec = tridag(sub_t, diag_t, sup_t, rhs)
    tbar = float(np.sum(uvec * grid.torb))

    result = MarkovResult(
        grid=grid, pcol=pcol, pejc=pejc, prip=prip,
        ftarg=ftarg, fcol1=fcol1, fcol2=fcol2,
        fejc1=fejc1, fejc2=fejc2, fejc3=fejc3, fout=fout, tbar=tbar,
        params_cgs=dict(
            mstar=mstar_cgs, ain=ain_cgs, atarg=atarg_cgs,
            atran=atran_cgs, aout=aout_cgs, mtarg=mtarg_cgs,
            msmall=msmall_cgs, mtran=mtran_cgs, mout=grid.mout,
            density=density, opt_target=opt_target, opt_system=opt_system,
        ),
    )

    if write_files:
        write_log_file(outdir / 'log.out', result, verbose=verbose)

    return result


# ======================================================================
# Output
# ======================================================================
def write_zones_file(path, grid: RadialGrid, pcol, pejc) -> None:
    header = ('     i        y          a (AU)     P (y)    '
              'Mplan (M_E)     Pcol        Pejc   \n')
    with open(path, 'w') as f:
        f.write(header)
        f.write('-' * 83 + '\n')
        for i in range(grid.n):
            f.write(
                ' {:6d} {:10.5f} {:13.5f} {:11.4e} {:11.4e} {:11.4e} '
                '{:11.4e}\n'.format(
                    i + 1, grid.y[i], grid.a[i] / AU,
                    grid.torb[i] / YEAR, grid.mplan[i] / M_EARTH,
                    pcol[i], pejc[i]))


def format_report(result: MarkovResult) -> str:
    p = result.params_cgs
    grid = result.grid

    target_desc = ('Target: absorbing boundary' if p['opt_target'] == 1
                    else 'Target: ordinary planet   ')
    system_desc = {1: 'small planets only        ',
                    2: 'giant planets only        '}.get(
        p['opt_system'], 'giant and small planets   ')

    lines = [
        '',
        f" Stellar mass (solar):              {p['mstar'] / M_SUN:7.3f}",
        f" Innermost planet distance (AU):    {p['ain'] / AU:11.4e}",
        f" Distance of target (AU):           {p['atarg'] / AU:11.4e}",
        f" Inner edge of target control (AU): {grid.atc_lo / AU:11.4e}",
        f" Outer edge of target control (AU): {grid.atc_hi / AU:11.4e}",
        f" Small/giant planet transition(AU): {p['atran'] / AU:11.4e}",
        f" Outermost giant-plan dist (AU):    {p['aout'] / AU:11.4e}",
        '',
        f" Target planet mass (Earth):        {p['mtarg'] / M_EARTH:11.4e}",
        f" Small-planet mass (Earth):         {p['msmall'] / M_EARTH:11.4e}",
        f" Innermost giant-plan mass (Earth): {p['mtran'] / M_EARTH:11.4e}",
        f" Outermost giant-plan mass (Earth): {p['mout'] / M_EARTH:11.4e}",
        f" Planetary density (g/cm^3):        {p['density']:11.4e}",
        f' {target_desc}',
        f' System: {system_desc}',
        '',
        f' Number of radial zones:           {grid.n:6d}',
        f' Starting zone for comets:         {grid.k:6d}',
        f' Small/giant plan transition zone: {grid.itran:6d}',
        f' Inner edge of target control:     {grid.itc_lo:6d}',
        f' Outer edge of target control:     {grid.itc_hi:6d}',
        '',
        f' Comets that hit the target:        {result.ftarg:10.3e}',
        f' Hit another small planet:          {result.fcol1:10.3e}',
        f' Hit another giant plannet:         {result.fcol2:10.3e}',
        f' Ejected by the target:             {result.fejc1:10.3e}',
        f' Otherwise ejected when a < a_out:  {result.fejc2:10.3e}',
        f' Otherwise ejected when a > a_out:  {result.fejc3:10.3e}',
        f' Reached outer edge of the grid:    {result.fout:10.3e}',
        f' Total probability:                 {result.total_probability:10.3e}',
        '',
        f' Mean comet lifetime (year): {result.tbar / YEAR:10.3e}',
        '',
    ]
    return '\n'.join(lines) + '\n'


def write_log_file(path, result: MarkovResult, verbose: bool = True) -> None:
    text = format_report(result)
    if verbose:
        print(text, end='')
    with open(path, 'w') as f:
        f.write(text)


# ======================================================================
# Command-line interface
# ======================================================================
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description='1D Markov-chain model of cometary orbital evolution.')
    parser.add_argument('infile', nargs='?', default='markov.in',
                         help="input file (default: 'markov.in')")
    parser.add_argument('-o', '--outdir', default='.',
                         help="directory for zones.out/log.out (default: '.')")
    parser.add_argument('-q', '--quiet', action='store_true',
                         help='suppress console output')
    args = parser.parse_args(argv)

    params = read_input(args.infile)
    run(**params, outdir=args.outdir, verbose=not args.quiet)
    return 0


if __name__ == '__main__':
    sys.exit(main())
