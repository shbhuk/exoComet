# exoComet: A 1D Markov-Chain Model of Cometary Orbital Evolution

This repository contains the code used to compute the analytic
Markov-chain estimates of cometary dynamical outcomes described in the paper `Comet Delivery in Compact Planetary Systems Is Stellar-Mass Independent` by John Chambers & Shubham Kanodia, which is currently under prep.


Given a star with a chain of planets, the model follows a population
of comets that start their dynamical lives just beyond the outermost
planet (e.g., scattered out of a Kuiper-belt-like reservoir) and
random-walk inward and outward in semi-major axis under repeated
planetary perturbations. Each state in the chain has some probability
of ejecting or colliding with its controlling planet. The full set of
absorption probabilities is obtained by solving a linear system for a
tridiagonal Markov-chain transition matrix, rather than by direct
N-body integration or Monte Carlo sampling.


## Requirements

- Python ≥ 3.9
- `numpy`
- `scipy`
- `astropy` (used for physical constants: `G`, `GM_sun`, `M_earth`,
  `au`, and the Julian year)

The code has been tested with numpy==1.21.0, scipy==1.13.1, astropy==6.0.1

```bash
pip install numpy scipy astropy
```

## Usage

### Command line

```bash
python markov.py ../examples/markov.in -o ../examples
```

This reads `markov.in`, runs the model, prints a summary to the
console, and writes two files to the working directory (or to
`--outdir`, if given):

- **`zones.out`** — one row per radial-grid state, giving its `y`
  coordinate, semi-major axis (AU), orbital period (yr), controlling
  planet mass (Earth masses), and per-state collision/ejection
  probabilities.
- **`log.out`** — the same human-readable summary printed to the
  console (input echo + outcome probabilities + mean lifetime).

```
usage: markov.py [-h] [-o OUTDIR] [-q] [infile]

positional arguments:
  infile                input file (default: 'markov.in')

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                         directory for zones.out/log.out (default: '.')
  -q, --quiet            suppress console output
```

### As a library

`run()` can be called directly with all parameters, without needing
an input file on disk — convenient for parameter sweeps, notebooks,
or embedding in other pipelines:

```python
import markov

result = markov.run(
    mstar=0.11,        # solar masses
    ain=0.052,         # AU, innermost planet distance
    atarg=0.052,       # AU, target planet distance
    atran=3.0,         # AU, small/giant planet transition
    aout=1.04,         # AU, outermost planet distance
    mtarg=1.0,         # Earth masses
    msmall=1.0,        # Earth masses
    mtran=1.0,         # Earth masses
    mout=1.0,          # Earth masses
    density=3.0,       # g/cm^3
    opt_target=0,      # 0 = ordinary planet, 1 = absorbing boundary
    opt_system=1,      # 1 = small planets only, 2 = giant only, 3 = both
    write_files=False, # skip writing zones.out/log.out
    verbose=False,
)

print(result.ftarg, result.fcol1, result.fout, result.tbar / markov.YEAR)
print(result.total_probability)  # should be ~1.0
```

`run()` returns a `MarkovResult` dataclass with the outcome
probabilities (`ftarg`, `fcol1`, `fcol2`, `fejc1`, `fejc2`, `fejc3`,
`fout`), the mean lifetime `tbar` (seconds), a `total_probability`
convenience property, and the underlying `RadialGrid` (state
positions, orbital periods, controlling-planet masses) plus the
per-state `pcol`/`pejc`/`prip` arrays, for further analysis or
plotting.

To parse an existing `.in` file into keyword arguments for `run()`:

```python
params = markov.read_input('markov.in')
result = markov.run(**params)
```

## Input file format (`markov.in`)

We include a sample input file in the `examples` directory. One value per line, in order, with everything after the first token
(e.g. `!` comments) ignored:

```
mstar        stellar mass (solar masses)
ain          distance of innermost planet (AU)
atarg        distance of target planet (AU)
atran        small-planet / giant-planet transition (AU)
aout         distance of outermost planet (AU)
mtarg        mass of target planet (Earth masses)
msmall       mass of each small planet (Earth masses)
mtran        mass of innermost giant planet (Earth masses)
mout         mass of outermost giant planet (Earth masses)
density      planetary bulk density (g/cm^3)
opt_target   0 = ordinary planet, 1 = absorbing boundary
opt_system   1 = small planets only, 2 = giant planets only, 3 = both
```

## Validation

This python implementation is converted from the original Fortran
(`markov.f90`) reference implementation written by John Chambers and also included with this repository. For the same input, the two
codes agree to numerical precision.

## Citation

If you use this code, please cite the accompanying paper 

FINDME

## License

Released under the [MIT License](LICENSE).
