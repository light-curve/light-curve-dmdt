# dm–dt map plotter

[![docs.rs badge](https://docs.rs/light-curve-dmdt/badge.svg)](https://docs.rs/light-curve-dmdt)
![testing](https://github.com/light-curve/light-curve-dmdt/actions/workflows/test.yml/badge.svg)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/light-curve/light-curve-dmdt/master.svg)](https://results.pre-commit.ci/latest/github/light-curve/light-curve-dmdt/master)

Rust library (`light-curve-dmdt`) and executable (`light-curve-dmdt-exec`) to transform light curve into dm–dt space, the implementation is based on papers
[Mahabal et al. 2011](https://ui.adsabs.harvard.edu/abs/2011BASI...39..387M), [Mahabal et al. 2017](https://arxiv.org/abs/1709.06257), [Soraisam et al. 2020](https://doi.org/10.3847/1538-4357/ab7b61).

If you are looking for Python bindings for this package, please see <https://github.com/light-curve/light-curve-python>

## Executable

The executable `dmdt` can be installed by running `cargo install light-curve-dmdt-exec`. You need rust toolchain to be
installed in your system, consider using your OS package manager or [rustup](https://rustup.rs) utility.

Example of conditional probability dm–dt map plotting for linear dm grid `[-1.5; 1.5)` with 64 cells and logarithmic dt
grid `[1; 100)` with 96 cells:

```sh
curl https://ztf.snad.space/dr4/csv/633207400004730 | # Get some ZTF data
tail +2 | # chomp CSV header
awk -F, '{print $3"	"$4"	"$5}' | # print needed columns and change separator to tab
dmdt \
  --max-abs-dm=1.5 --height=64 \
  --min-lgdt=0 --max-lgdt=2 --width=96 \
  --smear --approx-smearing \
  --norm=lgdt --norm=max \
  --output=example.png
```

![Example dm-dt map][example_png]

[example_png]: example.png

### `dmdt --help`

<details><summary>expand</summary>

```text
Program for dm-dt maps generator from light curves
Usage: dmdt [OPTIONS] --min-lgdt <FLOAT> --max-lgdt <FLOAT> --max-abs-dm <FLOAT>
Options:
  -i, --input <FILE>
          Path of the input file, should be built of space-separated columns of time, magnitude and
          magnitude error (required for --smare only). If '-' is given (the default), then the input
          is taken from the stdin
          [default: -]
  -o, --output <FILE>
          Path of the output PNG file. If '-' is given (the default), then outputs to the stdout
          [default: -]
  -s, --smear
          Produce dm-``smeared'' output using observation errors, which must be the third column of
          the input. Instead of just adding some value to the lg(dt)-dm cell, the whole lg(dt) =
          const row is filled by normally distributed dm-probabilities
      --min-lgdt <FLOAT>
          Left border of the lg(dt) grid, note that decimal logarithm is required, i.e. -1.0 input
          means 0.1 time units
      --max-lgdt <FLOAT>
          Right border of the lg(dt) grid, note that decimal logarithm is required, i.e. 2.0 input
          means 100.0 time units
      --max-abs-dm <FLOAT>
          Maximum dm value, the considered dm interval would be [-max-abs-dm, +max-abs-dm)
      --width <INT>
          number of lg(dt) cells, width of the output image
          [default: 128]
      --height <INT>
          number of dm cells, height of the output image
          [default: 128]
      --approx-smearing
          speed up smearing using approximate error function
  -n, --norm <normalisation>
          Normalisation to do after dmdt map building. The order of operations is:1) build dmdt map,
          each dm-lgdt pair brings a unity value to dmdt space;2) if --norm=lgdt, then divide each
          cell value by the total number of the corresponding lgdt pairs, i.e. divide each cell of
          some column by the integral value in the column (including values out of the interval of
          [-max_abs_dm; max_abs_dm)); 3) if --norm=max, then divide each cell by the overall maximum
          value; 4) if any of --norm=lgdt or --norm=max is specified, then all values should be in
          [0; 1] interval, so they are multiplied by 255 and casted to uint8 to make it possible to
          save dmdt map as a PNG file.
          [possible values: lgdt, max]
  -h, --help
          Print help information (use `-h` for a summary)
  -V, --version
          Print version information
```

</details>

## Rust crate

```rust
use light_curve_dmdt::{DmDt, Eps1Over1e3Erf};
use ndarray::Array1;

let dmdt = DmDt::from_lgdt_dm_limits(0.0, 2.0, 96, 1.5, 64);

let t = Array1::linspace(0.0, 100.0, 101);
let m = t.mapv(|x| 2.0 * f64::sin(x));
let err2 = Array1::ones(t.len()) * 0.01;

let prob = dmdt.cond_prob::<Eps1Over1e3Erf>(t.as_slice().unwrap(), m.as_slice().unwrap(), err2.as_slice().unwrap());
```
