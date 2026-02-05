# CLAUDE.md - AI Assistant Guide for light-curve-dmdt

This document provides comprehensive guidance for AI assistants working with the `light-curve-dmdt` codebase.

## Project Overview

This repository contains a Rust library and executable for transforming astronomical light curves into dm-dt (magnitude difference vs time difference) space. The implementation is based on peer-reviewed papers:
- [Mahabal et al. 2011](https://ui.adsabs.harvard.edu/abs/2011BASI...39..387M)
- [Mahabal et al. 2017](https://arxiv.org/abs/1709.06257)
- [Soraisam et al. 2020](https://doi.org/10.3847/1538-4357/ab7b61)

**Related Projects:**
- Python bindings: https://github.com/light-curve/light-curve-python
- Documentation: https://docs.rs/light-curve-dmdt

## Repository Structure

```
light-curve-dmdt/
├── Cargo.toml                    # Workspace root + binary crate config
├── Cargo.lock                    # Lock file for reproducible builds
├── README.md                     # Main documentation (auto-generated)
├── CHANGELOG.md                  # Changelog for executable (light-curve-dmdt-exec)
├── example.png                   # Example output image (auto-generated)
│
├── .ci/                          # CI/CD scripts
│   ├── compose-readme.py         # Regenerates README.md and example.png
│   └── readme-example.sh         # Example command for documentation
│
├── .github/
│   ├── workflows/test.yml        # GitHub Actions CI workflow
│   ├── copilot-instructions.md   # GitHub Copilot instructions
│   └── dependabot.yml            # Dependency update configuration
│
├── light-curve-dmdt/             # Library crate
│   ├── Cargo.toml
│   ├── CHANGELOG.md              # Changelog for library
│   ├── README.md                 # Symlinked to root README
│   ├── src/
│   │   ├── lib.rs                # Module re-exports and crate docs
│   │   ├── dmdt.rs               # Core DmDt struct (main API)
│   │   ├── grid.rs               # Grid abstractions (ArrayGrid, LinearGrid, LgGrid)
│   │   ├── erf.rs                # Error function implementations
│   │   ├── float_trait.rs        # Float trait definition
│   │   ├── images.rs             # PNG export (optional)
│   │   └── util.rs               # Utility functions
│   └── benches/                  # Criterion benchmarks
│       ├── lib.rs
│       ├── cond_prob.rs
│       ├── erf.rs
│       ├── grid.rs
│       └── gausses.rs
│
└── light-curve-dmdt-exec/
    └── main.rs                   # CLI executable source
```

## Quick Reference Commands

```bash
# Build all targets with all features
cargo build --all-targets --all-features --workspace

# Run all tests
cargo test --all-features --workspace

# Check formatting
cargo fmt -- --check

# Run clippy lints (treat warnings as errors)
cargo clippy --all-targets --all-features --workspace -- -D warnings

# Run benchmarks
cargo bench --package light-curve-dmdt

# Regenerate README and example.png
./.ci/compose-readme.py
```

## Key Technical Details

### Rust Version Requirements
- **Edition:** 2024
- **MSRV (Minimum Supported Rust Version):** 1.85

### Cargo Features (library crate)
| Feature | Description |
|---------|-------------|
| `doc-images` | Embeds example image in HTML docs (for docs.rs) |
| `png` | Adds `to_png()` function for PNG export |
| `serde` | Enables serialization/deserialization for `DmDt` |
| `full` | Enables all features |
| `default` | No features enabled by default |

### Key Dependencies
- **ndarray 0.17** - N-dimensional arrays (publicly re-exported)
- **libm** - Error function implementations
- **enum_dispatch** - Zero-cost trait dispatch
- **clap 4** - CLI argument parsing (executable only)
- **png** - PNG encoding (optional)
- **serde** - Serialization (optional)

## Core Types and Architecture

### DmDt (dmdt.rs)
The primary struct for dm-dt map operations.

```rust
pub struct DmDt<T> {
    pub dt_grid: Grid<T>,  // Time difference grid
    pub dm_grid: Grid<T>,  // Magnitude difference grid
}
```

**Key methods:**
- `from_lgdt_dm_limits()` - Create with logarithmic dt and linear dm grids
- `from_grids()` - Create with custom grids
- `points()` - Map (t, m) pairs to unity values in 2D grid
- `gausses::<Erf>()` - Apply Gaussian smearing with error propagation
- `cond_prob::<Erf>()` - Calculate conditional probability
- `dt_points()` - Count dt values per cell

### Grid System (grid.rs)
Abstraction for axis grids with enum-dispatch optimization.

```rust
#[enum_dispatch(GridTrait<T>)]
pub enum Grid<T> {
    Array(ArrayGrid<T>),   // Custom borders, O(log n) lookup
    Linear(LinearGrid<T>), // Uniform spacing, O(1) lookup
    Lg(LgGrid<T>),         // Logarithmic spacing, O(1) lookup
}
```

### Error Functions (erf.rs)
Two implementations for Gaussian smearing:

- **`ExactErf`** - Uses `libm` for exact erf() computation
- **`Eps1Over1e3Erf`** - Approximate erf with ~1e-3 max error (faster, uses interpolation)

### Float Trait (float_trait.rs)
Composite trait for numeric operations, implemented for `f32` and `f64`.

## Code Style Guidelines

1. **Follow Rust 2024 edition idioms**
2. **Use `cargo fmt` for formatting** - CI enforces this
3. **Address all clippy warnings** - CI runs with `-D warnings`
4. **Pre-commit hooks are configured** for fmt, check, and clippy

### API Design Patterns
- Use `ArrayView` for return types (not `&Array`)
- Use `ArrayRef` for function parameters accepting arrays
- Generic over float type `T` where appropriate
- Use `thiserror` for custom error types

## Testing

### Running Tests
```bash
cargo test --all-features --workspace
```

### Test Locations
- `dmdt.rs` - Tests for dm-dt operations consistency
- `erf.rs` - Tests for error function approximation accuracy

### Assertions
- Use `approx` crate for floating-point comparisons
- Use `static_assertions` for compile-time trait checks

## CI/CD Pipeline

The GitHub Actions workflow (`.github/workflows/test.yml`) runs:

| Job | Description |
|-----|-------------|
| `test` | Runs tests on ubuntu-latest and ARM |
| `fmt` | Checks code formatting |
| `clippy` | Runs lints with `-D warnings` |
| `readme` | Verifies README.md and example.png are up-to-date |
| `build` | Standard build validation |
| `msrv-build` | Builds against minimum supported Rust version |

## Changelog Management

Two separate changelogs following [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format:

- **`CHANGELOG.md`** (root) - For `light-curve-dmdt-exec` binary
- **`light-curve-dmdt/CHANGELOG.md`** - For `light-curve-dmdt` library

When making changes, update the appropriate changelog under the `## [Unreleased]` section.

## Regenerating README and Example

The README.md and example.png are auto-generated. When modifying:
- The example in `.ci/readme-example.sh`
- The `--help` output of the CLI

Run:
```bash
./.ci/compose-readme.py
```

This script:
1. Regenerates `example.png` from the shell script
2. Updates `README.md` with current `--help` output
3. Updates `README.md` with the example shell script

**CI verifies these are up-to-date** by running the script and checking for uncommitted changes.

## Common Tasks

### Adding a New Feature to the Library
1. Implement in appropriate module under `light-curve-dmdt/src/`
2. Export in `lib.rs` if public
3. Add tests in the same file or a test module
4. Update `light-curve-dmdt/CHANGELOG.md`
5. Run full CI checks locally before committing

### Modifying the CLI
1. Edit `light-curve-dmdt-exec/main.rs`
2. Regenerate README: `./.ci/compose-readme.py`
3. Update root `CHANGELOG.md`

### Version Bumps
- Library: Update `light-curve-dmdt/Cargo.toml` version
- Executable: Update root `Cargo.toml` version
- Update respective CHANGELOG.md files

## Pre-commit Hooks

Configured in `.pre-commit-config.yaml`:
- Trailing whitespace trimming
- EOF fixing
- YAML/TOML validation
- `cargo fmt`
- `cargo check`
- `cargo clippy`
- README composition

Install with:
```bash
pre-commit install
```

## Performance Considerations

- The library uses `enum_dispatch` for zero-cost grid abstraction
- `Eps1Over1e3Erf` provides ~10x speedup over `ExactErf` with acceptable accuracy
- Release profile uses LTO and single codegen unit for optimization
- Benchmarks available in `light-curve-dmdt/benches/`

## Important Notes for AI Assistants

1. **Always run the full CI checks** before suggesting changes are complete:
   ```bash
   cargo fmt -- --check && cargo clippy --all-targets --all-features --workspace -- -D warnings && cargo test --all-features --workspace
   ```

2. **Regenerate README after CLI changes** - CI will fail otherwise

3. **Use appropriate changelog** - library vs executable have separate changelogs

4. **Respect MSRV** - Don't use Rust features newer than 1.85

5. **ndarray is re-exported** - Users should use `light_curve_dmdt::ndarray` rather than adding ndarray as a direct dependency

6. **Feature gates matter** - PNG and serde are optional; ensure code compiles without them
