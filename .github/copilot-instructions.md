# Copilot Instructions for light-curve-dmdt

## Project Overview

This repository contains a Rust library and executable for transforming astronomical light curves into dm-dt (magnitude difference vs time difference) space. The implementation is based on papers by Mahabal et al. (2011, 2017) and Soraisam et al. (2020).

## Repository Structure

This is a Cargo workspace with two crates:

- `light-curve-dmdt` (library crate in `light-curve-dmdt/`): Core library for dm-dt map generation
- `light-curve-dmdt-exec` (binary crate, root `Cargo.toml`): Command-line executable `dmdt`

## Rust Version Requirements

- Minimum Supported Rust Version (MSRV): 1.67
- Edition: 2018

## Building and Testing

```bash
# Build all targets with all features
cargo build --all-targets --all-features --workspace

# Run all tests
cargo test --all-features --workspace

# Check formatting
cargo fmt -- --check

# Run clippy lints
cargo clippy --all-targets --all-features --workspace -- -D warnings
```

## Cargo Features (light-curve-dmdt crate)

- `doc-images`: Adds example image to HTML docs (for docs.rs)
- `png`: Adds `to_png()` function to save dm-dt map as PNG
- `serde`: Enables serde serialization for `DmDt` struct
- `full`: Enables all features
- `default`: No features enabled by default

## Code Style

- Follow Rust 2018 edition idioms
- Use `cargo fmt` for formatting
- All clippy warnings should be addressed
- Pre-commit hooks are configured for fmt, check, and clippy

## Key Types

- `DmDt`: Main struct for dm-dt map configuration and computation
- `Grid` / `GridTrait`: Grid abstraction for dt and dm axes
- `Erf` trait: Error function implementations (exact and approximate)

## Dependencies

Key dependencies include:
- `ndarray`: N-dimensional arrays (publicly re-exported)
- `png`: PNG encoding (optional, re-exported when enabled)
- `clap`: CLI argument parsing (for executable)

## Changelog

Changes should be documented in:
- `CHANGELOG.md` (root): For `light-curve-dmdt-exec` binary
- `light-curve-dmdt/CHANGELOG.md`: For `light-curve-dmdt` library

Follow [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format.
