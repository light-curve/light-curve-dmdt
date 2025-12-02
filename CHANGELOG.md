# Changelog

All notable changes to `light-curve-dmdt-exec` binary crates will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Bump minimum supported Rust version to 1.85
- Update Rust edition to 2024
- Update `clap` to 4.5 (removed pinned version constraint)
- Update dependencies to latest compatible versions

## [0.6.2] 2024-12-04

- Explicitly set the minimum supported Rust version to 1.67

## [0.6.1] 2022-12-30

### Changed

- Bump `light-curve-dmdt` 0.6.0 -> 0.7.0

## [0.6.0] 2022-11-30

### Added

- The first stand-alone release of the `dmdt` binary produced by `light-curve-dmdt-exec` crate, see
  `light-curve-dmdt/Changelog.md` for previous releases as a part of `light-curve-dmdt` crate

### Changed

- **braking** The crate was split into two: the executable was moved into a separate crate `light-curve-dmdt-exec`
- **braking** Update to `clap` v4 caused removing of shot arguments `-w`, `-h` and reformatted help message
- The project repository was split from other 'light-curve*' crates and moved
  into <https://gituhb.com/light-curve/light-curve-dmdt>
