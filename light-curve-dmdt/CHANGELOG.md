# Changelog

All notable changes to `light-curve-dmdt` library crate will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

--

### Changed

- **breaking** `to_png()` function and `png` crate re-export are now under non-default `png` Cargo feature
- **Braking** `Grid` trait is renamed to `GridTrait`, while `Grid` enum is introduced to handle all grid variants https://github.com/light-curve/light-curve-dmdt/pull/17

### Deprecated

--

### Removed

--

### Fixed

--

### Security

--

## [0.6.1] 2022-12-28

### Fixed

- [docs.rs](https://docs.rs/light-curve-dmdt) build

## [0.6.0] 2022-11-30

### Added

- `light-curve-dmdt` publicly re-exports `ndarray` and `png` dependencies

### Changed

- **braking** The crate was split into two: the executable was moved into a separate crate `light-curve-dmdt-exec`
- The project repository was split from other 'light-curve*' crates and moved into <https://gituhb.com/light-curve/light-curve-dmdt>

## [0.5.0]

### Added

- README.md
- Documentation comments

### Changed

- Remove unused `Normalisable` trait
- Make `is_sorted()` private

## [0.4.0]

### Added

- This `CHANGELOG.md` file

### Changed

- `ndarray` is updated to 0.15.3, it is a breaking change


## [0.3.x]

—

## [0.2.x]

—

## [0.1.x]

—
