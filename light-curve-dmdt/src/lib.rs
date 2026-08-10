// Licu, the project's mascot, in the rustdoc sidebar and the browser tab.
// These have to be URLs: rustdoc drops them straight into <img src> and
// <link rel="icon">, and copies nothing next to the generated pages, so a path
// relative to the repository would resolve against docs.rs and 404. They point
// at this repository's own copy under assets/logo/ rather than at the branding
// repository, so the two cannot drift apart.
//
// The adaptive mark is the right file for both: it is loaded as an isolated
// document that the surrounding page cannot style, so the only way it can suit
// a light and a dark theme is by switching on prefers-color-scheme itself.
#![doc(
    html_logo_url = "https://raw.githubusercontent.com/light-curve/light-curve-dmdt/master/assets/logo/mark-adaptive.svg",
    html_favicon_url = "https://raw.githubusercontent.com/light-curve/light-curve-dmdt/master/assets/logo/mark-adaptive.svg"
)]
#![cfg_attr(feature = "doc-images",
cfg_attr(all(),
doc = ::embed_doc_image::embed_image!("example_png", "example.png")))]
#![doc = include_str!("../README.md")]

pub use crate::dmdt::*;
pub use crate::erf::*;
pub use crate::float_trait::Float;
pub use crate::grid::*;
#[cfg(feature = "png")]
pub use crate::images::{png, to_png};

pub use ndarray;

mod dmdt;
mod erf;
mod float_trait;
mod grid;
#[cfg(feature = "png")]
mod images;
mod util;
