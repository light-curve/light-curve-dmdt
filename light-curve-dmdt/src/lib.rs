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
