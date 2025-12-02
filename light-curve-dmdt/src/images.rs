use ndarray::{Array2, ArrayRef2};
pub use png;
use std::io::Write;

/// Convert [u8] dm–dt map into PNG image
pub fn to_png<W>(w: W, a: &ArrayRef2<u8>) -> Result<(), png::EncodingError>
where
    W: Write,
{
    let transposed = {
        let mut b = Array2::zeros((a.ncols(), a.nrows()));
        b.assign(&a.t());
        b
    };
    let mut encoder = png::Encoder::new(w, transposed.ncols() as u32, transposed.nrows() as u32);
    encoder.set_color(png::ColorType::Grayscale);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header()?;
    writer.write_image_data(transposed.as_slice().unwrap())?;
    Ok(())
}
