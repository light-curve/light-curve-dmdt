use crate::Float;

use conv::{ConvAsUtil, ConvUtil, RoundToZero};
use enum_dispatch::enum_dispatch;
use ndarray::Array1;
use std::fmt::Debug;
use thiserror::Error;

/// Grid trait for dm or dt axis
#[enum_dispatch]
pub trait GridTrait<T>: Clone + Debug + Send + Sync
where
    T: Copy,
{
    /// Cell borders coordinates, [cell_count()](GridTrait::cell_count) + 1 length [Array1]
    fn get_borders(&self) -> &Array1<T>;

    /// Number of cells
    fn cell_count(&self) -> usize {
        self.get_borders().len() - 1
    }

    /// Coordinate of the left border of the leftmost cell
    fn get_start(&self) -> T {
        self.get_borders()[0]
    }

    /// Coordinate of the right border of the rightmost cell
    fn get_end(&self) -> T {
        self.get_borders()[self.cell_count()]
    }

    /// Get index of the cell containing given value
    ///
    /// Note that cells include their left borders but doesn't include right borders
    fn idx(&self, x: T) -> CellIndex;
}

/// Grid for dm or dt axis
#[enum_dispatch(GridTrait<T>)]
#[derive(Clone, Debug)]
#[non_exhaustive]
pub enum Grid<T>
where
    T: Float,
    Array1<T>: Clone + Debug,
{
    Array(ArrayGrid<T>),
    Linear(LinearGrid<T>),
    Lg(LgGrid<T>),
}

impl<T> Grid<T>
where
    T: Float,
    Array1<T>: Clone + Debug,
{
    pub fn array(borders: Array1<T>) -> Result<Self, ArrayGridError> {
        ArrayGrid::new(borders).map(Into::into)
    }

    pub fn linear(start: T, end: T, n: usize) -> Self {
        LinearGrid::new(start, end, n).into()
    }

    pub fn log_from_start_end(start: T, end: T, n: usize) -> Self {
        LgGrid::from_start_end(start, end, n).into()
    }

    pub fn log_from_lg_start_end(lg_start: T, lg_end: T, n: usize) -> Self {
        LgGrid::from_lg_start_end(lg_start, lg_end, n).into()
    }
}

/// An error to be returned from grid constructors
#[derive(Error, Debug)]
pub enum ArrayGridError {
    #[error("given grid is empty")]
    ArrayIsEmpty,
    #[error("given grid is not ascending")]
    ArrayIsNotAscending,
}

/// Grid which cell borders are defined by an ascending array
///
/// Lookup time is O(lb n)
#[derive(Clone, Debug)]
pub struct ArrayGrid<T> {
    borders: Array1<T>,
}

impl<T> ArrayGrid<T>
where
    Array1<T>: Clone + Debug,
    T: Float,
{
    /// Wraps given array into [ArrayGrid] or return an error
    ///
    /// Note that array describes cell borders, not center or whatever else
    pub fn new(borders: Array1<T>) -> Result<Self, ArrayGridError> {
        if borders.is_empty() {
            return Err(ArrayGridError::ArrayIsEmpty);
        }
        if !crate::util::is_sorted(borders.as_slice().unwrap()) {
            return Err(ArrayGridError::ArrayIsNotAscending);
        }
        Ok(Self { borders })
    }
}

impl<T> GridTrait<T> for ArrayGrid<T>
where
    Array1<T>: Clone + Debug,
    T: Float,
{
    #[inline]
    fn get_borders(&self) -> &Array1<T> {
        &self.borders
    }

    fn idx(&self, x: T) -> CellIndex {
        let i = self
            .borders
            .as_slice()
            .unwrap()
            .partition_point(|&b| b <= x);
        match i {
            0 => CellIndex::LowerMin,
            _ if i == self.borders.len() => CellIndex::GreaterMax,
            _ => CellIndex::Value(i - 1),
        }
    }
}

/// Linear grid defined by its start, end and number of cells
///
/// Lookup time is O(1)
#[derive(Clone, Debug)]
pub struct LinearGrid<T> {
    start: T,
    end: T,
    n: usize,
    cell_size: T,
    borders: Array1<T>,
}

impl<T> LinearGrid<T>
where
    T: Float,
{
    /// Create [LinearGrid] from borders and number of cells
    ///
    /// `start` is the left border of the leftmost cell, `end` is the right border of the rightmost
    /// cell, `n` is the number of cells. This means that the number of borders is `n + 1`, `start`
    /// border has zero index and `end` border has index `n`.
    pub fn new(start: T, end: T, n: usize) -> Self {
        assert!(end > start);
        let cell_size = (end - start) / n.value_as::<T>().unwrap();
        let borders = Array1::linspace(start, end, n + 1);
        Self {
            start,
            end,
            n,
            cell_size,
            borders,
        }
    }

    /// Cell size
    #[inline]
    pub fn get_cell_size(&self) -> T {
        self.cell_size
    }
}

impl<T> GridTrait<T> for LinearGrid<T>
where
    T: Float,
{
    #[inline]
    fn get_borders(&self) -> &Array1<T> {
        &self.borders
    }

    #[inline]
    fn cell_count(&self) -> usize {
        self.n
    }

    #[inline]
    fn get_start(&self) -> T {
        self.start
    }

    #[inline]
    fn get_end(&self) -> T {
        self.end
    }

    fn idx(&self, x: T) -> CellIndex {
        if x < self.start {
            return CellIndex::LowerMin;
        }
        if x >= self.end {
            return CellIndex::GreaterMax;
        }
        let i = ((x - self.start) / self.cell_size)
            .approx_by::<RoundToZero>()
            .unwrap();
        if i < self.n {
            CellIndex::Value(i)
        } else {
            // x is a bit smaller self.end + float rounding
            CellIndex::Value(self.n - 1)
        }
    }
}

/// Logarithmic grid defined by its start, end and number of cells
///
/// Lookup time is O(1)
#[derive(Clone, Debug)]
pub struct LgGrid<T>
where
    T: Copy,
{
    start: T,
    end: T,
    lg_start: T,
    lg_end: T,
    n: usize,
    cell_lg_size: T,
    borders: Array1<T>,
}

impl<T> LgGrid<T>
where
    T: Float,
{
    /// Create [LinearGrid] from borders and number of cells
    ///
    /// `start` is the left border of the leftmost cell, `end` is the right border of the rightmost
    /// cell, `n` is the number of cells. This means that the number of borders is `n + 1`, `start`
    /// border has zero index and `end` border has index `n`.
    pub fn from_start_end(start: T, end: T, n: usize) -> Self {
        assert!(end > start);
        assert!(start.is_positive());
        let lg_start = start.log10();
        let lg_end = end.log10();
        let cell_lg_size = (lg_end - lg_start) / n.value_as::<T>().unwrap();
        let mut borders = Array1::logspace(T::ten(), lg_start, lg_end, n + 1);
        borders[0] = start;
        borders[n] = end;
        Self {
            start,
            end,
            lg_start,
            lg_end,
            n,
            cell_lg_size,
            borders,
        }
    }

    /// Create [LinearGrid] from decimal logarithms of borders and number of cells
    ///
    /// `lg_start` is the decimal logarithm of the left border of the leftmost cell, `lg_end` is the
    /// decimal logarithm of the right border of the rightmost cell, `n` is the number of cells.
    /// This means that the number of borders is `n + 1`, `lg_start` border has zero index and
    /// `lg_end` border has index `n`.
    pub fn from_lg_start_end(lg_start: T, lg_end: T, n: usize) -> Self {
        Self::from_start_end(T::powf(T::ten(), lg_start), T::powf(T::ten(), lg_end), n)
    }

    /// Logarithmic size of cell
    #[inline]
    pub fn get_cell_lg_size(&self) -> T {
        self.cell_lg_size
    }

    /// Logarithm of the leftmost border
    #[inline]
    pub fn get_lg_start(&self) -> T {
        self.lg_start
    }

    /// Logarithm of the rightmost border
    #[inline]
    pub fn get_lg_end(&self) -> T {
        self.lg_end
    }
}

impl<T> GridTrait<T> for LgGrid<T>
where
    T: Float,
{
    #[inline]
    fn get_borders(&self) -> &Array1<T> {
        &self.borders
    }

    #[inline]
    fn cell_count(&self) -> usize {
        self.n
    }

    #[inline]
    fn get_start(&self) -> T {
        self.start
    }

    #[inline]
    fn get_end(&self) -> T {
        self.end
    }

    fn idx(&self, x: T) -> CellIndex {
        if x < self.start {
            return CellIndex::LowerMin;
        }
        if x >= self.end {
            return CellIndex::GreaterMax;
        }
        let i = ((x.log10() - self.lg_start) / self.cell_lg_size)
            .approx_by::<RoundToZero>()
            .unwrap();
        if i < self.n {
            CellIndex::Value(i)
        } else {
            // x is a bit smaller self.end + float rounding
            CellIndex::Value(self.n - 1)
        }
    }
}

/// Value to return from [GridTrait::idx]
pub enum CellIndex {
    /// Bellow the leftmost border
    LowerMin,
    /// Equal or greater the rightmost border
    GreaterMax,
    /// Cell index
    Value(usize),
}
