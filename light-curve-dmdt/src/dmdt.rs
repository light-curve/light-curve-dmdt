use crate::{CellIndex, ErfFloat, ErrorFunction, Float, Grid, GridTrait, LgGrid, LinearGrid};

use itertools::Itertools;
use ndarray::{s, Array1, Array2};
use std::fmt::Debug;

/// dm–dt map plotter
#[derive(Clone, Debug)]
pub struct DmDt<T>
where
    T: Float,
    Array1<T>: Clone + Debug,
{
    pub dt_grid: Grid<T>,
    pub dm_grid: Grid<T>,
}

impl<T> DmDt<T>
where
    T: Float,
    Array1<T>: Clone + Debug,
{
    /// Create new [DmDt]
    pub fn from_grids<Gdt, Gdm>(dt_grid: Gdt, dm_grid: Gdm) -> Self
    where
        Gdt: Into<Grid<T>>,
        Gdm: Into<Grid<T>>,
    {
        Self {
            dt_grid: dt_grid.into(),
            dm_grid: dm_grid.into(),
        }
    }

    /// Create new [DmDt] with logarithmic dt grid and linear dm grid
    ///
    /// dt grid will have borders `[10^min_lgdt, 10^max_lgdt)`, dm grid will have borders
    /// `[-max_abs_dm, max_abs_dm)`
    pub fn from_lgdt_dm_limits(
        min_lgdt: T,
        max_lgdt: T,
        lgdt_size: usize,
        max_abs_dm: T,
        dm_size: usize,
    ) -> Self {
        Self::from_grids(
            LgGrid::from_lg_start_end(min_lgdt, max_lgdt, lgdt_size),
            LinearGrid::new(-max_abs_dm, max_abs_dm, dm_size),
        )
    }

    /// N dt by N dm
    pub fn shape(&self) -> (usize, usize) {
        (self.dt_grid.cell_count(), self.dm_grid.cell_count())
    }

    /// Represents each pair of (t, m) points as a unity value in dm-dt map
    ///
    /// `t` must be an ascending slice
    pub fn points(&self, t: &[T], m: &[T]) -> Array2<u64> {
        let mut a = Array2::zeros(self.shape());
        for (i1, (&x1, &y1)) in t.iter().zip(m.iter()).enumerate() {
            for (&x2, &y2) in t[i1 + 1..].iter().zip(m[i1 + 1..].iter()) {
                let dt = x2 - x1;
                let idx_dt = match self.dt_grid.idx(dt) {
                    CellIndex::LowerMin => continue,
                    CellIndex::GreaterMax => break,
                    CellIndex::Value(idx_dt) => idx_dt,
                };
                let dm = y2 - y1;
                let idx_dm = match self.dm_grid.idx(dm) {
                    CellIndex::Value(idx_dm) => idx_dm,
                    CellIndex::LowerMin | CellIndex::GreaterMax => continue,
                };
                a[(idx_dt, idx_dm)] += 1;
            }
        }
        a
    }

    fn update_gausses_helper<Erf>(
        &self,
        a: &mut Array2<T>,
        idx_dt: usize,
        y1: T,
        y2: T,
        d1: T,
        d2: T,
    ) where
        T: ErfFloat,
        Erf: ErrorFunction<T>,
    {
        let dm = y2 - y1;
        let dm_err = T::sqrt(d1 + d2);

        let min_idx_dm = match self
            .dm_grid
            .idx(dm + Erf::min_dx_nonzero_normal_cdf(dm_err))
        {
            CellIndex::LowerMin => 0,
            CellIndex::GreaterMax => return,
            CellIndex::Value(min_idx_dm) => min_idx_dm,
        };
        let max_idx_dm = match self
            .dm_grid
            .idx(dm + Erf::max_dx_nonunity_normal_cdf(dm_err))
        {
            CellIndex::LowerMin => return,
            CellIndex::GreaterMax => self.dm_grid.cell_count(),
            CellIndex::Value(i) => usize::min(i + 1, self.dm_grid.cell_count()),
        };

        a.slice_mut(s![idx_dt, min_idx_dm..max_idx_dm])
            .iter_mut()
            .zip(
                self.dm_grid
                    .get_borders()
                    .slice(s![min_idx_dm..max_idx_dm + 1])
                    .iter()
                    .map(|&dm_border| Erf::normal_cdf(dm_border, dm, dm_err))
                    .tuple_windows()
                    .map(|(a, b)| b - a),
            )
            .for_each(|(cell, value)| *cell += value);
    }

    /// Represents each pair of (t, m, err2) points as a Gaussian distribution in dm-dt map
    ///
    /// `t` must be an ascending slice.
    ///
    /// Each observation is assumed to happen at time moment `t_i` and have Gaussian distribution of
    /// its magnitude `N(m_i, err2_i)`. Each pair of observations
    /// `(t_1, m_1, err2_1), (t_2, m_2, err2_2)` is represented by 1-D Gaussian in the dm-dt space
    /// having constant `dt` and `dm ~ N(m2-m1, err2_1 + err2_2)`. This distribution is integrated
    /// over each cell using `Erf` struct implementing [ErrorFunction].
    pub fn gausses<Erf>(&self, t: &[T], m: &[T], err2: &[T]) -> Array2<T>
    where
        T: ErfFloat,
        Erf: ErrorFunction<T>,
    {
        let mut a = Array2::zeros(self.shape());
        for (i1, ((&x1, &y1), &d1)) in t.iter().zip(m.iter()).zip(err2.iter()).enumerate() {
            for ((&x2, &y2), &d2) in t[i1 + 1..]
                .iter()
                .zip(m[i1 + 1..].iter())
                .zip(err2[i1 + 1..].iter())
            {
                let dt = x2 - x1;
                let idx_dt = match self.dt_grid.idx(dt) {
                    CellIndex::LowerMin => continue,
                    CellIndex::GreaterMax => break,
                    CellIndex::Value(idx_dt) => idx_dt,
                };
                self.update_gausses_helper::<Erf>(&mut a, idx_dt, y1, y2, d1, d2);
            }
        }
        a
    }

    /// Count dt in the each dt grid cell
    pub fn dt_points(&self, t: &[T]) -> Array1<u64> {
        let mut a = Array1::zeros(self.dt_grid.cell_count());
        for (i1, &x1) in t.iter().enumerate() {
            for &x2 in t[i1 + 1..].iter() {
                let dt = x2 - x1;
                let idx_dt = match self.dt_grid.idx(dt) {
                    CellIndex::LowerMin => continue,
                    CellIndex::GreaterMax => break,
                    CellIndex::Value(idx_dt) => idx_dt,
                };
                a[idx_dt] += 1;
            }
        }
        a
    }

    /// Conditional probability `p(m2-m1|t2-t1)`
    ///
    /// Technically it is optimized version of [DmDt::gausses()] normalized by [DmDt::dt_points] but
    /// with better performance. Mathematically it represents the distribution of conditional
    /// probability `p(m2-m1|t2-t1)`, see
    /// [Soraisam et al. 2020](https://doi.org/10.3847/1538-4357/ab7b61) for details.
    pub fn cond_prob<Erf>(&self, t: &[T], m: &[T], err2: &[T]) -> Array2<T>
    where
        T: ErfFloat,
        Erf: ErrorFunction<T>,
    {
        let mut a: Array2<T> = Array2::zeros(self.shape());
        let mut dt_points: Array1<u64> = Array1::zeros(self.dt_grid.cell_count());
        for (i1, ((&x1, &y1), &d1)) in t.iter().zip(m.iter()).zip(err2.iter()).enumerate() {
            for ((&x2, &y2), &d2) in t[i1 + 1..]
                .iter()
                .zip(m[i1 + 1..].iter())
                .zip(err2[i1 + 1..].iter())
            {
                let dt = x2 - x1;
                let idx_dt = match self.dt_grid.idx(dt) {
                    CellIndex::LowerMin => continue,
                    CellIndex::GreaterMax => break,
                    CellIndex::Value(idx_dt) => idx_dt,
                };

                dt_points[idx_dt] += 1;

                self.update_gausses_helper::<Erf>(&mut a, idx_dt, y1, y2, d1, d2);
            }
        }
        ndarray::Zip::from(a.rows_mut())
            .and(&dt_points)
            .for_each(|mut row, &count| {
                if count == 0 {
                    return;
                }
                row /= T::approx_from(count).unwrap();
            });
        a
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::dmdt::DmDt;
    use crate::erf::{Eps1Over1e3Erf, ExactErf};

    use approx::assert_abs_diff_eq;
    use ndarray::Axis;
    use static_assertions::assert_impl_all;

    assert_impl_all!(DmDt<f32>: Clone, Debug, Send, Sync);
    assert_impl_all!(DmDt<f64>: Clone, Debug, Send, Sync);

    #[test]
    fn dt_points_vs_points() {
        let dmdt = DmDt::from_lgdt_dm_limits(0.0_f32, 2.0_f32, 32, 3.0_f32, 32);
        let t = Array1::linspace(0.0, 100.0, 101);
        // dm is within map borders
        let m = t.mapv(f32::sin);

        let points = dmdt.points(t.as_slice().unwrap(), m.as_slice().unwrap());
        let dt_points = dmdt.dt_points(t.as_slice().unwrap());

        assert_eq!(points.sum_axis(Axis(1)), dt_points,);
    }

    #[test]
    fn dt_points_vs_gausses() {
        let dmdt = DmDt::from_lgdt_dm_limits(0.0_f32, 2.0_f32, 32, 3.0_f32, 32);
        let t = Array1::linspace(0.0, 100.0, 101);
        // dm is within map borders
        let m = t.mapv(f32::sin);
        // err is ~0.03
        let err2 = Array1::from_elem(101, 0.001_f32);

        let gausses = dmdt.gausses::<ExactErf>(
            t.as_slice().unwrap(),
            m.as_slice().unwrap(),
            err2.as_slice().unwrap(),
        );
        let sum_gausses = gausses.sum_axis(Axis(1));
        let dt_points = dmdt.dt_points(t.as_slice().unwrap()).mapv(|x| x as f32);

        assert_abs_diff_eq!(
            sum_gausses.as_slice().unwrap(),
            dt_points.as_slice().unwrap(),
            epsilon = 1e-4,
        );
    }

    #[test]
    fn cond_prob() {
        let dmdt = DmDt::from_lgdt_dm_limits(0.0_f32, 2.0_f32, 32, 1.25_f32, 32);

        let t = Array1::linspace(0.0, 100.0, 101);
        let m = t.mapv(f32::sin);
        // err is ~0.03
        let err2 = Array1::from_elem(101, 0.001);

        let from_gausses_dt_points = {
            let mut map = dmdt.gausses::<Eps1Over1e3Erf>(
                t.as_slice().unwrap(),
                m.as_slice().unwrap(),
                err2.as_slice_memory_order().unwrap(),
            );
            let dt_points = dmdt.dt_points(t.as_slice().unwrap());
            let dt_non_zero_points = dt_points.mapv(|x| if x == 0 { 1.0 } else { x as f32 });
            map /= &dt_non_zero_points.into_shape((map.nrows(), 1)).unwrap();
            map
        };

        let from_cond_prob = dmdt.cond_prob::<Eps1Over1e3Erf>(
            t.as_slice().unwrap(),
            m.as_slice().unwrap(),
            err2.as_slice().unwrap(),
        );

        assert_abs_diff_eq!(
            from_gausses_dt_points.as_slice().unwrap(),
            from_cond_prob.as_slice().unwrap(),
            epsilon = std::f32::EPSILON,
        );
    }
}
