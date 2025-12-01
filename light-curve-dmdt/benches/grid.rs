use conv::*;
use criterion::{Criterion, black_box};
use light_curve_dmdt::{DmDt, Float, Grid, GridTrait, LinearGrid};
use ndarray::Array1;

pub fn bench_linear_grid_idx(c: &mut Criterion) {
    let values = [-0.1_f32, 0.0, 0.5, 0.9, 1.0, 1.1];

    let linear_grid = LinearGrid::new(0.0_f32, 1.0, 33);
    let linear_grid_enum = Grid::Linear(linear_grid.clone());

    c.bench_function("LinearGrid::idx", |b| {
        b.iter(|| {
            for &x in values.iter() {
                black_box(black_box(&linear_grid).idx(black_box(x)));
            }
        })
    });
    c.bench_function("wrapped LinearGrid::idx", |b| {
        b.iter(|| {
            for &x in values.iter() {
                black_box(black_box(&linear_grid_enum).idx(black_box(x)));
            }
        })
    });
}

pub fn bench_log_linear_grids<T>(c: &mut Criterion)
where
    T: Float + ValueFrom<f32>,
{
    let dmdt_lg_linear = DmDt::from_lgdt_dm_limits(
        T::zero(),
        2.0_f32.value_as::<T>().unwrap(),
        32,
        1.25_f32.value_as::<T>().unwrap(),
        32,
    );
    let dmdt_arrays: DmDt<T> = DmDt {
        dt_grid: Grid::array(Array1::logspace(
            10.0_f32.value_as::<T>().unwrap(),
            0.0_f32.value_as::<T>().unwrap(),
            2.0_f32.value_as::<T>().unwrap(),
            33,
        ))
        .unwrap(),
        dm_grid: Grid::array(Array1::linspace(
            -(1.25_f32.value_as::<T>().unwrap()),
            1.25_f32.value_as::<T>().unwrap(),
            33,
        ))
        .unwrap(),
    };

    let t = Array1::linspace(T::zero(), 100.0_f32.value_as::<T>().unwrap(), 101);
    let m = t.mapv(T::sin);

    c.bench_function(
        format!(
            "DmDt<{}>::points, log and linear grids",
            std::any::type_name::<T>()
        )
        .as_str(),
        |b| {
            b.iter(|| {
                black_box(
                    dmdt_lg_linear.points(black_box(t.as_slice().unwrap()), m.as_slice().unwrap()),
                );
            })
        },
    );

    c.bench_function(
        format!("DmDt<{}>::points, array grids", std::any::type_name::<T>()).as_str(),
        |b| {
            b.iter(|| {
                black_box(
                    dmdt_arrays.points(black_box(t.as_slice().unwrap()), m.as_slice().unwrap()),
                );
            })
        },
    );
}
