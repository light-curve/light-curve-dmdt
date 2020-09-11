use crate::evaluator::*;

#[derive(Clone, Default)]
pub struct TimeStandardDeviation {}

impl TimeStandardDeviation {
    pub fn new() -> Self {
        Self {}
    }
}

lazy_info!(
    TIME_STANDARD_DEVIATION_INFO,
    size: 1,
    min_ts_length: 2,
    t_required: true,
    m_required: false,
    w_required: false,
    sorting_required: false,
);

impl<T> FeatureEvaluator<T> for TimeStandardDeviation
where
    T: Float,
{
    fn eval(&self, ts: &mut TimeSeries<T>) -> Result<Vec<T>, EvaluatorError> {
        self.check_ts_length(ts)?;
        Ok(vec![ts.t.get_std()])
    }

    fn get_info(&self) -> &EvaluatorInfo {
        &TIME_STANDARD_DEVIATION_INFO
    }

    fn get_names(&self) -> Vec<&str> {
        vec!["ANTIFEATURE_time_standard_deviation"]
    }
}

#[cfg(test)]
#[allow(clippy::unreadable_literal)]
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;
    use crate::tests::*;

    eval_info_test!(
        time_standard_deviation_info,
        TimeStandardDeviation::default()
    );

    feature_test!(
        time_standard_deviation,
        [Box::new(TimeStandardDeviation::new())],
        [1.5811388300841898],
        [0.0_f32, 1.0, 2.0, 3.0, 4.0],
    );
}
