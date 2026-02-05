use crate::Float;

use ndarray::Array1;
use thiserror::Error;

/// Error types for Bayesian Blocks algorithm
#[derive(Error, Debug)]
pub enum BayesianBlocksError {
    #[error("input data is empty")]
    EmptyData,
    #[error("input data must have at least 2 points")]
    InsufficientData,
    #[error("prior parameter p0 must be in (0, 1), got {0}")]
    InvalidP0(f64),
    #[error("prior parameter gamma must be positive, got {0}")]
    InvalidGamma(f64),
    #[error("ncp_prior must be non-negative, got {0}")]
    InvalidNcpPrior(f64),
}

/// Prior specification for Bayesian Blocks
///
/// Controls the penalty for adding change points. Higher values = fewer blocks.
#[derive(Clone, Debug)]
pub enum Prior<T> {
    /// False alarm probability. Typical value: 0.05
    /// ncp_prior = 4 - ln(73.53 * p0 * N^(-0.478))
    P0(T),
    /// Geometric prior on number of blocks: P(N_blocks) ∝ gamma^N_blocks
    /// ncp_prior = -ln(gamma)
    Gamma(T),
    /// Direct specification of the prior penalty per change point
    NcpPrior(T),
}

impl<T: Float> Default for Prior<T> {
    fn default() -> Self {
        // Default p0 = 0.05 is a common choice
        Prior::P0(T::from(0.05).unwrap())
    }
}

/// Fitness function type for different data
#[derive(Clone, Copy, Debug, Default)]
pub enum FitnessFunc {
    /// For event data (time series), uses N_k * ln(N_k / T_k)
    /// Based on Eq. 19 from Scargle 2013
    #[default]
    Events,
    /// For point measurements with Gaussian errors
    /// Based on Eq. 41 from Scargle 2013
    PointMeasures,
}

/// Bayesian Blocks algorithm for optimal histogram binning
///
/// Implements the algorithm from Scargle et al. 2013
/// "Studies in Astronomical Time Series Analysis. VI. Bayesian Block Representations"
/// <https://doi.org/10.1088/0004-637X/764/2/167>
///
/// This finds optimal adaptive-width bin edges for histogramming data.
#[derive(Clone, Debug)]
pub struct BayesianBlocks<T> {
    prior: Prior<T>,
    fitness: FitnessFunc,
}

impl<T: Float> Default for BayesianBlocks<T> {
    fn default() -> Self {
        Self {
            prior: Prior::default(),
            fitness: FitnessFunc::default(),
        }
    }
}

impl<T: Float> BayesianBlocks<T> {
    /// Create a new BayesianBlocks with default settings
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the prior parameter
    pub fn with_prior(mut self, prior: Prior<T>) -> Self {
        self.prior = prior;
        self
    }

    /// Set the fitness function type
    pub fn with_fitness(mut self, fitness: FitnessFunc) -> Self {
        self.fitness = fitness;
        self
    }

    /// Compute ncp_prior from Prior specification
    fn compute_ncp_prior(&self, n: usize) -> Result<T, BayesianBlocksError> {
        match self.prior {
            Prior::P0(p0) => {
                let p0_f64 = p0.to_f64().unwrap();
                if p0_f64 <= 0.0 || p0_f64 >= 1.0 {
                    return Err(BayesianBlocksError::InvalidP0(p0_f64));
                }
                // ncp_prior = 4 - ln(73.53 * p0 * N^(-0.478))
                let n_f64 = n as f64;
                let ncp = 4.0 - (73.53 * p0_f64 * n_f64.powf(-0.478)).ln();
                Ok(T::from(ncp).unwrap())
            }
            Prior::Gamma(gamma) => {
                let gamma_f64 = gamma.to_f64().unwrap();
                if gamma_f64 <= 0.0 {
                    return Err(BayesianBlocksError::InvalidGamma(gamma_f64));
                }
                Ok(-gamma.ln())
            }
            Prior::NcpPrior(ncp) => {
                let ncp_f64 = ncp.to_f64().unwrap();
                if ncp_f64 < 0.0 {
                    return Err(BayesianBlocksError::InvalidNcpPrior(ncp_f64));
                }
                Ok(ncp)
            }
        }
    }

    /// Find optimal bin edges for event data (time series)
    ///
    /// # Arguments
    /// * `t` - Time values (must be sorted in ascending order)
    ///
    /// # Returns
    /// Array of bin edges defining optimal segmentation
    pub fn find_bins(&self, t: &[T]) -> Result<Array1<T>, BayesianBlocksError> {
        match self.fitness {
            FitnessFunc::Events => self.find_bins_events(t),
            FitnessFunc::PointMeasures => {
                // For point measures without explicit values/errors, treat as events
                self.find_bins_events(t)
            }
        }
    }

    /// Find optimal bin edges for point measurements with values and errors
    ///
    /// # Arguments
    /// * `t` - Time values (must be sorted in ascending order)
    /// * `x` - Measured values at each time
    /// * `sigma` - Measurement errors (standard deviations)
    ///
    /// # Returns
    /// Array of bin edges defining optimal segmentation
    pub fn find_bins_with_errors(
        &self,
        t: &[T],
        x: &[T],
        sigma: &[T],
    ) -> Result<Array1<T>, BayesianBlocksError> {
        if t.is_empty() {
            return Err(BayesianBlocksError::EmptyData);
        }
        if t.len() < 2 {
            return Err(BayesianBlocksError::InsufficientData);
        }

        let n = t.len();
        let ncp_prior = self.compute_ncp_prior(n)?;

        // Compute edges (including endpoints)
        let mut edges = Array1::zeros(n + 1);
        edges[0] = t[0];
        for i in 1..n {
            edges[i] = T::half() * (t[i - 1] + t[i]);
        }
        edges[n] = t[n - 1];

        // Precompute variance (sigma^2)
        let var: Vec<T> = sigma.iter().map(|&s| s * s).collect();

        // For point measures: fitness = (b_k^2) / (4 * a_k)
        // where a_k = sum(1/sigma_i^2) and b_k = sum(x_i/sigma_i^2)
        // Precompute cumulative sums for efficiency
        let mut a_sum = vec![T::zero(); n + 1]; // cumsum of 1/var
        let mut b_sum = vec![T::zero(); n + 1]; // cumsum of x/var
        for i in 0..n {
            let inv_var = T::one() / var[i];
            a_sum[i + 1] = a_sum[i] + inv_var;
            b_sum[i + 1] = b_sum[i] + x[i] * inv_var;
        }

        // Dynamic programming arrays
        let mut best = vec![T::neg_infinity(); n];
        let mut last = vec![0usize; n];

        // Main loop
        for r in 0..n {
            // Compute fitness for blocks [k..=r] for all k in 0..=r
            let mut fit_vec = vec![T::zero(); r + 1];
            for k in 0..=r {
                let a_k = a_sum[r + 1] - a_sum[k];
                let b_k = b_sum[r + 1] - b_sum[k];
                if a_k > T::zero() {
                    fit_vec[k] = (b_k * b_k) / (T::from(4.0).unwrap() * a_k);
                }
            }

            // Apply prior penalty and add previous best
            let mut a_r = vec![T::zero(); r + 1];
            for k in 0..=r {
                a_r[k] = fit_vec[k] - ncp_prior;
                if k > 0 {
                    a_r[k] += best[k - 1];
                }
            }

            // Find best configuration
            let (i_max, max_val) = a_r
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .unwrap();

            last[r] = i_max;
            best[r] = *max_val;
        }

        // Backtrack to find change points
        let mut change_points = vec![n];
        let mut i = n - 1;
        loop {
            change_points.push(last[i]);
            if last[i] == 0 {
                break;
            }
            i = last[i] - 1;
        }
        change_points.reverse();

        // Extract edges at change points
        let result: Array1<T> = change_points.iter().map(|&i| edges[i]).collect();
        Ok(result)
    }

    /// Internal: Find bins for event data
    fn find_bins_events(&self, t: &[T]) -> Result<Array1<T>, BayesianBlocksError> {
        if t.is_empty() {
            return Err(BayesianBlocksError::EmptyData);
        }
        if t.len() < 2 {
            return Err(BayesianBlocksError::InsufficientData);
        }

        let n = t.len();
        let ncp_prior = self.compute_ncp_prior(n)?;

        // Compute edges (cell boundaries between data points)
        let mut edges = Array1::zeros(n + 1);
        edges[0] = t[0];
        for i in 1..n {
            edges[i] = T::half() * (t[i - 1] + t[i]);
        }
        edges[n] = t[n - 1];

        // Block lengths from each edge to the end
        let last_edge = edges[n];
        let block_length: Vec<T> = edges.iter().map(|&e| last_edge - e).collect();

        // Dynamic programming arrays
        let mut best = vec![T::neg_infinity(); n];
        let mut last_idx = vec![0usize; n];

        // Main loop: iterate through each position
        for r in 0..n {
            // For events fitness: N_k * ln(N_k / T_k)
            // N_k = number of points in block [k..=r] = r - k + 1
            // T_k = block width = block_length[k] - block_length[r+1]

            let mut fit_vec = vec![T::zero(); r + 1];
            for k in 0..=r {
                let n_k = T::from(r - k + 1).unwrap();
                let t_k = block_length[k] - block_length[r + 1];
                if t_k > T::zero() && n_k > T::zero() {
                    fit_vec[k] = n_k * (n_k / t_k).ln();
                }
            }

            // Apply prior penalty and add previous best
            let mut a_r = vec![T::zero(); r + 1];
            for k in 0..=r {
                a_r[k] = fit_vec[k] - ncp_prior;
                if k > 0 {
                    a_r[k] += best[k - 1];
                }
            }

            // Find best configuration
            let (i_max, max_val) = a_r
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .unwrap();

            last_idx[r] = i_max;
            best[r] = *max_val;
        }

        // Backtrack to find change points
        let mut change_points = vec![n];
        let mut i = n - 1;
        loop {
            change_points.push(last_idx[i]);
            if last_idx[i] == 0 {
                break;
            }
            i = last_idx[i] - 1;
        }
        change_points.reverse();

        // Extract edges at change points
        let result: Array1<T> = change_points.iter().map(|&i| edges[i]).collect();
        Ok(result)
    }
}

/// Convenience function to compute Bayesian Blocks bin edges for event data
///
/// # Arguments
/// * `t` - Time values (must be sorted in ascending order)
/// * `p0` - False alarm probability (typical: 0.05). Lower values = fewer blocks.
///
/// # Returns
/// Array of bin edges defining optimal segmentation
///
/// # Example
/// ```
/// use light_curve_dmdt::bayesian_blocks;
///
/// let times: Vec<f64> = (0..100).map(|i| i as f64).collect();
/// let edges = bayesian_blocks(&times, 0.05).unwrap();
/// ```
pub fn bayesian_blocks<T: Float>(t: &[T], p0: f64) -> Result<Array1<T>, BayesianBlocksError> {
    BayesianBlocks::new()
        .with_prior(Prior::P0(T::from(p0).unwrap()))
        .find_bins(t)
}

/// Convenience function to compute Bayesian Blocks bin edges for point measurements
///
/// # Arguments
/// * `t` - Time values (must be sorted in ascending order)
/// * `x` - Measured values
/// * `sigma` - Measurement errors
/// * `p0` - False alarm probability (typical: 0.05)
///
/// # Returns
/// Array of bin edges defining optimal segmentation
pub fn bayesian_blocks_with_errors<T: Float>(
    t: &[T],
    x: &[T],
    sigma: &[T],
    p0: f64,
) -> Result<Array1<T>, BayesianBlocksError> {
    BayesianBlocks::new()
        .with_prior(Prior::P0(T::from(p0).unwrap()))
        .with_fitness(FitnessFunc::PointMeasures)
        .find_bins_with_errors(t, x, sigma)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_data() {
        // Uniformly spaced data should produce few blocks
        let t: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let edges = bayesian_blocks(&t, 0.05).unwrap();

        // With uniform data, we expect relatively few change points
        assert!(edges.len() >= 2); // At least start and end
        assert!(edges.len() <= 10); // Shouldn't have too many blocks
    }

    #[test]
    fn test_rate_change() {
        // Data with rate change should detect it
        let mut t: Vec<f64> = Vec::new();
        // First half: sparse (step = 2)
        for i in 0..50 {
            t.push(i as f64 * 2.0);
        }
        // Second half: dense (step = 0.5)
        let offset = 100.0;
        for i in 0..100 {
            t.push(offset + i as f64 * 0.5);
        }

        let edges = bayesian_blocks(&t, 0.05).unwrap();

        // Should detect at least one change point around the rate change
        assert!(edges.len() >= 3); // Start, change point, end

        // The change should be detected somewhere around t=100
        let has_change_near_100 = edges
            .iter()
            .skip(1)
            .take(edges.len() - 2)
            .any(|&e| e > 90.0 && e < 110.0);
        assert!(has_change_near_100);
    }

    #[test]
    fn test_empty_data() {
        let t: Vec<f64> = vec![];
        let result = bayesian_blocks(&t, 0.05);
        assert!(matches!(result, Err(BayesianBlocksError::EmptyData)));
    }

    #[test]
    fn test_insufficient_data() {
        let t: Vec<f64> = vec![1.0];
        let result = bayesian_blocks(&t, 0.05);
        assert!(matches!(result, Err(BayesianBlocksError::InsufficientData)));
    }

    #[test]
    fn test_invalid_p0() {
        let t: Vec<f64> = (0..10).map(|i| i as f64).collect();

        let result = BayesianBlocks::new()
            .with_prior(Prior::P0(0.0))
            .find_bins(&t);
        assert!(matches!(result, Err(BayesianBlocksError::InvalidP0(_))));

        let result = BayesianBlocks::new()
            .with_prior(Prior::P0(1.0))
            .find_bins(&t);
        assert!(matches!(result, Err(BayesianBlocksError::InvalidP0(_))));
    }

    #[test]
    fn test_point_measures() {
        // Test with point measurements
        let t: Vec<f64> = (0..50).map(|i| i as f64).collect();
        let x: Vec<f64> = t.iter().map(|&ti| (ti * 0.1).sin()).collect();
        let sigma: Vec<f64> = vec![0.1; 50];

        let edges = bayesian_blocks_with_errors(&t, &x, &sigma, 0.05).unwrap();
        assert!(edges.len() >= 2);
    }

    #[test]
    fn test_edges_sorted() {
        let t: Vec<f64> = (0..100).map(|i| i as f64 * 0.5).collect();
        let edges = bayesian_blocks(&t, 0.05).unwrap();

        // Edges should be sorted
        for i in 1..edges.len() {
            assert!(edges[i] > edges[i - 1], "Edges must be strictly ascending");
        }
    }

    #[test]
    fn test_edges_contain_bounds() {
        let t: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let edges = bayesian_blocks(&t, 0.05).unwrap();

        // First edge should be at or before first data point
        assert!(edges[0] <= t[0] + 1e-10);
        // Last edge should be at or after last data point
        assert!(edges[edges.len() - 1] >= t[t.len() - 1] - 1e-10);
    }

    #[test]
    fn test_gamma_prior() {
        let t: Vec<f64> = (0..100).map(|i| i as f64).collect();

        // gamma = 1 means no penalty for change points (many blocks)
        // gamma << 1 means strong penalty (few blocks)
        let edges_low_gamma = BayesianBlocks::new()
            .with_prior(Prior::Gamma(0.01))
            .find_bins(&t)
            .unwrap();

        let edges_high_gamma = BayesianBlocks::new()
            .with_prior(Prior::Gamma(0.9))
            .find_bins(&t)
            .unwrap();

        // Higher gamma should give more blocks (or equal)
        assert!(edges_high_gamma.len() >= edges_low_gamma.len());
    }

    #[test]
    fn test_f32_support() {
        let t: Vec<f32> = (0..50).map(|i| i as f32).collect();
        let edges = bayesian_blocks(&t, 0.05).unwrap();
        assert!(edges.len() >= 2);
    }
}
