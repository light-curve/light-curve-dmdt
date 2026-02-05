use crate::Float;

use ndarray::Array1;
use thiserror::Error;

/// Error types for Bayesian Blocks algorithm
#[derive(Error, Debug)]
pub enum BayesianBlocksError {
    #[error("input data must have at least {0} points, got {1}")]
    InsufficientData(usize, usize),
    #[error("prior parameter p0 must be in (0, 1), got {0}")]
    InvalidP0(f64),
    #[error("prior parameter gamma must be positive, got {0}")]
    InvalidGamma(f64),
    #[error("ncp_prior must be non-negative, got {0}")]
    InvalidNcpPrior(f64),
}

/// Constants for the p0-to-ncp_prior conversion from Scargle et al. 2013.
///
/// From Section 3.3, Eq. 21 of "Studies in Astronomical Time Series Analysis. VI.
/// Bayesian Block Representations" (Scargle et al. 2013, ApJ 764:167):
/// <https://doi.org/10.1088/0004-637X/764/2/167>
///
/// The formula is: ncp_prior = 4 - ln(73.53 * p0 * N^(-0.478))
///
/// These constants were determined empirically by the authors through simulations
/// to calibrate the false positive rate for the events fitness function.
mod p0_prior_constants {
    /// Additive constant in the ncp_prior formula
    pub const ADDITIVE: f64 = 4.0;
    /// Multiplicative factor inside the logarithm
    pub const MULTIPLIER: f64 = 73.53;
    /// Exponent for the data size N
    pub const EXPONENT: f64 = -0.478;
}

/// Prior specification for Bayesian Blocks
///
/// Controls the penalty for adding change points. Higher values = fewer blocks.
#[derive(Clone, Debug)]
pub enum Prior<T> {
    /// False alarm probability (typical: 0.05).
    ///
    /// Converted to ncp_prior using empirical calibration from Scargle et al. 2013, Eq. 21:
    /// `ncp_prior = 4 - ln(73.53 * p0 * N^(-0.478))`
    P0(T),
    /// Geometric prior on number of blocks: P(N_blocks) ∝ gamma^N_blocks
    ///
    /// Converted as: `ncp_prior = -ln(gamma)`
    Gamma(T),
    /// Direct specification of the prior penalty per change point
    NcpPrior(T),
}

impl<T: Float> Default for Prior<T> {
    fn default() -> Self {
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
    const MIN_DATA_POINTS: usize = 2;

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
                if !(0.0..1.0).contains(&p0_f64) || p0_f64 == 0.0 {
                    return Err(BayesianBlocksError::InvalidP0(p0_f64));
                }
                // Scargle et al. 2013, Eq. 21
                let ncp = p0_prior_constants::ADDITIVE
                    - (p0_prior_constants::MULTIPLIER
                        * p0_f64
                        * (n as f64).powf(p0_prior_constants::EXPONENT))
                    .ln();
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
            FitnessFunc::PointMeasures => self.find_bins_events(t),
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
        if t.len() < Self::MIN_DATA_POINTS {
            return Err(BayesianBlocksError::InsufficientData(
                Self::MIN_DATA_POINTS,
                t.len(),
            ));
        }

        let n = t.len();
        let ncp_prior = self.compute_ncp_prior(n)?;

        // Compute edges: first point, midpoints between consecutive points, last point
        let edges: Array1<T> = std::iter::once(t[0])
            .chain(t.windows(2).map(|w| T::half() * (w[0] + w[1])))
            .chain(std::iter::once(t[n - 1]))
            .collect();

        // Precompute cumulative sums for fitness calculation
        // a_sum[i] = sum(1/sigma_j^2) for j in 0..i
        // b_sum[i] = sum(x_j/sigma_j^2) for j in 0..i
        let (a_sum, b_sum): (Vec<T>, Vec<T>) = std::iter::once((T::zero(), T::zero()))
            .chain(sigma.iter().zip(x.iter()).scan(
                (T::zero(), T::zero()),
                |(a_acc, b_acc), (&s, &xi)| {
                    let inv_var = T::one() / (s * s);
                    *a_acc += inv_var;
                    *b_acc += xi * inv_var;
                    Some((*a_acc, *b_acc))
                },
            ))
            .unzip();

        // Dynamic programming: best[r] = best fitness ending at position r
        //                      last[r] = start index of best block ending at r
        let (_best, last) = (0..n).fold(
            (vec![T::neg_infinity(); n], vec![0usize; n]),
            |(mut best, mut last), r| {
                // Fitness for block [k..=r]: (b_k^2) / (4 * a_k) where
                // a_k = a_sum[r+1] - a_sum[k], b_k = b_sum[r+1] - b_sum[k]
                let (i_max, max_val) = (0..=r)
                    .map(|k| {
                        let a_k = a_sum[r + 1] - a_sum[k];
                        let b_k = b_sum[r + 1] - b_sum[k];
                        let fitness = if a_k > T::zero() {
                            (b_k * b_k) / (T::from(4.0).unwrap() * a_k)
                        } else {
                            T::zero()
                        };
                        let score =
                            fitness - ncp_prior + if k > 0 { best[k - 1] } else { T::zero() };
                        (k, score)
                    })
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .unwrap();

                last[r] = i_max;
                best[r] = max_val;
                (best, last)
            },
        );

        Ok(Self::backtrack_edges(&edges, &last, n))
    }

    /// Internal: Find bins for event data using dynamic programming
    fn find_bins_events(&self, t: &[T]) -> Result<Array1<T>, BayesianBlocksError> {
        if t.len() < Self::MIN_DATA_POINTS {
            return Err(BayesianBlocksError::InsufficientData(
                Self::MIN_DATA_POINTS,
                t.len(),
            ));
        }

        let n = t.len();
        let ncp_prior = self.compute_ncp_prior(n)?;

        // Compute edges: first point, midpoints between consecutive points, last point
        let edges: Array1<T> = std::iter::once(t[0])
            .chain(t.windows(2).map(|w| T::half() * (w[0] + w[1])))
            .chain(std::iter::once(t[n - 1]))
            .collect();

        // Block length from edge i to the last edge
        let last_edge = edges[n];
        let block_length: Vec<T> = edges.iter().map(|&e| last_edge - e).collect();

        // Dynamic programming
        let (_best, last) = (0..n).fold(
            (vec![T::neg_infinity(); n], vec![0usize; n]),
            |(mut best, mut last), r| {
                // Fitness for block [k..=r]: N_k * ln(N_k / T_k)
                // N_k = r - k + 1, T_k = block_length[k] - block_length[r+1]
                let (i_max, max_val) = (0..=r)
                    .map(|k| {
                        let n_k = T::from(r - k + 1).unwrap();
                        let t_k = block_length[k] - block_length[r + 1];
                        let fitness = if t_k > T::zero() {
                            n_k * (n_k / t_k).ln()
                        } else {
                            T::zero()
                        };
                        let score =
                            fitness - ncp_prior + if k > 0 { best[k - 1] } else { T::zero() };
                        (k, score)
                    })
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .unwrap();

                last[r] = i_max;
                best[r] = max_val;
                (best, last)
            },
        );

        Ok(Self::backtrack_edges(&edges, &last, n))
    }

    /// Backtrack through the `last` array to find change points and extract edges
    fn backtrack_edges(edges: &Array1<T>, last: &[usize], n: usize) -> Array1<T> {
        // Collect change points by backtracking: start at n-1, follow last[] chain
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

        change_points.into_iter().map(|idx| edges[idx]).collect()
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
    use approx::assert_abs_diff_eq;

    /// Test against astropy.stats.bayesian_blocks reference implementation
    ///
    /// Generated with:
    /// ```python
    /// from astropy.stats import bayesian_blocks
    /// import numpy as np
    /// t = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    /// edges = bayesian_blocks(t, fitness='events', p0=0.05)
    /// # Result: [0.0, 4.0]
    /// ```
    #[test]
    fn test_astropy_simple_uniform() {
        let t: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let edges = bayesian_blocks(&t, 0.05).unwrap();

        // astropy returns [0.0, 4.0] for this uniform data
        assert_eq!(edges.len(), 2);
        assert_abs_diff_eq!(edges[0], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(edges[1], 4.0, epsilon = 1e-10);
    }

    /// Test rate change detection against astropy reference
    ///
    /// Generated with:
    /// ```python
    /// from astropy.stats import bayesian_blocks
    /// import numpy as np
    /// t_sparse = np.arange(0, 50, 2.0)  # 25 points
    /// t_dense = 50.0 + np.arange(0, 25, 0.5)  # 50 points
    /// t = np.concatenate([t_sparse, t_dense])
    /// edges = bayesian_blocks(t, fitness='events', p0=0.05)
    /// # Result: [0.0, 50.25, 74.5]
    /// ```
    #[test]
    fn test_astropy_rate_change() {
        let t_sparse: Vec<f64> = (0..25).map(|i| i as f64 * 2.0).collect();
        let t_dense: Vec<f64> = (0..50).map(|i| 50.0 + i as f64 * 0.5).collect();
        let t: Vec<f64> = t_sparse.into_iter().chain(t_dense).collect();

        let edges = bayesian_blocks(&t, 0.05).unwrap();

        // astropy returns [0.0, 50.25, 74.5]
        assert_eq!(edges.len(), 3);
        assert_abs_diff_eq!(edges[0], 0.0, epsilon = 1e-10);
        assert_abs_diff_eq!(edges[1], 50.25, epsilon = 1e-10);
        assert_abs_diff_eq!(edges[2], 74.5, epsilon = 1e-10);
    }

    #[test]
    fn test_insufficient_data() {
        let empty: Vec<f64> = vec![];
        let result = bayesian_blocks(&empty, 0.05);
        assert!(matches!(
            result,
            Err(BayesianBlocksError::InsufficientData(2, 0))
        ));

        let single: Vec<f64> = vec![1.0];
        let result = bayesian_blocks(&single, 0.05);
        assert!(matches!(
            result,
            Err(BayesianBlocksError::InsufficientData(2, 1))
        ));
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

        let result = BayesianBlocks::new()
            .with_prior(Prior::P0(-0.1))
            .find_bins(&t);
        assert!(matches!(result, Err(BayesianBlocksError::InvalidP0(_))));
    }

    #[test]
    fn test_point_measures() {
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

        for i in 1..edges.len() {
            assert!(edges[i] > edges[i - 1], "Edges must be strictly ascending");
        }
    }

    #[test]
    fn test_edges_contain_bounds() {
        let t: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let edges = bayesian_blocks(&t, 0.05).unwrap();

        assert_abs_diff_eq!(edges[0], t[0], epsilon = 1e-10);
        assert_abs_diff_eq!(edges[edges.len() - 1], t[t.len() - 1], epsilon = 1e-10);
    }

    #[test]
    fn test_gamma_prior() {
        let t: Vec<f64> = (0..100).map(|i| i as f64).collect();

        // Lower gamma = stronger penalty = fewer blocks
        let edges_low_gamma = BayesianBlocks::new()
            .with_prior(Prior::Gamma(0.01))
            .find_bins(&t)
            .unwrap();

        // Higher gamma = weaker penalty = more blocks
        let edges_high_gamma = BayesianBlocks::new()
            .with_prior(Prior::Gamma(0.9))
            .find_bins(&t)
            .unwrap();

        assert!(edges_high_gamma.len() >= edges_low_gamma.len());
    }

    #[test]
    fn test_f32_support() {
        let t: Vec<f32> = (0..50).map(|i| i as f32).collect();
        let edges = bayesian_blocks(&t, 0.05).unwrap();
        assert!(edges.len() >= 2);
    }

    #[test]
    fn test_ncp_prior_direct() {
        let t: Vec<f64> = (0..50).map(|i| i as f64).collect();

        // Direct ncp_prior specification should work
        let edges = BayesianBlocks::new()
            .with_prior(Prior::NcpPrior(4.0))
            .find_bins(&t)
            .unwrap();
        assert!(edges.len() >= 2);
    }
}
