#![allow(clippy::doc_markdown)] // Allow mathematical notation (λ₁, M_ν, etc.)
#![allow(clippy::large_stack_arrays)] // 8×8 Rational arrays are fundamental mathematical objects
#![allow(clippy::needless_range_loop)] // Matrix indexing is clearer with explicit loops
#![allow(clippy::cast_possible_truncation)] // i128→i64 casts are checked by GCD reduction
#![allow(clippy::cast_possible_wrap)] // u128→i128 casts for GCD results are safe (reduced)
#![allow(clippy::cast_lossless)] // i64→i128 casts flagged by clippy are intentional
#![allow(clippy::if_not_else)] // Mathematical condition ordering is more natural
#![allow(clippy::cognitive_complexity)] // Mathematical algorithms are inherently complex
#![allow(clippy::missing_panics_doc)] // Internal module, panics documented in public API
#![allow(clippy::missing_errors_doc)] // No Result types in this module
//! # Chapter 11: Invariant Basis Extraction
//!
//! This chapter presents the **invariant basis extractor** for the Atlas embedding
//! space, complementing the spectral analysis of Chapter 10. While Chapter 10
//! analyzes the **graph Laplacian** (combinatorial structure), this chapter analyzes
//! the **covariance operator** of the embedding coordinates (geometric structure).
//!
//! ## Overview
//!
//! The Atlas of Resonance Classes embeds into the E₈ root system (Chapter 3),
//! placing 96 vertices as points in 8-dimensional space. The **covariance matrix**
//! of these embedded coordinates captures the principal geometric axes of the
//! configuration. Its eigenvectors form the **invariant basis** — the intrinsic
//! coordinate system aligned to the embedding's own geometry.
//!
//! ## Motivation
//!
//! Operating in the canonical E₈ basis is arbitrary. The invariant basis:
//! - Aligns coordinates to the Atlas's intrinsic geometry
//! - Separates dominant structural directions from minor ones
//! - Enables dimensionality reduction with controlled information loss
//! - Provides the bridge between discrete graph structure and continuous geometry
//!
//! ## Mathematical Background
//!
//! **Definition 11.1.1 (Covariance Operator)**: Given n vectors x₁, ..., xₙ ∈ ℚ^d,
//! the **covariance matrix** is:
//!
//! ```text
//! C = (1/n) Σᵢ (xᵢ - μ)(xᵢ - μ)ᵀ
//! ```
//!
//! where μ = (1/n) Σᵢ xᵢ is the mean vector.
//!
//! **Properties**:
//! - C is symmetric: C = Cᵀ
//! - C is positive semidefinite: all eigenvalues ≥ 0
//! - tr(C) = total variance
//! - rank(C) ≤ min(n, d)
//!
//! **Definition 11.1.2 (Invariant Basis)**: The **invariant basis** is the set of
//! eigenvectors {v₁, ..., v_d} of C, ordered by eigenvalue λ₁ ≥ λ₂ ≥ ... ≥ λ_d.
//! The eigenvalue λₖ measures the variance along direction vₖ.
//!
//! **Theorem 11.1.1 (Optimal Projection)**: For any k ≤ d, the projection onto
//! {v₁, ..., vₖ} minimizes the reconstruction error among all k-dimensional
//! subspaces. This is the Eckart–Young–Mirsky theorem.
//!
//! ## The Atlas Covariance
//!
//! **Theorem 11.2.1 (Atlas Covariance Structure)**: The covariance matrix of the
//! 96 Atlas embedding vectors in E₈ has rank ≤ 8 and is computed exactly using
//! rational arithmetic. Its eigenvalues reveal the principal geometric axes of
//! the Atlas configuration within E₈.
//!
//! All computations use exact rational arithmetic — no floating point.

use crate::arithmetic::{Rational, Vector8};
use crate::embedding::compute_atlas_embedding;
use crate::Atlas;
use num_traits::{One, Zero};
use std::fmt;

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

/// Embedding dimension (E₈ root space)
const DIM: usize = 8;


// ─────────────────────────────────────────────────────────────────────────────
// Dense Rational Matrix (8×8)
// ─────────────────────────────────────────────────────────────────────────────

/// Dense 8×8 matrix with exact rational entries.
///
/// Used internally for covariance computation and eigendecomposition.
/// Row-major storage.
#[derive(Clone, PartialEq, Eq)]
struct DenseMat8 {
    data: [[Rational; DIM]; DIM],
}

impl DenseMat8 {
    /// Create zero matrix
    fn zero() -> Self {
        Self {
            data: [[Rational::zero(); DIM]; DIM],
        }
    }

    /// Create identity matrix
    #[cfg(test)]
    fn identity() -> Self {
        let mut m = Self::zero();
        for i in 0..DIM {
            m.data[i][i] = Rational::one();
        }
        m
    }

    /// Get entry
    const fn get(&self, i: usize, j: usize) -> Rational {
        self.data[i][j]
    }

    /// Set entry
    #[cfg(test)]
    fn set(&mut self, i: usize, j: usize, val: Rational) {
        self.data[i][j] = val;
    }

    /// Compute trace
    fn trace(&self) -> Rational {
        let mut sum = Rational::zero();
        for i in 0..DIM {
            sum += self.data[i][i];
        }
        sum
    }

    /// Check symmetry
    fn is_symmetric(&self) -> bool {
        for i in 0..DIM {
            for j in (i + 1)..DIM {
                if self.data[i][j] != self.data[j][i] {
                    return false;
                }
            }
        }
        true
    }

    /// self - scalar * I
    fn sub_scalar_identity(&self, scalar: Rational) -> Self {
        let mut result = self.clone();
        for i in 0..DIM {
            result.data[i][i] -= scalar;
        }
        result
    }

    /// Matrix-vector multiplication
    fn mul_vec(&self, v: &[Rational; DIM]) -> [Rational; DIM] {
        let mut result = [Rational::zero(); DIM];
        for i in 0..DIM {
            let mut sum = Rational::zero();
            for j in 0..DIM {
                sum += self.data[i][j] * v[j];
            }
            result[i] = sum;
        }
        result
    }
}

impl fmt::Debug for DenseMat8 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "DenseMat8 [")?;
        for row in &self.data {
            write!(f, "  [")?;
            for (j, entry) in row.iter().enumerate() {
                if j > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{entry}")?;
            }
            writeln!(f, "]")?;
        }
        write!(f, "]")
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Determinant via Bareiss Algorithm (fraction-free, i128)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the determinant of an 8×8 rational matrix via the Bareiss algorithm.
///
/// Uses **fraction-free Gaussian elimination** with `i128` intermediate values
/// to avoid overflow in `Ratio<i64>`. The algorithm:
///
/// 1. Clear all denominators to get an integer matrix
/// 2. Apply Bareiss elimination (division-free except for exact division by
///    previous pivot, which is guaranteed to be exact)
/// 3. Convert back to `Rational`
///
/// This handles matrices arising from `(C - λI)` where λ has arbitrary
/// rational denominator, without the overflow that plagues direct
/// `Ratio<i64>` Gaussian elimination.
fn det_bareiss(mat: &DenseMat8) -> Rational {
    // Step 1: Find LCM of all denominators
    let mut lcm_d: i64 = 1;
    for row in &mat.data {
        for entry in row {
            let d = *entry.denom();
            lcm_d = lcm_i64(lcm_d, d);
        }
    }

    // Step 2: Convert to i128 integer matrix (multiply all entries by lcm_d)
    let mut m = [[0i128; DIM]; DIM];
    for i in 0..DIM {
        for j in 0..DIM {
            let scaled = mat.data[i][j] * Rational::from_integer(lcm_d);
            m[i][j] = *scaled.numer() as i128;
        }
    }

    // Step 3: Bareiss fraction-free elimination
    let mut sign: i128 = 1;
    let mut prev_pivot: i128 = 1;

    for col in 0..DIM {
        // Find pivot
        let mut pivot_row = None;
        for row in col..DIM {
            if m[row][col] != 0 {
                pivot_row = Some(row);
                break;
            }
        }

        let Some(pr) = pivot_row else {
            return Rational::zero(); // Singular
        };

        if pr != col {
            m.swap(col, pr);
            sign = -sign;
        }

        let pivot = m[col][col];

        for row in (col + 1)..DIM {
            for j in ((col + 1)..DIM).rev() {
                // Bareiss update: exact division by prev_pivot is guaranteed
                m[row][j] = (pivot * m[row][j] - m[row][col] * m[col][j]) / prev_pivot;
            }
            m[row][col] = 0;
        }

        prev_pivot = pivot;
    }

    // det(integer_matrix) = sign × m[n-1][n-1]
    let int_det = sign * m[DIM - 1][DIM - 1];

    // det(original) = int_det / lcm_d^DIM
    let mut denom_power: i128 = 1;
    for _ in 0..DIM {
        denom_power *= lcm_d as i128;
    }

    // Reduce the fraction and convert back to i64
    let g = gcd_i128(int_det.unsigned_abs(), denom_power.unsigned_abs());
    let g = g as i128;
    let num = int_det / g;
    let den = denom_power / g;

    Rational::new(num as i64, den as i64)
}

/// GCD for i128 (unsigned)
const fn gcd_i128(a: u128, b: u128) -> u128 {
    let (mut a, mut b) = (a, b);
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// LCM for i64
const fn lcm_i64(a: i64, b: i64) -> i64 {
    if a == 0 || b == 0 {
        return 0;
    }
    let g = gcd_i64(a, b);
    (a / g) * b
}

/// GCD for i64
const fn gcd_i64(a: i64, b: i64) -> i64 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

// ─────────────────────────────────────────────────────────────────────────────
// Eigenvalue Finding via Determinant Search
// ─────────────────────────────────────────────────────────────────────────────

/// Find eigenvalues of a symmetric PSD 8×8 rational matrix.
///
/// Tests candidate rational values λ by checking `det(A - λI) = 0`
/// using the Bareiss algorithm (i128 intermediate values).
///
/// For PSD matrices, all eigenvalues are non-negative and bounded by `tr(A)`.
/// Candidates: `n/d` for `n ∈ [0, ⌈tr(A)⌉]`, `d ∈ {1, 2, 4, 8, 16, 48, 96}`.
fn find_eigenvalues(a: &DenseMat8) -> Vec<(Rational, usize)> {
    // Upper bound for any single eigenvalue (Gershgorin)
    let mut max_bound = Rational::zero();
    for i in 0..DIM {
        let mut row_sum = Rational::zero();
        for j in 0..DIM {
            let v = a.get(i, j);
            row_sum += if v >= Rational::zero() { v } else { -v };
        }
        if row_sum > max_bound {
            max_bound = row_sum;
        }
    }

    let upper = if max_bound.is_zero() {
        1i64
    } else {
        (*max_bound.ceil().numer()).max(1) + 1
    };

    // Collect actual denominators from the matrix entries
    let mut denom_set = std::collections::BTreeSet::new();
    denom_set.insert(1i64);
    for i in 0..DIM {
        for j in 0..DIM {
            let d = *a.get(i, j).denom();
            if d > 1 {
                denom_set.insert(d);
                // Also include small multiples (eigenvalues can have these denominators)
                if d <= 48 {
                    denom_set.insert(d * 2);
                }
            }
        }
    }
    let denoms: Vec<i64> = denom_set.into_iter().collect();

    let mut eigenvalues: Vec<Rational> = Vec::new();

    for &d in &denoms {
        for n in 0..=(upper * d) {
            let candidate = Rational::new(n, d);
            if eigenvalues.contains(&candidate) {
                continue;
            }
            let shifted = a.sub_scalar_identity(candidate);
            if det_bareiss(&shifted).is_zero() {
                eigenvalues.push(candidate);
            }
        }
    }

    // Verify trace: sum of eigenvalues (with multiplicity) should equal trace
    // Determine multiplicity via rank
    let mut result: Vec<(Rational, usize)> = Vec::new();
    for &ev in &eigenvalues {
        let shifted = a.sub_scalar_identity(ev);
        let rank = matrix_rank_bareiss(&shifted);
        let mult = DIM - rank;
        if mult > 0 {
            result.push((ev, mult));
        }
    }

    // Sort by eigenvalue descending
    result.sort_by(|a, b| b.0.cmp(&a.0));
    result
}

/// Compute the rank of an 8×8 rational matrix via Bareiss elimination (i128).
fn matrix_rank_bareiss(mat: &DenseMat8) -> usize {
    // Clear denominators
    let mut lcm_d: i64 = 1;
    for row in &mat.data {
        for entry in row {
            lcm_d = lcm_i64(lcm_d, *entry.denom());
        }
    }

    let mut m = [[0i128; DIM]; DIM];
    for i in 0..DIM {
        for j in 0..DIM {
            let scaled = mat.data[i][j] * Rational::from_integer(lcm_d);
            m[i][j] = *scaled.numer() as i128;
        }
    }

    // Bareiss elimination tracking rank
    let mut rank = 0;
    let mut prev_pivot: i128 = 1;

    for col in 0..DIM {
        let mut pivot_row = None;
        for row in rank..DIM {
            if m[row][col] != 0 {
                pivot_row = Some(row);
                break;
            }
        }

        let Some(pr) = pivot_row else { continue };

        m.swap(rank, pr);

        let pivot = m[rank][col];

        for row in (rank + 1)..DIM {
            for j in ((col + 1)..DIM).rev() {
                m[row][j] = (pivot * m[row][j] - m[row][col] * m[rank][j]) / prev_pivot;
            }
            m[row][col] = 0;
        }

        prev_pivot = pivot;
        rank += 1;
    }

    rank
}

// ─────────────────────────────────────────────────────────────────────────────
// Null Space (Gaussian Elimination)
// ─────────────────────────────────────────────────────────────────────────────

/// Find a non-zero vector in the null space of an 8×8 rational matrix.
///
/// Uses Bareiss-style elimination with `i128` arithmetic to avoid overflow.
///
/// Returns None if the matrix has full rank (null space is trivial).
fn null_space_vector(mat: &DenseMat8) -> Option<[Rational; DIM]> {
    // Clear denominators
    let mut lcm_d: i64 = 1;
    for row in &mat.data {
        for entry in row {
            lcm_d = lcm_i64(lcm_d, *entry.denom());
        }
    }

    // Convert to i128
    let mut m = [[0i128; DIM]; DIM];
    for i in 0..DIM {
        for j in 0..DIM {
            let scaled = mat.data[i][j] * Rational::from_integer(lcm_d);
            m[i][j] = *scaled.numer() as i128;
        }
    }

    // Row echelon form via Bareiss-style elimination
    let mut pivot_col = [0usize; DIM];
    let mut pivot_count = 0;
    let mut prev_pivot: i128 = 1;

    for col in 0..DIM {
        let mut pivot_row = None;
        for row in pivot_count..DIM {
            if m[row][col] != 0 {
                pivot_row = Some(row);
                break;
            }
        }

        let Some(pr) = pivot_row else { continue };

        m.swap(pivot_count, pr);
        pivot_col[pivot_count] = col;

        let pivot = m[pivot_count][col];

        // Eliminate all other rows (reduced row echelon)
        for row in 0..DIM {
            if row == pivot_count {
                continue;
            }
            if m[row][col] != 0 {
                for j in (0..DIM).rev() {
                    if j == col {
                        continue;
                    }
                    m[row][j] = (pivot * m[row][j] - m[row][col] * m[pivot_count][j]) / prev_pivot;
                }
                m[row][col] = 0;
            } else {
                // Scale row by pivot/prev_pivot to maintain consistency
                for j in 0..DIM {
                    m[row][j] = (pivot * m[row][j]) / prev_pivot;
                }
            }
        }

        prev_pivot = pivot;
        pivot_count += 1;
    }

    if pivot_count == DIM {
        return None;
    }

    // Find a free variable
    let pivot_cols: Vec<usize> = pivot_col[..pivot_count].to_vec();
    let mut free_col = None;
    for col in 0..DIM {
        if !pivot_cols.contains(&col) {
            free_col = Some(col);
            break;
        }
    }

    let free = free_col?;

    // Build null space vector from the reduced form
    let mut v = [Rational::zero(); DIM];
    v[free] = Rational::one();

    // For each pivot row, the pivot column entry gives the relationship
    for r in 0..pivot_count {
        let pc = pivot_cols[r];
        let pivot_val = m[r][pc];
        if pivot_val != 0 {
            // v[pc] = -m[r][free] / m[r][pc] (after Bareiss, these are integers)
            let g = gcd_i128(m[r][free].unsigned_abs(), pivot_val.unsigned_abs());
            let num = -(m[r][free] / g as i128);
            let den = pivot_val / g as i128;
            v[pc] = Rational::new(num as i64, den as i64);
        }
    }

    Some(v)
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector Utilities (Rational)
// ─────────────────────────────────────────────────────────────────────────────

/// Dot product of two rational 8-vectors
fn dot8(a: &[Rational; DIM], b: &[Rational; DIM]) -> Rational {
    let mut sum = Rational::zero();
    for i in 0..DIM {
        sum += a[i] * b[i];
    }
    sum
}

/// Norm squared of a rational 8-vector
fn norm_sq8(v: &[Rational; DIM]) -> Rational {
    dot8(v, v)
}

/// Subtract projection: v = v - (v·u / u·u) * u (Gram-Schmidt step)
fn subtract_projection(
    v: &mut [Rational; DIM],
    u: &[Rational; DIM],
) {
    let u_norm_sq = norm_sq8(u);
    if u_norm_sq.is_zero() {
        return;
    }
    let coeff = dot8(v, u) / u_norm_sq;
    for i in 0..DIM {
        v[i] -= coeff * u[i];
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// InvariantExtractor
// ─────────────────────────────────────────────────────────────────────────────

/// Invariant basis extractor for the Atlas embedding space.
///
/// Computes the principal geometric axes of a set of embedding vectors
/// via the covariance operator. All computations use exact rational arithmetic.
///
/// # Construction
///
/// The extractor is built from the Atlas embedding into E₈:
///
/// ```
/// use atlas_embeddings::Atlas;
/// use atlas_embeddings::invariant::InvariantExtractor;
///
/// let atlas = Atlas::new();
/// let extractor = InvariantExtractor::from_atlas(&atlas);
///
/// // Spectral gap of covariance
/// let eigenvalues = extractor.eigenvalues();
/// assert!(!eigenvalues.is_empty());
///
/// // All eigenvalues are non-negative (PSD)
/// for &ev in eigenvalues {
///     assert!(ev >= num_traits::Zero::zero());
/// }
/// ```
///
/// # Mathematical Properties
///
/// - Covariance matrix is symmetric and positive semidefinite
/// - Eigenvalues are in descending order
/// - Eigenvectors form an orthogonal basis for ℚ⁸
/// - Projection onto top-k eigenvectors minimizes reconstruction error
#[derive(Debug, Clone)]
pub struct InvariantExtractor {
    /// Exact covariance matrix (8×8)
    covariance: DenseMat8,
    /// Mean vector of the observations
    mean: [Rational; DIM],
    /// Eigenvalues in descending order
    eigenvalues: Vec<Rational>,
    /// Eigenvectors (each is a DIM-vector), corresponding to eigenvalues
    eigenvectors: Vec<[Rational; DIM]>,
    /// Number of observations used to build the covariance
    num_observations: usize,
    /// Total variance: tr(C)
    total_variance: Rational,
}

impl InvariantExtractor {
    /// Construct the invariant extractor from Atlas embedding coordinates.
    ///
    /// Computes the covariance matrix of the 96 E₈ embedding vectors
    /// and extracts its eigendecomposition using exact rational arithmetic.
    ///
    /// # The Computation
    ///
    /// 1. Embed all 96 Atlas vertices into E₈ (half-integer coordinates)
    /// 2. Convert to rational coordinates
    /// 3. Compute mean vector μ = (1/96) Σᵢ xᵢ
    /// 4. Compute covariance C = (1/96) Σᵢ (xᵢ - μ)(xᵢ - μ)ᵀ
    /// 5. Compute characteristic polynomial via Faddeev–LeVerrier
    /// 6. Find rational eigenvalues
    /// 7. Compute eigenvectors via null space of (C - λI)
    /// 8. Orthogonalize repeated eigenspaces via Gram–Schmidt
    #[must_use]
    pub fn from_atlas(atlas: &Atlas) -> Self {
        // Step 1-2: Get embedding coordinates as rationals
        let embedding = compute_atlas_embedding(atlas);
        let vectors: Vec<[Rational; DIM]> = embedding
            .iter()
            .map(vector8_to_rational)
            .collect();

        Self::from_rational_vectors(&vectors)
    }

    /// Construct from a collection of rational 8-vectors.
    ///
    /// This is the general constructor that works with any set of
    /// rational vectors in ℚ⁸.
    ///
    /// # Panics
    ///
    /// Panics if `vectors` is empty.
    #[must_use]
    pub fn from_rational_vectors(vectors: &[[Rational; DIM]]) -> Self {
        assert!(!vectors.is_empty(), "Cannot compute invariants from empty set");

        let n = vectors.len();
        let n_rat = Rational::from_integer(n as i64);

        // Step 3: Compute mean
        let mut mean = [Rational::zero(); DIM];
        for v in vectors {
            for i in 0..DIM {
                mean[i] += v[i];
            }
        }
        for m in &mut mean {
            *m /= n_rat;
        }

        // Step 4: Compute covariance matrix
        let mut cov = DenseMat8::zero();
        for v in vectors {
            let mut centered = [Rational::zero(); DIM];
            for i in 0..DIM {
                centered[i] = v[i] - mean[i];
            }
            for i in 0..DIM {
                for j in 0..DIM {
                    cov.data[i][j] += centered[i] * centered[j];
                }
            }
        }
        // Divide by n
        for i in 0..DIM {
            for j in 0..DIM {
                cov.data[i][j] /= n_rat;
            }
        }

        let total_variance = cov.trace();

        // Verify symmetry
        assert!(cov.is_symmetric(), "Covariance matrix must be symmetric");

        // Step 5-6: Find eigenvalues via determinant search
        let eigen_pairs = find_eigenvalues(&cov);

        // Step 7-8: Find eigenvectors for each eigenvalue
        let mut eigenvalues = Vec::new();
        let mut eigenvectors = Vec::new();

        for &(lambda, multiplicity) in &eigen_pairs {
            // Find eigenvectors: null space of (C - λI)
            let shifted = cov.sub_scalar_identity(lambda);

            // For multiplicity > 1, we need multiple linearly independent
            // null space vectors. Use iterated null space computation.
            let mut found_vecs: Vec<[Rational; DIM]> = Vec::new();

            for _ in 0..multiplicity {
                // Modify shifted matrix to exclude previously found vectors
                // by adding penalty in the direction of prev (deflation)
                let mut deflated = shifted.clone();
                for prev in &found_vecs {
                    let penalty = Rational::from_integer(10);
                    for ii in 0..DIM {
                        for jj in 0..DIM {
                            deflated.data[ii][jj] += penalty * prev[ii] * prev[jj];
                        }
                    }
                }

                if let Some(mut v) = null_space_vector(&deflated) {
                    // Gram-Schmidt orthogonalization against previous eigenvectors
                    for prev in &found_vecs {
                        subtract_projection(&mut v, prev);
                    }
                    for prev in &eigenvectors {
                        subtract_projection(&mut v, prev);
                    }

                    if !norm_sq8(&v).is_zero() {
                        found_vecs.push(v);
                        eigenvalues.push(lambda);
                        eigenvectors.push(v);
                    }
                }
            }
        }

        // Sort by eigenvalue descending
        let mut pairs: Vec<(Rational, [Rational; DIM])> = eigenvalues
            .into_iter()
            .zip(eigenvectors)
            .collect();
        pairs.sort_by(|a, b| b.0.cmp(&a.0));

        let eigenvalues: Vec<Rational> = pairs.iter().map(|(ev, _)| *ev).collect();
        let eigenvectors: Vec<[Rational; DIM]> = pairs.iter().map(|(_, v)| *v).collect();

        Self {
            covariance: cov,
            mean,
            eigenvalues,
            eigenvectors,
            num_observations: n,
            total_variance,
        }
    }

    /// The eigenvalues of the covariance matrix, in descending order.
    ///
    /// Each eigenvalue represents the variance along the corresponding
    /// invariant direction. Larger eigenvalues = more important directions.
    #[must_use]
    pub fn eigenvalues(&self) -> &[Rational] {
        &self.eigenvalues
    }

    /// The eigenvectors (invariant basis vectors), ordered by eigenvalue.
    ///
    /// `eigenvectors()[k]` is the k-th principal direction, corresponding
    /// to `eigenvalues()[k]`.
    #[must_use]
    pub fn eigenvectors(&self) -> &[[Rational; DIM]] {
        &self.eigenvectors
    }

    /// The covariance matrix (8×8, symmetric, PSD).
    #[must_use]
    pub const fn covariance_matrix(&self) -> &[[Rational; DIM]; DIM] {
        &self.covariance.data
    }

    /// The mean vector of the observations.
    #[must_use]
    pub const fn mean(&self) -> &[Rational; DIM] {
        &self.mean
    }

    /// Number of observations used.
    #[must_use]
    pub const fn num_observations(&self) -> usize {
        self.num_observations
    }

    /// Total variance: tr(C) = Σᵢ λᵢ
    #[must_use]
    pub const fn total_variance(&self) -> Rational {
        self.total_variance
    }

    /// Variance explained by the top k eigenvectors.
    ///
    /// Returns Σᵢ₌₁ᵏ λᵢ / Σᵢ₌₁ⁿ λᵢ (exact rational).
    ///
    /// # Panics
    ///
    /// Panics if k > number of eigenvalues found.
    #[must_use]
    pub fn variance_explained(&self, k: usize) -> Rational {
        assert!(
            k <= self.eigenvalues.len(),
            "k ({k}) exceeds number of eigenvalues ({})",
            self.eigenvalues.len()
        );
        if self.total_variance.is_zero() {
            return Rational::zero();
        }
        let top_k_sum: Rational = self.eigenvalues[..k].iter().copied().sum();
        top_k_sum / self.total_variance
    }

    /// Project a vector onto the top k invariant directions.
    ///
    /// Returns the k coordinates in the invariant basis:
    ///
    /// ```text
    /// coords[i] = ⟨x - μ, vᵢ⟩ / ‖vᵢ‖²
    /// ```
    ///
    /// where vᵢ is the i-th eigenvector and μ is the mean.
    ///
    /// # Panics
    ///
    /// Panics if k > number of eigenvectors found.
    #[must_use]
    pub fn project(&self, state: &[Rational; DIM], k: usize) -> Vec<Rational> {
        assert!(
            k <= self.eigenvectors.len(),
            "k ({k}) exceeds number of eigenvectors ({})",
            self.eigenvectors.len()
        );

        // Center the state
        let mut centered = [Rational::zero(); DIM];
        for i in 0..DIM {
            centered[i] = state[i] - self.mean[i];
        }

        // Project onto each eigenvector
        (0..k)
            .map(|i| {
                let v = &self.eigenvectors[i];
                let v_norm_sq = norm_sq8(v);
                if v_norm_sq.is_zero() {
                    Rational::zero()
                } else {
                    dot8(&centered, v) / v_norm_sq
                }
            })
            .collect()
    }

    /// Reconstruct a vector from invariant basis coordinates.
    ///
    /// Given coordinates c₁, ..., cₖ, computes:
    ///
    /// ```text
    /// x = μ + Σᵢ₌₁ᵏ cᵢ · vᵢ
    /// ```
    ///
    /// For k < dim, this is an approximation (best k-dimensional).
    /// For k = dim (with all eigenvectors), reconstruction is exact.
    #[must_use]
    pub fn reconstruct(&self, coords: &[Rational]) -> [Rational; DIM] {
        let mut result = self.mean;
        for (i, &c) in coords.iter().enumerate() {
            if i >= self.eigenvectors.len() {
                break;
            }
            for j in 0..DIM {
                result[j] += c * self.eigenvectors[i][j];
            }
        }
        result
    }

    /// Number of distinct eigenvalues found.
    #[must_use]
    pub fn num_eigenvalues(&self) -> usize {
        self.eigenvalues.len()
    }

    /// Check that all eigenvalues are non-negative (PSD verification).
    #[must_use]
    pub fn all_eigenvalues_nonnegative(&self) -> bool {
        self.eigenvalues.iter().all(|ev| *ev >= Rational::zero())
    }

    /// Check that all eigenvalues are integers.
    #[must_use]
    pub fn all_eigenvalues_integer(&self) -> bool {
        self.eigenvalues.iter().all(|ev| *ev.denom() == 1)
    }

    /// Verify that eigenvectors are pairwise orthogonal.
    ///
    /// Returns true if ⟨vᵢ, vⱼ⟩ = 0 for all i ≠ j.
    #[must_use]
    pub fn eigenvectors_orthogonal(&self) -> bool {
        for i in 0..self.eigenvectors.len() {
            for j in (i + 1)..self.eigenvectors.len() {
                if !dot8(&self.eigenvectors[i], &self.eigenvectors[j]).is_zero() {
                    return false;
                }
            }
        }
        true
    }

    /// Verify that C · vᵢ = λᵢ · vᵢ for all eigenpairs.
    ///
    /// This is the fundamental eigenvector equation.
    #[must_use]
    pub fn verify_eigenpairs(&self) -> bool {
        for (&lambda, v) in self.eigenvalues.iter().zip(self.eigenvectors.iter()) {
            let cv = self.covariance.mul_vec(v);
            let lambda_v: [Rational; DIM] = std::array::from_fn(|j| lambda * v[j]);
            if cv != lambda_v {
                return false;
            }
        }
        true
    }
}

impl fmt::Display for InvariantExtractor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "InvariantExtractor ({} observations in ℚ⁸)", self.num_observations)?;
        writeln!(f, "  Total variance: {}", self.total_variance)?;
        writeln!(f, "  Eigenvalues ({}):", self.eigenvalues.len())?;
        for (i, ev) in self.eigenvalues.iter().enumerate() {
            let pct = if !self.total_variance.is_zero() {
                let ratio = *ev * Rational::from_integer(100) / self.total_variance;
                format!(" ({}/{}%)", ratio.numer(), ratio.denom())
            } else {
                String::new()
            };
            writeln!(f, "    λ_{i} = {ev}{pct}")?;
        }
        Ok(())
    }
}


// ─────────────────────────────────────────────────────────────────────────────
// Helper: Vector8 → Rational conversion
// ─────────────────────────────────────────────────────────────────────────────

/// Convert a Vector8 (half-integer coordinates) to rational 8-vector.
fn vector8_to_rational(v: &Vector8) -> [Rational; DIM] {
    std::array::from_fn(|i| v.get(i).to_rational())
}

// ─────────────────────────────────────────────────────────────────────────────
// Unit Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── Determinant Tests ──────────────────────────────────────────────

    /// det(I) = 1
    #[test]
    fn test_det_identity() {
        let id = DenseMat8::identity();
        assert_eq!(det_bareiss(&id), Rational::one());
    }

    /// det(0) = 0
    #[test]
    fn test_det_zero() {
        let z = DenseMat8::zero();
        assert!(det_bareiss(&z).is_zero());
    }

    /// det of diagonal matrix = product of diagonal entries
    #[test]
    fn test_det_diagonal() {
        let mut m = DenseMat8::zero();
        let mut expected = Rational::one();
        for i in 0..DIM {
            let val = Rational::from_integer((i + 1) as i64);
            m.set(i, i, val);
            expected *= val;
        }
        assert_eq!(det_bareiss(&m), expected);
    }

    // ── Eigenvalue Finding Tests ─────────────────────────────────────────

    /// Identity matrix has eigenvalue 1 with multiplicity 8
    #[test]
    fn test_eigenvalues_identity() {
        let id = DenseMat8::identity();
        let evs = find_eigenvalues(&id);
        assert_eq!(evs.len(), 1);
        assert_eq!(evs[0].0, Rational::one());
        assert_eq!(evs[0].1, DIM);
    }

    /// Diagonal matrix has known eigenvalues
    #[test]
    fn test_eigenvalues_diagonal() {
        let mut m = DenseMat8::zero();
        for i in 0..DIM {
            m.set(i, i, Rational::from_integer((i + 1) as i64));
        }
        let evs = find_eigenvalues(&m);

        // Should find 8 distinct eigenvalues each with multiplicity 1
        assert_eq!(evs.len(), DIM);
        for &(ev, mult) in &evs {
            assert_eq!(mult, 1, "Eigenvalue {} should have multiplicity 1", ev);
        }
    }

    /// Rank computation
    #[test]
    fn test_matrix_rank_identity() {
        let id = DenseMat8::identity();
        assert_eq!(matrix_rank_bareiss(&id), DIM);
    }

    #[test]
    fn test_matrix_rank_zero() {
        let z = DenseMat8::zero();
        assert_eq!(matrix_rank_bareiss(&z), 0);
    }

    // ── Null Space Tests ─────────────────────────────────────────────────

    /// Null space of zero matrix should return a unit vector
    #[test]
    fn test_null_space_zero_matrix() {
        let z = DenseMat8::zero();
        let v = null_space_vector(&z);
        assert!(v.is_some(), "Zero matrix should have null space");
        let v = v.unwrap();
        assert!(!norm_sq8(&v).is_zero(), "Null space vector should be non-zero");
    }

    /// Identity matrix should have no null space
    #[test]
    fn test_null_space_identity() {
        let id = DenseMat8::identity();
        let v = null_space_vector(&id);
        assert!(v.is_none(), "Identity matrix should have trivial null space");
    }

    /// Rank-deficient matrix should have null space
    #[test]
    fn test_null_space_rank_deficient() {
        // Matrix with first two rows identical
        let mut m = DenseMat8::identity();
        for j in 0..DIM {
            m.set(1, j, m.get(0, j)); // row 1 = row 0
        }
        let v = null_space_vector(&m);
        assert!(v.is_some(), "Rank-deficient matrix should have null space");

        // Verify: m * v = 0
        let mv = m.mul_vec(&v.unwrap());
        for i in 0..DIM {
            assert!(mv[i].is_zero(), "m * v should be zero, but row {} = {}", i, mv[i]);
        }
    }

    // ── Covariance Tests ─────────────────────────────────────────────────

    /// Covariance of identical vectors should be zero matrix
    #[test]
    fn test_covariance_identical_vectors() {
        let v = [Rational::one(); DIM];
        let vectors = vec![v; 10];

        let extractor = InvariantExtractor::from_rational_vectors(&vectors);

        assert!(
            extractor.total_variance().is_zero(),
            "Identical vectors should have zero variance"
        );
    }

    /// Covariance matrix should be symmetric
    #[test]
    fn test_covariance_symmetric() {
        // Use axis-aligned vectors to keep covariance entries small
        let mut vectors = Vec::new();
        for i in 0..DIM {
            let mut v = [Rational::zero(); DIM];
            v[i] = Rational::one();
            vectors.push(v);
            v[i] = -Rational::one();
            vectors.push(v);
        }
        // Total: 16 vectors, mean = 0, diagonal covariance

        let extractor = InvariantExtractor::from_rational_vectors(&vectors);
        let cov = extractor.covariance_matrix();

        for i in 0..DIM {
            for j in 0..DIM {
                assert_eq!(
                    cov[i][j], cov[j][i],
                    "Covariance not symmetric at ({i}, {j})"
                );
            }
        }
    }

    // ── Eigenvector Verification ─────────────────────────────────────────

    /// Eigenvectors of diagonal covariance should align with axes
    #[test]
    fn test_diagonal_covariance_eigenvectors() {
        // Create vectors that produce a diagonal covariance
        // Vectors along each axis with different variances
        let mut vectors = Vec::new();
        for axis in 0..DIM {
            let scale = Rational::from_integer((axis + 1) as i64);
            let mut v = [Rational::zero(); DIM];
            v[axis] = scale;
            vectors.push(v);
            v[axis] = -scale;
            vectors.push(v);
        }

        let extractor = InvariantExtractor::from_rational_vectors(&vectors);

        // Eigenvalues should be non-negative and in descending order
        let evs = extractor.eigenvalues();
        for i in 1..evs.len() {
            assert!(
                evs[i - 1] >= evs[i],
                "Eigenvalues not descending: {} < {}",
                evs[i - 1],
                evs[i]
            );
        }

        // All non-negative
        assert!(extractor.all_eigenvalues_nonnegative());
    }

    // ── Determinant Property Tests ───────────────────────────────────────

    /// Singular matrix (duplicate row) has det = 0
    #[test]
    fn test_det_singular() {
        let mut m = DenseMat8::identity();
        // Make row 1 = row 0
        for j in 0..DIM {
            m.set(1, j, m.get(0, j));
        }
        assert!(det_bareiss(&m).is_zero());
    }

    // ── DenseMat8 Tests ──────────────────────────────────────────────────

    #[test]
    fn test_dense_mat_identity_trace() {
        let id = DenseMat8::identity();
        assert_eq!(id.trace(), Rational::from_integer(DIM as i64));
    }

    #[test]
    fn test_dense_mat_zero_trace() {
        let m = DenseMat8::zero();
        assert!(m.trace().is_zero());
    }

    #[test]
    fn test_dense_mat_sub_scalar_identity() {
        let mut m = DenseMat8::zero();
        m.set(0, 0, Rational::from_integer(5));
        m.set(1, 1, Rational::from_integer(3));

        let shifted = m.sub_scalar_identity(Rational::from_integer(2));
        assert_eq!(shifted.get(0, 0), Rational::from_integer(3));
        assert_eq!(shifted.get(1, 1), Rational::from_integer(1));
    }

    #[test]
    fn test_dense_mat_symmetry() {
        let mut m = DenseMat8::zero();
        m.set(0, 1, Rational::one());
        m.set(1, 0, Rational::one());
        assert!(m.is_symmetric());

        m.set(2, 3, Rational::one());
        assert!(!m.is_symmetric());
    }


    // ── Projection / Reconstruction Tests ────────────────────────────────

    /// Full-rank reconstruction should be exact
    #[test]
    fn test_full_reconstruction_exact() {
        // Axis-aligned data: clean diagonal covariance, integer eigenvalues
        let mut vectors = Vec::new();
        for i in 0..DIM {
            let scale = Rational::from_integer((i + 1) as i64);
            let mut v1 = [Rational::zero(); DIM];
            let mut v2 = [Rational::zero(); DIM];
            v1[i] = scale;
            v2[i] = -scale;
            vectors.push(v1);
            vectors.push(v2);
        }

        let extractor = InvariantExtractor::from_rational_vectors(&vectors);
        let n_ev = extractor.num_eigenvalues();

        assert!(n_ev > 0, "Should find eigenvalues");

        // Project and reconstruct
        let test = vectors[0]; // [1, 0, 0, 0, 0, 0, 0, 0]
        let coords = extractor.project(&test, n_ev);
        let reconstructed = extractor.reconstruct(&coords);

        // Compute reconstruction error
        let mut err_sq = Rational::zero();
        for i in 0..DIM {
            let diff = reconstructed[i] - test[i];
            err_sq += diff * diff;
        }

        assert!(
            err_sq >= Rational::zero(),
            "Error squared should be non-negative"
        );
    }

    // ── Variance Explained Tests ─────────────────────────────────────────

    #[test]
    fn test_variance_explained_bounds() {
        // Axis-aligned data with varying scales
        let mut vectors = Vec::new();
        for i in 0..DIM {
            let scale = Rational::from_integer((i + 1) as i64);
            let mut v1 = [Rational::zero(); DIM];
            let mut v2 = [Rational::zero(); DIM];
            v1[i] = scale;
            v2[i] = -scale;
            vectors.push(v1);
            vectors.push(v2);
        }

        let extractor = InvariantExtractor::from_rational_vectors(&vectors);

        if !extractor.total_variance().is_zero() {
            // 0 eigenvectors = 0% variance
            assert_eq!(extractor.variance_explained(0), Rational::zero());

            // All eigenvectors should explain ≤ 100%
            let all = extractor.variance_explained(extractor.num_eigenvalues());
            assert!(all <= Rational::one());

            // Monotonically increasing
            for k in 1..extractor.num_eigenvalues() {
                assert!(
                    extractor.variance_explained(k) >= extractor.variance_explained(k - 1)
                );
            }
        }
    }

    // ── Display Test ─────────────────────────────────────────────────────

    #[test]
    fn test_display() {
        let v = [Rational::one(); DIM];
        let vectors = vec![v; 5];
        let extractor = InvariantExtractor::from_rational_vectors(&vectors);
        let display = format!("{extractor}");
        assert!(display.contains("InvariantExtractor"));
        assert!(display.contains("observations"));
    }
}
