"""
Invariant basis extractor for atlas-embeddings.

Extracts the stable invariant basis from a stream of evolving embedding
states using the covariance operator. The dominant eigenmodes define the
intrinsic geometric axes of the system, separating meaningful structure
from transient components.

This acts as a coordinate-locking layer: instead of operating in an
arbitrary embedding space, the system aligns itself to its own intrinsic
geometry.
"""

import numpy as np
from numpy.linalg import eigh


class InvariantExtractor:
    """Extract invariant basis from evolving embedding states.

    Maintains a rolling buffer of observed states and computes the
    covariance operator. The eigenvectors of this operator define the
    invariant basis — the directions where the system's geometry lives.

    Parameters
    ----------
    dim : int
        Dimensionality of the state vectors (must be >= 1).
    memory : int, optional
        Number of recent states to retain in the rolling buffer.
        Must be >= 2 to compute a meaningful covariance. Default is 50.

    Examples
    --------
    >>> extractor = InvariantExtractor(dim=4, memory=100)
    >>> for state in state_stream:
    ...     extractor.update(state)
    >>> eigenvalues, basis = extractor.invariants()
    >>> coords = extractor.project(state, k=2)
    >>> reconstructed = extractor.reconstruct(coords)
    """

    def __init__(self, dim: int, memory: int = 50):
        if dim < 1:
            raise ValueError(f"dim must be >= 1, got {dim}")
        if memory < 2:
            raise ValueError(f"memory must be >= 2, got {memory}")

        self.dim = dim
        self.memory = memory
        self._buffer = np.zeros((memory, dim))
        self._count = 0
        self._basis = None
        self._eigenvalues = None

    @property
    def full(self) -> bool:
        """Whether the rolling buffer has been completely filled."""
        return self._count >= self.memory

    def update(self, state) -> None:
        """Add a state observation to the rolling buffer.

        Parameters
        ----------
        state : array_like
            State vector of length `dim`.

        Raises
        ------
        ValueError
            If `state` has wrong dimensionality.
        """
        state = np.asarray(state, dtype=float)
        if state.shape != (self.dim,):
            raise ValueError(
                f"Expected state of shape ({self.dim},), got {state.shape}"
            )

        idx = self._count % self.memory
        self._buffer[idx] = state
        self._count += 1
        # Invalidate cached decomposition
        self._basis = None
        self._eigenvalues = None

    def _compute(self) -> None:
        """Compute eigendecomposition of the covariance operator."""
        n = min(self._count, self.memory)
        data = self._buffer[:n]
        mean = data.mean(axis=0)
        centered = data - mean
        cov = (centered.T @ centered) / n
        eigenvalues, eigenvectors = eigh(cov)
        # Sort descending
        order = np.argsort(eigenvalues)[::-1]
        self._eigenvalues = eigenvalues[order]
        self._basis = eigenvectors[:, order]

    def invariants(self):
        """Return the invariant basis and associated eigenvalues.

        Returns
        -------
        eigenvalues : ndarray, shape (dim,)
            Eigenvalues in descending order (variance along each axis).
        basis : ndarray, shape (dim, dim)
            Columns are the invariant basis vectors, ordered by eigenvalue.
        """
        if self._eigenvalues is None:
            self._compute()
        return self._eigenvalues.copy(), self._basis.copy()

    def project(self, state, k: int = None):
        """Project a state onto the top-k invariant directions.

        Parameters
        ----------
        state : array_like
            State vector of length `dim`.
        k : int, optional
            Number of dominant directions to keep. Defaults to `dim`.

        Returns
        -------
        coords : ndarray, shape (k,)
            Coordinates in the invariant basis (top-k components).
        """
        if k is None:
            k = self.dim
        if self._basis is None:
            self._compute()
        state = np.asarray(state, dtype=float)
        return (self._basis[:, :k].T @ state)

    def reconstruct(self, coords):
        """Reconstruct a state from invariant-basis coordinates.

        Parameters
        ----------
        coords : array_like
            Coordinates from `project()`.

        Returns
        -------
        state : ndarray, shape (dim,)
            Reconstructed state vector.
        """
        coords = np.asarray(coords, dtype=float)
        k = len(coords)
        if self._basis is None:
            self._compute()
        return self._basis[:, :k] @ coords
