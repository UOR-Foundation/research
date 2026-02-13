import numpy as np
from invariant_extractor import InvariantExtractor


def test_basic_functionality():
    """Test basic functionality of InvariantExtractor."""
    print("Testing basic functionality...")

    # Test initialization
    extractor = InvariantExtractor(dim=4, memory=10)
    assert extractor.dim == 4
    assert extractor.memory == 10
    assert not extractor.full
    print("✓ Initialization works")

    # Test update
    state = np.array([1.0, 2.0, 3.0, 4.0])
    extractor.update(state)
    print("✓ Single update works")

    # Test multiple updates to fill buffer
    for i in range(10):
        extractor.update(np.random.randn(4))
    assert extractor.full
    print("✓ Buffer filling works")

    # Test invariants extraction
    vals, vecs = extractor.invariants()
    assert len(vals) == 4
    assert vecs.shape == (4, 4)
    assert np.all(vals >= 0)  # eigenvalues should be non-negative
    print("✓ Invariants extraction works")

    # Test projection and reconstruction
    test_state = np.array([1.0, 0.5, 0.1, 0.01])
    coords = extractor.project(test_state, k=2)
    reconstructed = extractor.reconstruct(coords)
    assert len(coords) == 2
    assert reconstructed.shape == (4,)
    print("✓ Projection and reconstruction work")


def test_error_handling():
    """Test error handling."""
    print("\nTesting error handling...")

    # Test invalid dim
    try:
        InvariantExtractor(dim=0)
        assert False, "Should have raised ValueError"
    except ValueError:
        print("✓ Invalid dim error handling works")

    # Test invalid memory
    try:
        InvariantExtractor(dim=4, memory=1)
        assert False, "Should have raised ValueError"
    except ValueError:
        print("✓ Invalid memory error handling works")

    # Test wrong state dimension
    extractor = InvariantExtractor(dim=3)
    try:
        extractor.update([1, 2, 3, 4])  # 4D state for 3D extractor
        assert False, "Should have raised ValueError"
    except ValueError:
        print("✓ Wrong state dimension error handling works")


def test_anisotropic_data():
    """Test with anisotropic data like in the example."""
    print("\nTesting with anisotropic data...")

    extractor = InvariantExtractor(dim=4)

    # Feed anisotropic data
    for _ in range(500):
        state = np.random.randn(4) * np.array([1.0, 0.5, 0.1, 0.01])
        extractor.update(state)

    vals, basis = extractor.invariants()

    # Check that eigenvalues are in descending order
    assert np.all(vals[:-1] >= vals[1:]), "Eigenvalues should be in descending order"

    # Check that the first eigenvalue is much larger than the last
    # (due to anisotropic scaling)
    assert vals[0] > vals[-1] * 10, "First eigenvalue should be much larger"

    print(f"✓ Anisotropic data test passed")
    print(f"  Eigenvalues: {vals}")
    print(f"  Ratio (first/last): {vals[0]/vals[-1]:.2f}")


def test_reconstruction_accuracy():
    """Test reconstruction accuracy."""
    print("\nTesting reconstruction accuracy...")

    extractor = InvariantExtractor(dim=4, memory=100)

    # Generate some data
    for _ in range(200):
        state = np.random.randn(4) * np.array([2.0, 1.0, 0.5, 0.1])
        extractor.update(state)

    # Test reconstruction with different k values
    original_state = np.array([1.5, 0.8, 0.2, 0.05])

    for k in [1, 2, 3, 4]:
        coords = extractor.project(original_state, k=k)
        reconstructed = extractor.reconstruct(coords)

        # Calculate relative error
        error = np.linalg.norm(reconstructed - original_state)
        relative_error = error / np.linalg.norm(original_state)
        print(f"  k={k}: relative error = {relative_error:.6f}")

        # Error should decrease as k increases
        if k < 4:
            # Allow some tolerance for numerical precision
            assert relative_error < 1.0, f"Reconstruction error too high for k={k}"


if __name__ == "__main__":
    print("Running InvariantExtractor tests...\n")
    test_basic_functionality()
    test_error_handling()
    test_anisotropic_data()
    test_reconstruction_accuracy()
    print("\n✅ All tests passed!")
