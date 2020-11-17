import numpy as np
import numpy.testing as npt
import pytest
import prism


def test_invalid_field():
    "Check if passing an invalid field raises an error"
    model = np.array([[-100, 100, -100, 100, 100, 200]])
    density = np.array([1000])
    coordinates = np.array([[0], [0], [0]])
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="invalid field")


def test_invalid_prism_boundaries():
    "Check if passing an invalid prism boundaries raises an error"
    density = np.array([1000])
    coordinates = np.array([[0], [0], [0]])
    field = "potential"
    # wrong x boundaries
    model = np.array([[100, -100, -100, 100, 100, 200]])
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="g_x")
    # wrong y boundaries
    model = np.array([[-100, 100, 100, -100, 100, 200]])
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="g_z")
    # wrong z boundaries
    model = np.array([[-100, 100, -100, 100, 200, 100]])
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="potential")


def test_invalid_prism():
    "Check if passing an invalid prism raises an error"
    density = np.array([1000])
    coordinates = np.array([[0], [0], [0]])
    field = "potential"
    # shape (1,)
    model = np.array([100, -100, -100, 100, 100, 200])
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="potential")
    # shape (2,4)
    model = np.empty((2, 4))
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="g_x")
    # shape (1,4)
    model = np.empty((1, 5))
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="g_z")


def test_invalid_coordinates():
    "Check if passing an invalid coordinates raises an error"
    model = np.array([[-100, 100, -100, 100, 100, 200]])
    density = np.array([1000])
    # shape (1,)
    coordinates = np.array([0, 0, 0])
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="g_z")
    # shape (4,3)
    coordinates = np.zeros((4,3))
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="g_z")


def test_invalid_density():
    "Check if density with shape[0] != 3 raises an error"
    model = np.array([[-100, 100, -100, 100, 100, 200]])
    density = np.array([1000])
    coordinates = np.array([0, 0, 0])
    with pytest.raises(ValueError):
        prism.gravitational(coordinates, model, density, field="g_z")


def test_safe_atan2():
    "Test the safe_atan2 function"
    # Test safe_atan2 for one point per quadrant
    # First quadrant
    x, y = 1, 1
    npt.assert_allclose(prism.safe_atan2(y, x), np.pi / 4)
    # Second quadrant
    x, y = -1, 1
    npt.assert_allclose(prism.safe_atan2(y, x), -np.pi / 4)
    # Third quadrant
    x, y = -1, -1
    npt.assert_allclose(prism.safe_atan2(y, x), np.pi / 4)
    # Forth quadrant
    x, y = 1, -1
    npt.assert_allclose(prism.safe_atan2(y, x), -np.pi / 4)
    # Test safe_atan2 if the denominator is equal to zero
    npt.assert_allclose(prism.safe_atan2(1, 0), np.pi / 2)
    npt.assert_allclose(prism.safe_atan2(-1, 0), -np.pi / 2)
    # Test safe_atan2 if both numerator and denominator are equal to zero
    npt.assert_allclose(prism.safe_atan2(0, 0), 0)


def test_safe_log():
    "Test the safe_log function"
    # Check if safe_log function satisfies safe_log(0) == 0
    npt.assert_allclose(prism.safe_log(0), 0)
    # Check if safe_log behaves like the natural logarithm in case that x != 0
    x = np.linspace(1, 100, 101)
    for x_i in x:
        npt.assert_allclose(prism.safe_log(x_i), np.log(x_i))


def test_field_decreases_with_distance():
    "Check if field decreases with distance"
    model = np.array([[-100, 100, -100, 100, 100, 200]])
    density = np.array([1000])
    close = np.array([[0], [20], [0]])
    far = np.array([[0], [20], [-100]])
    # potentia
    potential_close = prism.gravitational(close, model, density, field="potential")
    potential_far = prism.gravitational(far, model, density, field="potential")
    # gz
    gz_close = prism.gravitational(close, model, density, field="g_z")
    gz_far = prism.gravitational(far, model, density, field="g_z")
    # gx
    gx_close = prism.gravitational(close, model, density, field="g_x")
    gx_far = prism.gravitational(far, model, density, field="g_x")
    diffs = np.array([np.abs(potential_far) < np.abs(potential_close),
                      np.abs(gz_far) < np.abs(gz_close),
                      np.abs(gx_far) < np.abs(gx_close)])
    npt.assert_allclose(diffs, np.ones((3,1), dtype=bool))
