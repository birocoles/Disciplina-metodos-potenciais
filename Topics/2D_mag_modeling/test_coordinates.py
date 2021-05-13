import numpy as np
import coordinates as coord
from numpy.testing import assert_almost_equal as aae
from pytest import raises
from functools import reduce


def test_rotation_NEU_known_values():
    'verify results obtained for known input'
    latitude = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    longitude = np.array([0, 0, 0, 90, 90, 90, 180, 180])

    sinlat = np.array([0, -(np.sqrt(6) - np.sqrt(2))/4,
                      (np.sqrt(2 - np.sqrt(2)))/2, -0.5,
                      np.sqrt(2)/2, -np.sqrt(3)/2,
                      (np.sqrt(6) + np.sqrt(2))/4, -1])
    coslat = np.array([1, (np.sqrt(6) + np.sqrt(2))/4,
                      (np.sqrt(2 + np.sqrt(2)))/2, np.sqrt(3)/2,
                      np.sqrt(2)/2, 0.5, (np.sqrt(6) - np.sqrt(2))/4, 0])
    sinlon = np.array([0, 0, 0, 1, 1, 1, 0, 0])
    coslon = np.array([1, 1, 1, 0, 0, 0, -1, -1])

    R11_true = -sinlat*coslon
    R12_true = -sinlon
    R13_true = coslat*coslon
    R21_true = -sinlat*sinlon
    R22_true = coslon
    R23_true = coslat*sinlon
    R31_true = coslat
    #R32_true = 0.
    R33_true = sinlat

    R_true = np.vstack([R11_true, R12_true, R13_true,
                        R21_true, R22_true, R23_true,
                        R31_true, R33_true]).T

    R = coord.rotation_NEU(latitude, longitude)

    aae(R_true, R, decimal=15)


def test_rotation_NEU_known_values_floats():
    'verify results obtained for known input'
    latitude = -15
    longitude = 0

    sinlat = -(np.sqrt(6) - np.sqrt(2))/4
    coslat = (np.sqrt(6) + np.sqrt(2))/4
    sinlon = 0
    coslon = 1

    R11_true = -sinlat*coslon
    R12_true = -sinlon
    R13_true = coslat*coslon
    R21_true = -sinlat*sinlon
    R22_true = coslon
    R23_true = coslat*sinlon
    R31_true = coslat
    #R32_true = 0.
    R33_true = sinlat

    R_true = np.array([[R11_true, R12_true, R13_true,
                        R21_true, R22_true, R23_true,
                        R31_true, R33_true]])

    R = coord.rotation_NEU(latitude, longitude)

    aae(R_true, R, decimal=15)


def test_rotation_NEU_orthogonality():
    'Rotation matrix must be mutually orthogonal'
    latitude = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    longitude = np.array([-17, 0, 30, 9, 90, 23, 180, 1])
    R = coord.rotation_NEU(latitude, longitude)
    for Ri in R:
        M = np.array([[Ri[0], Ri[1], Ri[2]],
                      [Ri[3], Ri[4], Ri[5]],
                      [Ri[6],     0, Ri[7]]])
        aae(np.dot(M.T, M), np.identity(3), decimal=15)
        aae(np.dot(M, M.T), np.identity(3), decimal=15)


def test_rotation_NEU_lines_bad_arguments():
    'latitude and longitude with different number of elements'
    latitude = np.ones(100)
    longitude = np.zeros(34)
    raises(AssertionError, coord.rotation_NEU, latitude, longitude)


def test_rotation_NED_orthogonality():
    'Rotation matrix must be mutually orthogonal'
    latitude = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    longitude = np.array([-17, 0, 30, 9, 90, 23, 180, 1])
    R = coord.rotation_NED(latitude, longitude)
    for Ri in R:
        M = np.array([[Ri[0], Ri[1], Ri[2]],
                      [Ri[3], Ri[4], Ri[5]],
                      [Ri[6],     0, Ri[7]]])
        aae(np.dot(M.T, M), np.identity(3), decimal=15)
        aae(np.dot(M, M.T), np.identity(3), decimal=15)


def test_rotation_NED_lines_bad_arguments():
    'latitude and longitude with different number of elements'
    latitude = np.ones(100)
    longitude = np.zeros(34)
    raises(AssertionError, coord.rotation_NED, latitude, longitude)


def test_rotation_NEU_versus_NED():
    'Two first columns must be the same and the last must be opposite'
    latitude = np.array([0, -15, 22.5, -30, 45, -60, 75, -90])
    longitude = np.array([-17, 0, 30, 9, 90, 23, 180, 1])
    R_NEU = coord.rotation_NEU(latitude, longitude)
    R_NED = coord.rotation_NED(latitude, longitude)
    for R1, R2 in zip(R_NEU, R_NED):
        aae(R1[0], R2[0], decimal=15)
        aae(R1[1], R2[1], decimal=15)
        aae(R1[2], -R2[2], decimal=15)
        aae(R1[3], R2[3], decimal=15)
        aae(R1[4], R2[4], decimal=15)
        aae(R1[5], -R2[5], decimal=15)
        aae(R1[6], R2[6], decimal=15)
        aae(R1[7], -R2[7], decimal=15)


def test_unit_vector_TCCS_bad_arguments():
    'inclination and declination with different number of elements'
    inclination = np.empty(12)
    declination = np.empty(10)
    raises(AssertionError, coord.unit_vector_TCCS, inclination, declination)


def test_unit_vector_TCCS_returns_unit_vector():
    'the computed vectos must be unitary'
    inclination, declination = np.meshgrid(np.linspace(-90, 90, 3),
                                           np.linspace(0, 180, 3))
    inclination = np.ravel(inclination)
    declination = np.ravel(declination)
    vx, vy, vz = coord.unit_vector_TCCS(inclination, declination)
    norm = np.sqrt(vx*vx + vy*vy + vz*vz)
    aae(norm, np.ones_like(norm), decimal=15)


def test_R1_R2_R3_orthonal():
    'Rotation matrices must be orthogonal'
    A = coord.R1(-19)
    B = coord.R2(34.71)
    C = coord.R3(28)

    aae(np.dot(A, A.T), np.dot(A.T, A), decimal=15)
    aae(np.dot(A, A.T), np.identity(3), decimal=15)
    aae(np.dot(A.T, A), np.identity(3), decimal=15)

    aae(np.dot(B, B.T), np.dot(B.T, B), decimal=15)
    aae(np.dot(B, B.T), np.identity(3), decimal=15)
    aae(np.dot(B.T, B), np.identity(3), decimal=15)

    aae(np.dot(C, C.T), np.dot(C.T, C), decimal=15)
    aae(np.dot(C, C.T), np.identity(3), decimal=15)
    aae(np.dot(C.T, C), np.identity(3), decimal=15)


def test_R1_R2_R3_transposition():
    'R(-alpha) must be equal to transposed R(alpha)'
    A1 = coord.R1(-67)
    A2 = coord.R1(67).T
    aae(A1, A2, decimal=15)

    A1 = coord.R2(-17)
    A2 = coord.R2(17).T
    aae(A1, A2, decimal=15)

    A1 = coord.R3(-39)
    A2 = coord.R3(39).T
    aae(A1, A2, decimal=15)

    A1 = coord.R1(13)
    A2 = coord.R1(-13).T
    aae(A1, A2, decimal=15)

    A1 = coord.R2(8)
    A2 = coord.R2(-8).T
    aae(A1, A2, decimal=15)

    A1 = coord.R3(40)
    A2 = coord.R3(-40).T
    aae(A1, A2, decimal=15)
