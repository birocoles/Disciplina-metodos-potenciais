import numpy as np

def rotation_NEU(latitude, longitude):
    '''
    Compute the elements of a rotation matrix whose columns are
    the unit vectors v (North), w (East) and u (Up). If latitude is geodetic,
    then the computed u is normal to the referrence elipsoid. If latitude is
    spherical, then u is normal to the sphere. The unit vectors v and w point
    to directions of increasing latitude (spherical or geodetic) and longitude,
    respectively. At a given point with coordinates (lat, lon), the rotation
    matrix can be written as follows:

        | -sin(lat)*cos(lon)   -sin(lon)    cos(lat)*cos(lon) |
    R = | -sin(lat)*sin(lon)    cos(lon)    cos(lat)*sin(lon) |
        |       cos(lat)           0             sin(lat)     |

    input

    latitude: numpy array 1D - vector containing the latitude (in degrees)
               of the computation points.
    longitude: numpy array 1D - vector containing the lonitude (in degrees)
               of the computation points.

    output

    R: numpy array 2D - matrix whose columns contain the elements
       11, 12, 13, 21, 22, 23, 31 and 33 of the rotation matrix evaluated
       at the computation points. Notice that the element 32 is null and,
       consequently, it is not computed.

    '''
    latitude = np.asarray(latitude)
    longitude = np.asarray(longitude)

    assert latitude.size == longitude.size, 'latitude and longitude must have \
the same numer of elements'

    # convert degrees to radian
    lat = np.deg2rad(latitude)
    lon = np.deg2rad(longitude)

    coslat = np.cos(lat)
    sinlat = np.sin(lat)
    coslon = np.cos(lon)
    sinlon = np.sin(lon)

    R11 = -sinlat*coslon
    R21 = -sinlat*sinlon
    R31 = coslat

    R12 = -sinlon
    R22 = coslon

    R13 = coslat*coslon
    R23 = coslat*sinlon
    R33 = sinlat

    R = np.vstack([R11, R12, R13, R21, R22, R23, R31, R33]).T

    return R

def rotation_NED(latitude, longitude):
    '''
    Compute the elements of a rotation matrix whose columns are
    the unit vectors v (North), w (East) and -u (Down). If latitude is geodetic,
    then u is normal to the referrence elipsoid. If latitude is
    spherical, then u is normal to the sphere. The unit vectors v and w point
    to directions of increasing latitude (spherical or geodetic) and longitude,
    respectively. At a given point with coordinates (lat, lon), the rotation
    matrix can be written as follows:

        | -sin(lat)*cos(lon)   -sin(lon)   -cos(lat)*cos(lon) |
    R = | -sin(lat)*sin(lon)    cos(lon)   -cos(lat)*sin(lon) |
        |       cos(lat)           0            -sin(lat)     |

    input

    latitude: numpy array 1D - vector containing the latitude (in degrees)
               of the computation points.
    longitude: numpy array 1D - vector containing the lonitude (in degrees)
               of the computation points.

    output

    R: numpy array 2D - matrix whose columns contain the elements
       11, 12, 13, 21, 22, 23, 31 and 33 of the rotation matrix evaluated
       at the computation points. Notice that the element 32 is null and,
       consequently, it is not computed.

    '''
    latitude = np.asarray(latitude)
    longitude = np.asarray(longitude)

    assert latitude.size == longitude.size, 'latitude and longitude must have \
the same numer of elements'

    # convert degrees to radian
    lat = np.deg2rad(latitude)
    lon = np.deg2rad(longitude)

    coslat = np.cos(lat)
    sinlat = np.sin(lat)
    coslon = np.cos(lon)
    sinlon = np.sin(lon)

    R11 = -sinlat*coslon
    R21 = -sinlat*sinlon
    R31 = coslat

    R12 = -sinlon
    R22 = coslon

    R13 = -coslat*coslon
    R23 = -coslat*sinlon
    R33 = -sinlat

    R = np.vstack([R11, R12, R13, R21, R22, R23, R31, R33]).T

    return R


def unit_vector_normal(latitude, longitude):
    '''
    Compute the elements of a unit vector u pointing to the direction of
    the outward normal. If latitude is geodetic, then the computed u is
    normal to the referrence elipsoid. If latitude is spherical, then u is
    normal to the sphere. The vector u is referred to the Geocentric
    Cartesian System. At a given point with (spherical or geodetic) coordinates
    (lat, lon), the components of u along the axes X, Y and Z of a Geocentric
    Cartesian System can be written as follows:

        | cos(lat)*cos(lon) |
    u = | cos(lat)*sin(lon) | .
        | sin(lat)          |

    input

    latitude: numpy array 1D - vector containing the latitude (in degrees)
               of the computation points.
    longitude: numpy array 1D - vector containing the lonitude (in degrees)
               of the computation points.

    output

    uX, uY, uZ: numpy arrays 1D - components of the vector u.

    '''
    latitude = np.asarray(latitude)
    longitude = np.asarray(longitude)

    assert latitude.size == longitude.size, 'latitude and longitude must have \
the same numer of elements'

    # convert degrees to radian
    lat = np.deg2rad(latitude)
    lon = np.deg2rad(longitude)

    coslat = np.cos(lat)
    sinlat = np.sin(lat)
    coslon = np.cos(lon)
    sinlon = np.sin(lon)

    uX = coslat*coslon
    uY = coslat*sinlon
    uZ = sinlat

    return uX, uY, uZ


def unit_vector_latitude(latitude, longitude):
    '''
    Compute the elements of a unit vector v pointing to the direction of
    increasing latitude. If latitude is geodetic, then the computed v is
    tangent to the referrence elipsoid. If latitude is spherical, then v is
    tangent to the sphere. This vector is referred to the Geocentric
    Cartesian System. At a given point with (spherical or geodetic) coordinates
    (lat, lon), the components of v along the axes X, Y and Z of a Geocentric
    Cartesian System can be written as follows:

        | -sin(lat)*cos(lon) |
    v = | -sin(lat)*sin(lon) | .
        |  cos(lat)          |

    input

    latitude: numpy array 1D - vector containing the latitude (in degrees)
               of the computation points.
    longitude: numpy array 1D - vector containing the lonitude (in degrees)
               of the computation points.

    output

    vX, vY, vZ: numpy arrays 1D - components of the vector v.

    '''
    latitude = np.asarray(latitude)
    longitude = np.asarray(longitude)

    assert latitude.size == longitude.size, 'latitude and longitude must have \
the same numer of elements'

    # convert degrees to radian
    lat = np.deg2rad(latitude)
    lon = np.deg2rad(longitude)

    coslat = np.cos(lat)
    sinlat = np.sin(lat)
    coslon = np.cos(lon)
    sinlon = np.sin(lon)

    vX = -sinlat*coslon
    vY = -sinlat*sinlon
    vZ = coslat

    return vX, vY, vZ


def unit_vector_longitude(longitude):
    '''
    Compute the elements of a unit vector w pointing to the direction of
    increasing longitude. This vector is referred to the Geocentric
    Cartesian System. At a given point with (spherical or geodetic) coordinates
    (lat, lon), the components of w along the axes X, Y and Z of a Geocentric
    Cartesian System can be written as follows:

        | -sin(lon) |
    w = |  cos(lon) | .
        |  0        |

    input

    longitude: numpy array 1D - vector containing the longitude (in degrees)
               of the computation points.

    output

    wX, wY: numpy arrays 1D - non-null components of the vector w.

    '''
    longitude = np.asarray(longitude)

    # convert degrees to radian
    lon = np.deg2rad(longitude)

    coslon = np.cos(lon)
    sinlon = np.sin(lon)

    wX = -sinlon
    wY = coslon
    #wZ = 0.

    return wX, wY


def unit_vector_TCCS(inclination, declination):
    '''
    Compute the components x, y and z of unit vectors v referred to a
    Topocentric Cartesian Coordinate System (TCCS). Each vector is defined by
    a pair (inclination, declination) as follows:
        | cos(inclination)*cos(declination) |
    v = | cos(inclination)*sin(declination) | .
        | sin(inclination)                  |

    input

    inclination: numpy array 1D or float - inclination values in degrees.
    declination: numpy array 1D or float - declination values in degrees.

    output

    vx, vy, vz: numpy arrays 1D or floats - components of the unit vector(s).

    '''
    inclination = np.asarray(inclination)
    declination = np.asarray(declination)

    assert inclination.size == declination.size, 'inclination and declination \
must have the same numer of elements'

    # convert degrees to radian
    inc = np.deg2rad(inclination)
    dec = np.deg2rad(declination)

    cosinc = np.cos(inc)
    sininc = np.sin(inc)
    cosdec = np.cos(dec)
    sindec = np.sin(dec)

    vx = cosinc*cosdec
    vy = cosinc*sindec
    vz = sininc

    return vx, vy, vz


# This is a copy of the Fatiando a Terra routine
# fatiando.gridder.regular
# https://www.fatiando.org/api/gridder.html#fatiando.gridder.regular
def regular_grid(area, shape, z=None):
    """
    Create a regular grid.

    The x directions is North-South and y East-West. Imagine the grid as a
    matrix with x varying in the lines and y in columns.

    Returned arrays will be flattened to 1D with ``numpy.ravel``.

    .. warning::

        As of version 0.4, the ``shape`` argument was corrected to be
        ``shape = (nx, ny)`` instead of ``shape = (ny, nx)``.


    Parameters:

    * area
        ``(x1, x2, y1, y2)``: Borders of the grid
    * shape
        Shape of the regular grid, ie ``(nx, ny)``.
    * z
        Optional. z coordinate of the grid points. If given, will return an
        array with the value *z*.

    Returns:

    * ``[x, y]``
        Numpy arrays with the x and y coordinates of the grid points
    * ``[x, y, z]``
        If *z* given. Numpy arrays with the x, y, and z coordinates of the grid
        points

    Examples::

        >>> x, y = regular((0, 10, 0, 5), (5, 3))
        >>> x
        array([  0. ,   0. ,   0. ,   2.5,   2.5,   2.5,   5. ,   5. ,   5. ,
                 7.5,   7.5,   7.5,  10. ,  10. ,  10. ])
        >>> x.reshape((5, 3))
        array([[  0. ,   0. ,   0. ],
               [  2.5,   2.5,   2.5],
               [  5. ,   5. ,   5. ],
               [  7.5,   7.5,   7.5],
               [ 10. ,  10. ,  10. ]])
        >>> y.reshape((5, 3))
        array([[ 0. ,  2.5,  5. ],
               [ 0. ,  2.5,  5. ],
               [ 0. ,  2.5,  5. ],
               [ 0. ,  2.5,  5. ],
               [ 0. ,  2.5,  5. ]])
        >>> x, y = regular((0, 0, 0, 5), (1, 3))
        >>> x.reshape((1, 3))
        array([[ 0.,  0.,  0.]])
        >>> y.reshape((1, 3))
        array([[ 0. ,  2.5,  5. ]])
        >>> x, y, z = regular((0, 10, 0, 5), (5, 3), z=-10)
        >>> z.reshape((5, 3))
        array([[-10., -10., -10.],
               [-10., -10., -10.],
               [-10., -10., -10.],
               [-10., -10., -10.],
               [-10., -10., -10.]])


    """
    nx, ny = shape
    x1, x2, y1, y2 = area
    assert x1 <= x2, \
        "Invalid area dimensions {}, {}. x1 must be < x2.".format(x1, x2)
    assert y1 <= y2, \
        "Invalid area dimensions {}, {}. y1 must be < y2.".format(y1, y2)
    xs = np.linspace(x1, x2, nx)
    ys = np.linspace(y1, y2, ny)
    # Must pass ys, xs in this order because meshgrid uses the first argument
    # for the columns
    arrays = np.meshgrid(ys, xs)[::-1]
    if z is not None:
        arrays.append(z*np.ones(nx*ny, dtype=np.float))
    return [i.ravel() for i in arrays]


def R1(angle):
    '''
    Orthogonal matrix performing a rotation around
    the x-axis of a Cartesian coordinate system.

    Parameters:
    * angle : float
        Rotation angle (in degrees).

    Returns:
    * R : 2D numpy array
        Rotation matrix.
    '''

    assert isinstance(1.*angle, float), 'angle must be a float'

    ang = np.deg2rad(angle)

    cos_angle = np.cos(ang)
    sin_angle = np.sin(ang)

    R = np.array([[1, 0, 0],
                  [0, cos_angle, sin_angle],
                  [0, -sin_angle, cos_angle]])

    return R


def R2(angle):
    '''
    Orthogonal matrix performing a rotation around
    the y-axis of a Cartesian coordinate system.

    Parameters:
    * angle : float
        Rotation angle (in degrees).

    Returns:
    * R : 2D numpy array
        Rotation matrix.
    '''

    assert isinstance(1.*angle, float), 'angle must be a float'

    ang = np.deg2rad(angle)

    cos_angle = np.cos(ang)
    sin_angle = np.sin(ang)

    R = np.array([[cos_angle, 0, -sin_angle],
                  [0, 1, 0],
                  [sin_angle, 0, cos_angle]])

    return R


def R3(angle):
    '''
    Orthogonal matrix performing a rotation around
    the z-axis of a Cartesian coordinate system.

    Parameters:
    * angle : float
        Rotation angle (in degrees).

    Returns:
    * R : 2D numpy array
        Rotation matrix.
    '''

    assert isinstance(1.*angle, float), 'angle must be a float'

    ang = np.deg2rad(angle)

    cos_angle = np.cos(ang)
    sin_angle = np.sin(ang)

    R = np.array([[cos_angle, sin_angle, 0],
                  [-sin_angle, cos_angle, 0],
                  [0, 0, 1]])

    return R
