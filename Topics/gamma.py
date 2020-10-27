import numpy as np

def WGS84():
    '''
    This function returns the following parameters defining the 
    reference elipsoid WGS84:
    a = semimajor axis [m]
    f = flattening
    GM = geocentric gravitational constant of the Earth 
         (including the atmosphere) [m**3/s**2]
    omega = angular velocity [rad/s]
    
    output:
    a, f, GM, omega
    '''
    a = 6378137.0
    f = 1.0/298.257223563
    GM = 3986004.418*(10**8)
    omega = 7292115*(10**-11)
    
    return a, f, GM, omega

def somigliana(a, f, GM, omega, phi):
    '''
    This function calculates the normal gravity by using
    the Somigliana's formula.
    
    input:
    a: float containing the semimajor axis [m]
    f: float containing the flattening
    GM: float containing the geocentric gravitational constant 
        of the Earth (including the atmosphere) [m**3/s**2]
    omega: float containing the angular velocity [rad/s]
    phi: array containing the geodetic latitudes [degree]
    
    output:
    gamma: array containing the values of normal gravity
           on the surface of the elipsoid for each geodetic
           latitude [mGal]
    '''
    b = a*(1.0-f)
    a2 = a**2
    b2 = b**2
    E = np.sqrt(a2 - b2)
    elinha = E/b
    bE = b/E
    Eb = E/b
    atg = np.arctan(Eb)
    q0 = 0.5*((1+3*(bE**2))*atg - (3*bE))
    q0linha = 3.0*(1+(bE**2))*(1-(bE*atg)) - 1
    m = (omega**2)*(a2)*b/GM
    aux = elinha*q0linha/q0
    gammaa = (GM/(a*b))*(1-m-(m/6.0)*aux)
    gammab = (GM/a2)*(1+(m/3.0)*aux)
    aux = np.deg2rad(phi)
    s2 = np.sin(aux)**2
    c2 = np.cos(aux)**2
    # the 10**5 converts from m/s**2 to mGal
    gamma = (10**5)*((a*gammaa*c2) + (b*gammab*s2))/np.sqrt((a2*c2) + (b2*s2))
    return gamma
    
def closedform(a, f, GM, omega, phi, h):
    '''
    This function calculates the normal gravity by using
    a closed-form formula.
    
    input:
    a: float containing the semimajor axis [m]
    f: float containing the flattening
    GM: float containing the geocentric gravitational constant 
        of the Earth (including the atmosphere) [m**3/s**-2]
    omega: float containing the angular velocity [rad/s]
    phi: array containing the geodetic latitudes [degree]
    h: array containing the normal heights [m]
    
    output:
    gamma: array containing the values of normal gravity
           on the surface of the elipsoid for each geodetic
           latitude [mGal]
    '''
    b = a*(1.0-f)
    a2 = a**2
    b2 = b**2
    E = np.sqrt(a2 - b2)
    E2 = E**2
    bE = b/E
    Eb = E/b
    atanEb = np.arctan(Eb)
    phirad = np.deg2rad(phi)
    tanphi = np.tan(phirad)
    cosphi = np.cos(phirad)
    sinphi = np.sin(phirad)
    beta = np.arctan(b*tanphi/a)
    sinbeta = np.sin(beta)
    cosbeta = np.cos(beta)
    zl = b*sinbeta+h*sinphi
    rl = a*cosbeta+h*cosphi
    zl2 = zl**2
    rl2 = rl**2
    dll2 = rl2-zl2
    rll2 = rl2+zl2
    D = dll2/E2
    R = rll2/E2
    cosbetal = np.sqrt(0.5*(1+R) - np.sqrt(0.25*(1+R**2) - 0.5*D))
    cosbetal2 = cosbetal**2
    sinbetal2 = 1-cosbetal2
    bl = np.sqrt(rll2 - E2*cosbetal2)
    bl2 = bl**2
    blE = bl/E
    Ebl = E/bl
    atanEbl = np.arctan(Ebl)
    q0 = 0.5*((1+3*(bE**2))*atanEb - (3*bE))
    q0l = 3.0*(1+(blE**2))*(1-(blE*atanEbl)) - 1
    W = np.sqrt((bl2+E2*sinbetal2)/(bl2+E2))

    gamma = GM/(bl2+E2) - cosbetal2*bl*omega**2
    gamma += (((omega**2)*a2*E*q0l)/((bl2+E2)*q0))*(0.5*sinbetal2 - 1./6.)
    # the 10**5 converts from m/s**2 to mGal
    gamma = (10**5)*gamma/W
    
    return gamma
