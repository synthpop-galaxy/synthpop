"""
functions for coordinates transformations
following Bovy 2011.
"""
__all__ = ["rotation_matrix", "get_trans_matrix", "getA", "lb_to_ad", "ad_to_lb", "dlb_to_xyz",
           "xyz_to_rphiz", "dlb_to_rphiz", "uvw_to_vrmulb", "uvw_to_vrmuad"]
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__license__ = "GPLv3"
__date__ = "2022-07-08"
__version__ = '1.0.0'

from typing import Tuple
import numpy as np

try:
    from ... import constants as const
except (ImportError, ValueError):
    import constants as const


def rotation_matrix(
        theta_rad: float or np.ndarray or None = None,
        st: float or np.ndarray or None = None, ct: float or np.ndarray or None = None,
        axis: str = ''
        ) -> np.ndarray:
    """
    creates the rotation matrix along a given axis.
    These are:

         | 1    0    0  |
    RX = | 0   cos -sin |
         | 0   sin  cos |

         | cos  0  -sin |
    RY = |  0   1    0  |
         | sin  0   cos |

         | cos -sin   0  |
    RZ = | sin  cos   0  |
         |  0    0    1  |

    Parameters
    ----------
    theta_rad: float or np.ndarray or None [rad]
        rotation angle
    st, ct: float or np.ndarray or None [rad]
        sin(theta) and cos(theta)
    axis: str
        string to specify the axis. either 'x', 'y', or 'z'

    Returns
    -------
    rot_matrix : np.ndarray
        rotation matrix
    """

    if axis not in ['x', 'y', 'z']:
        raise ValueError("axis has to be 'x', 'y', or 'z'!")
    if theta_rad is None and (st is None or ct is None):
        raise ValueError("either theta or ct and st has to be specified")

    if st is None or ct is None:
        st = np.sin(theta_rad)
        ct = np.cos(theta_rad)

    if isinstance(st, np.ndarray):
        one = np.ones(st.shape)
        zero = np.zeros(st.shape)
    else:
        one = 1
        zero = 0

    if axis == 'x':
        return np.array([[one, zero, zero], [zero, ct, -st], [zero, st, ct]])
    if axis == 'y':
        return np.array([[ct, zero, -st], [zero, one, zero], [st, zero, ct]])
    if axis == 'z':
        return np.array([[ct, -st, zero], [st, ct, zero], [zero, zero, one]])


def get_trans_matrix() -> np.ndarray:
    """
    creates the transformation matrix between galactic and equatorial coordinates

    Returns
    -------
    trans_matrix :
        transformation matrix between galactic and equatorial coordinates
    """
    rot_1 = rotation_matrix(const.L_NCP_DEG * np.pi / 180, axis='z')
    flip = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    rot_2 = rotation_matrix(const.D_NGP_DEG * np.pi / 180 + np.pi / 2, axis='y')
    rot_3 = rotation_matrix(-const.A_NGP_DEG * np.pi / 180, axis='z')
    return np.matmul(np.matmul(np.matmul(rot_1, rot_2), flip), rot_3)


trans_matrix = get_trans_matrix()


def getA(longitude_rad: np.ndarray or float, latitude_rad: np.ndarray or float) \
        -> np.ndarray:
    """
    determines the Position Matrix A from Bovy 2011
    Note that R^T is equivalent to A

    Parameters
    ----------
    longitude_rad: float, ndarray
        galactic or ecliptic longitude
    latitude_rad: float, ndarray
        galactic or ecliptic latitude

    Returns
    -------
    A : ndarray
        transformation matrix
    """
    clong, slong = np.cos(longitude_rad), np.sin(longitude_rad)
    clat, slat = np.cos(latitude_rad), np.sin(latitude_rad)

    r_z_long = rotation_matrix(ct=clong, st=slong, axis='z')
    r_y_lat = rotation_matrix(ct=clat, st=slat, axis='y')
    A = np.matmul(r_z_long, r_y_lat, axes=[(0, 1), (0, 1), (0, 1)])
    return A


def lb_to_ad(l_deg: np.ndarray, b_deg: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Conversion from Galactic coordinates to Equatorial coordinates

    Parameters
    ----------
    l_deg : float, ndarray [degree]
        Galactic longitude
    b_deg : float, ndarray [degree]
        Galactic latitude

    Returns
    -------
    a_deg : float, ndarray [degree]
        right ascension
    d_deg : float, ndarray [degree]
        declination
    """
    # translate variables to radians
    l_rad = l_deg * np.pi / 180.
    b_rad = b_deg * np.pi / 180.

    # calculate trig functions using coordinates in radians
    cl, sl, cb, sb = np.cos(l_rad), np.sin(l_rad), np.cos(b_rad), np.sin(b_rad)

    # rotate coordinate system
    init_array = [cb * cl, cb * sl, sb]
    new_mat = np.dot(trans_matrix.T, init_array)

    # estimate ra and dec
    d_rad = np.arcsin(new_mat[2])
    a_rad = np.arctan2(new_mat[1], new_mat[0])

    # use angle fom 0 to 2pi
    a_rad += 2 * np.pi * (a_rad < 0)  # only if a_rad < 0
    # this way it works with both ndarray and floats

    # convert  to degree
    a_deg = a_rad * 180 / np.pi
    d_deg = d_rad * 180 / np.pi
    return a_deg, d_deg


def ad_to_lb(a_deg: np.ndarray, d_deg: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Conversion from Equatorial coordinates to Galactic coordinates

    Parameters
    ----------
    a_deg : float, ndarray [degree]
        right ascension
    d_deg : float, ndarray [degree]
        declination

    Returns
    -------
    l_deg : float, ndarray [degree]
        Galactic longitude
    b_deg : float, ndarray [degree]
        Galactic latitude
    """
    # translate variables to radians
    a_rad = a_deg * np.pi / 180
    d_rad = d_deg * np.pi / 180

    # estimate cartesian coordinates
    ca, sa, cd, sd = np.cos(a_rad), np.sin(a_rad), np.cos(d_rad), np.sin(d_rad)
    init_mat = [cd * ca, cd * sa, sd]

    # rotate coordinate system
    new_mat = np.dot(trans_matrix, init_mat)

    # estimate l and b
    b_rad = np.arcsin(new_mat[2])
    l_rad = np.arctan2(new_mat[1], new_mat[0])

    # convert  to degree
    l_deg = l_rad * 180 / np.pi
    b_deg = b_rad * 180 / np.pi
    return l_deg, b_deg


def dlb_to_xyz(d_kpc: np.ndarray, l_deg: np.ndarray, b_deg: np.ndarray) \
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert from galactic coordinates an distance
    to galactocentric cartesian coordinates.

    Parameters
    ----------
    l_deg : float, ndarray, [degrees]
        galactic longitude
    b_deg : float, ndarray, [degrees]
        galactic latitude
    d_kpc : float ndarray, [kpc]
        distance
    Returns
    -------
    x,y,z : float nd_array [kpc]
        galactocentric cartesian coordinates
    """

    # convert to radian
    l_rad = l_deg * np.pi / 180.
    b_rad = b_deg * np.pi / 180.

    # convert to heliocentric cartesian
    cl = np.cos(l_rad)
    cb = np.cos(b_rad)
    sl = np.sin(l_rad)
    sb = np.sin(b_rad)
    hc = np.array([d_kpc * cb * cl, d_kpc * cb * sl, d_kpc * sb])

    # shift zero point to galactic center
    hc[0] -= const.GAL_DIST

    # correction for tilt between galactic plane and direction to galactic center
    H = rotation_matrix(st=-const.ST_SUN, ct=const.CT_SUN, axis='y')
    x_kpc, y_kpc, z_kpc = np.dot(H, hc)
    return x_kpc, y_kpc, z_kpc


def xyz_to_rphiz(x_kpc: np.ndarray, y_kpc: np.ndarray, z_kpc: np.ndarray) \
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Translate from galactocentric cartesian coordinates into
    cylindrical cartesian coordinates,
    the zero point of z follows the warp of the galaxy.

    Parameters
    ----------
    x_kpc,y_kpc,z_kpc : float, ndarray [kpc]
        galactocentric cartesian coordinates

    Returns
    -------
    r_kpc :  float, ndarray [kpc]
    phi_rad : float, ndarray [rad]
        polar angle follows the rotation of the galaxy
        zero point is at the sun.
    z : float, ndarray [kpc]
        hight above/below the galactic plane
    """

    # estimate r and theta
    r_kpc = np.sqrt(x_kpc ** 2 + y_kpc ** 2)
    phi_rad = np.arctan2(y_kpc, -x_kpc)  # zero point is at the sun (y=0,x~=-8)

    # adjust for warp
    z_kpc -= warp_correction(r_kpc, phi_rad)  # subtract position of the warp

    return r_kpc, phi_rad, z_kpc


def warp_correction(r_kpc, phi_rad):
    """Correction for the WARP following Chen X. et al. 2019 """
    # only if r >= rWarp, else hw=0
    hw = const.AMP_WARP * np.maximum(r_kpc - const.R_WARP, 0) ** const.ALPHA_WARP
    return hw * np.sin(phi_rad - const.PHI_WARP)


def dlb_to_rphiz(d_kpc: np.ndarray, l_deg: np.ndarray, b_deg: np.ndarray) \
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    translates d, l, b  into  r, theta, z

    phi increases along the galactic rotation with a zero point at the position of the sun


    Parameters
    ----------
    d_kpc : float, ndarray [kpc]
    l_deg :  float, ndarray [degree]
    b_deg :  float, ndarray [degree]

    Returns
    -------
    r :  float, ndarray [kpc]
    phi_rad : float, ndarray [rad]
        polar angle follows the rotation of the galaxy
        zero point is at the sun.
    z : float, ndarray [kpc]
        hight above/below the galactic plane
    """
    # convert int galactocentric cartesian coordinates
    x, y, z = dlb_to_xyz(d_kpc, l_deg, b_deg)
    # convert into cylindrical coordinates including warp of the galaxy
    r_kpc, phi_rad, z_kpc = xyz_to_rphiz(x, y, z)
    return r_kpc, phi_rad, z_kpc


def uvw_to_vrmulb(
        l_deg: np.ndarray, b_deg: np.ndarray, dist_kpc: np.ndarray,
        u_kmps: np.ndarray, v_kmps: np.ndarray, w_kmps: np.ndarray
        ) \
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """

    Conversion from u,v,w to v_r, mu_l, and mu_b
    Using (vr, vl, vb) = R * (u,v,w)  from Bovy 2011

    Parameters
    ----------
    l_deg : float, ndarray [degrees]
        galactic longitude
    b_deg :  float, ndarray [degrees]
        galactic latitude
    dist_kpc : float, ndarray [kpc]
        distance
    u_kmps : float, ndarray [km/s]
        rectangular velocity with_out correction for the motion of the sun
    v_kmps : float, ndarray [km/s]
        rectangular velocity with_out correction for the motion of the sun
    w_kmps : float, ndarray [km/s]
         rectangular velocity with_out correction for the motion of the sun

    Returns
    -------
    vr_kmps : float, ndarray[km/s]
        radial velocity
    mu_l_maspyr : float, ndarray [mas/yr]
        proper motion in galactic longitude  cos(b) is applied
    mu_b_maspyr : float, ndarray [mas/yr]
        proper motion in galactic latitude
    """
    # convert to radian
    l_rad = l_deg * np.pi / 180.
    b_rad = b_deg * np.pi / 180.

    # subtract motion of the galaxy
    uvw = [u_kmps - const.U_SUN, v_kmps - const.V_SUN, w_kmps - const.W_SUN]

    # correction for tilt between galactic plane and direction to galactic center
    H = [[const.CT_SUN, 0, -const.ST_SUN], [0, 1, 0], [const.ST_SUN, 0, const.CT_SUN]]
    uvw = np.dot(H, uvw)

    # get transformation matrix
    A = getA(l_rad, b_rad)
    iA = np.swapaxes(A, 0, 1)  # I.E A^T for all stars separately

    if len(iA.shape) == 3:
        vr_vt = np.sum(iA * uvw, axis=1)
    else:
        vr_vt = np.dot(iA, uvw)

    vr_kmps = vr_vt[0]
    # transfer tangential velocity to proper motion
    mu_l_maspyr, mu_b_maspyr = vr_vt[1:] / dist_kpc / const.MUxD_to_VT

    return vr_kmps, mu_l_maspyr, mu_b_maspyr


def uvw_to_vrmuad(
        l_deg: np.ndarray, b_deg: np.ndarray, dist_kpc: np.ndarray,
        u_kmps: np.ndarray, v_kmps: np.ndarray, w_kmps: np.ndarray
        ) \
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """

        Conversion from u,v,w
        to v_r, mu_l, and mu_b

        Parameters
        ----------
        l_deg : float, ndarray [degrees]
            galactic longitude
        b_deg :  float, ndarray [degrees]
            galactic latitude
        dist_kpc : float, ndarray [kpc]
            distance
        u_kmps : float, ndarray [km/s]
            rectangular velocity with_out correction for the motion of the sun
        v_kmps : float, ndarray [km/s]
            rectangular velocity with_out correction for the motion of the sun
        w_kmps : float, ndarray [km/s]
             rectangular velocity with_out correction for the motion of the sun

        Returns
        -------
        vr_kmps : float, ndarray[km/s]
            radial velocity
        mu_a_maspyr : float, ndarray [mas/yr]
            proper motion in galactic longitude cos(delta) is applied
        mu_d_maspyr : float, ndarray [mas/yr]
            proper motion in galactic latitude
        """

    # Galactic to Equatorial coordinate transform.
    a_deg, d_deg = lb_to_ad(l_deg, b_deg)
    a_rad = a_deg * np.pi / 180
    d_rad = d_deg * np.pi / 180

    # subtract motion of the galaxy
    uvw = [u_kmps - const.U_SUN, v_kmps - const.V_SUN, w_kmps - const.W_SUN]

    # correction for tilt between galactic plane and direction to galactic center
    H = rotation_matrix(st=-const.ST_SUN, ct=const.CT_SUN, axis='y')
    uvw = np.dot(H, uvw)

    # Make coordinate transform
    A = getA(a_rad, d_rad)
    iA = np.swapaxes(A, 0, 1)

    vr_vt = np.dot(trans_matrix.T, uvw)
    if len(iA.shape) == 3:
        vr_vt = np.sum(iA * vr_vt, axis=1)
    else:
        vr_vt = np.dot(iA, vr_vt)

    vr = vr_vt[0]
    mu_a_maspyr, mu_d_maspyr = vr_vt[1:] / dist_kpc / const.MUxD_to_VT

    return vr, mu_a_maspyr, mu_d_maspyr
