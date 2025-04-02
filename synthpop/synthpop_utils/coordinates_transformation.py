"""
Functions for coordinates transformations following Bovy 2011.
"""
__all__ = ["get_trans_matrix", "getA", "lb_to_ad", "ad_to_lb", "dlb_to_xyz",
           "xyz_to_rphiz", "dlb_to_rphiz", "uvw_to_vrmulb", "uvw_to_vrmuad", "CoordTrans"]
__credits__ = ["J. Klüter", "S. Johnson", "M.J. Huston", "A. Aronica", "M. Penny"]
__date__ = "2022-07-08"

from typing import Tuple, Any
import numpy as np
from numpy import ndarray

try:
    from .. import constants as const
except (ImportError):
    import constants as const
from .utils_functions import rotation_matrix
from .sun_info import default_sun, SunInfo



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


trans_matrix = get_trans_matrix()


class CoordTrans:
    def __init__(self, sun: SunInfo = None,
            amp_warp: float = 0, amp_warp_pos: float = None, amp_warp_neg: float = None,
            r_warp: float = 0, alpha_warp: float = 0,
            phi_warp_deg: float = 0, phi_warp_rad: float = 0):

        self.sun = sun if sun is not None else default_sun
        self.amp_warp_pos = amp_warp_pos if amp_warp_pos is not None else amp_warp
        self.amp_warp_neg = amp_warp_neg if amp_warp_neg is not None else amp_warp
        self.r_warp = r_warp
        self.alpha_warp = alpha_warp
        self.phi_warp_rad = phi_warp_rad if phi_warp_rad else phi_warp_deg * np.pi / 180

    def warp_correction(self, r_kpc, phi_rad):
        """Correction for the WARP following Chen X. et al. 2019 """
        # only if r >= rWarp, else hw=0
        height_warp_pos = self.amp_warp_pos * np.maximum(r_kpc - self.r_warp, 0) ** self.alpha_warp
        height_warp_neg = self.amp_warp_neg * np.maximum(r_kpc - self.r_warp, 0) ** self.alpha_warp

        sin_warp = np.sin(phi_rad - self.phi_warp_rad)
        # split sin_warp into two components
        pos_dir = np.maximum(0, sin_warp)
        neg_dir = np.minimum(0, sin_warp)

        return height_warp_pos * pos_dir + height_warp_neg * neg_dir

    def dlb_to_rphiz(self, d_kpc: np.ndarray, l_deg: np.ndarray, b_deg: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        translates d, l, b  into  r, phi, z

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
        x, y, z = self.dlb_to_xyz(d_kpc, l_deg, b_deg)
        # convert into cylindrical coordinates including warp of the galaxy
        r_kpc, phi_rad, z_kpc = self.xyz_to_rphiz(x, y, z)
        return r_kpc, phi_rad, z_kpc

    def dlb_to_xyz(self, d_kpc: np.ndarray, l_deg: np.ndarray, b_deg: np.ndarray) \
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
        hc[0] -= self.sun.gal_dist

        # correction for tilt between galactic plane and direction to galactic center
        H = rotation_matrix(st=-self.sun.sin_theta, ct=self.sun.cos_theta, axis='y')
        x_kpc, y_kpc, z_kpc = np.dot(H, hc)
        return x_kpc, y_kpc, z_kpc

    def rphiz_to_xyz(self, r_kpc: np.ndarray, phi_rad: np.ndarray, z_kpc: np.ndarray) \
            -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
            Translate from cylindrical cartesian coordinates to
            galactocentric cartesian coordinates
            the zero point of z follows the warp of the galaxy.

            Parameters
            ----------
            r_kpc :  float, ndarray [kpc]
            phi_rad : float, ndarray [rad]
                polar angle follows the rotation of the galaxy
                zero point is at the sun.
            z_kpc : float, ndarray [kpc]
                hight above/below the galactic plane

            Returns
            -------
            x_kpc,y_kpc,z_kpc : float, ndarray [kpc]
                galactocentric cartesian coordinates
        """
        # remove correction of warp
        z_kpc += self.warp_correction(r_kpc, phi_rad)

        # estimate x and y
        x_kpc = -r_kpc * np.cos(phi_rad)
        y_kpc = r_kpc * np.sin(phi_rad)

        return x_kpc, y_kpc, z_kpc

    def xyz_to_rphiz(self, x_kpc: np.ndarray, y_kpc: np.ndarray, z_kpc: np.ndarray) \
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
        z_kpc : float, ndarray [kpc]
            hight above/below the galactic plane
        """

        # estimate r and theta
        r_kpc = np.sqrt(x_kpc ** 2 + y_kpc ** 2)
        phi_rad = np.arctan2(y_kpc, -x_kpc)  # zero point is at the sun (y=0,x~=-8)

        # adjust for warp
        z_kpc -= self.warp_correction(r_kpc, phi_rad)  # subtract position of the warp

        return r_kpc, phi_rad, z_kpc

    @staticmethod
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

    @staticmethod
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

    def uvw_to_vrmulb(self,
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
        uvw = [u_kmps - self.sun.u, v_kmps - self.sun.v, w_kmps - self.sun.w]

        # correction for tilt between galactic plane and direction to galactic center
        H = [[self.sun.cos_theta, 0, -self.sun.sin_theta],
             [0, 1, 0],
             [self.sun.sin_theta, 0, self.sun.cos_theta]]
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
            self, l_deg: np.ndarray, b_deg: np.ndarray, dist_kpc: np.ndarray,
            u_kmps: np.ndarray, v_kmps: np.ndarray, w_kmps: np.ndarray
            ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
        a_deg, d_deg = self.lb_to_ad(l_deg, b_deg)
        a_rad = a_deg * np.pi / 180
        d_rad = d_deg * np.pi / 180

        # subtract motion of the galaxy
        uvw = [u_kmps - self.sun.u, v_kmps - self.sun.v, w_kmps - self.sun.w]

        # correction for tilt between galactic plane and direction to galactic center
        H = rotation_matrix(st=-self.sun.sin_theta, ct=self.sun.cos_theta, axis='y')
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


# create wrappers for the default instance.
_coord_trans = CoordTrans()


def warp_correction(r_kpc, phi_rad):
    """Correction for the WARP following Chen X. et al. 2019 """
    return _coord_trans.warp_correction(r_kpc, phi_rad)


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
    return _coord_trans.dlb_to_rphiz(d_kpc, l_deg, b_deg)


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
    return _coord_trans.dlb_to_xyz(d_kpc, l_deg, b_deg)


def rphiz_to_xyz(r_kpc: np.ndarray, phi_rad: np.ndarray, z_kpc: np.ndarray)\
        -> ndarray:
    """
        Translate from cylindrical cartesian coordinates to
        galactocentric cartesian coordinates
        the zero point of z follows the warp of the galaxy.

        Parameters
        ----------
        r_kpc :  float, ndarray [kpc]
        phi_rad : float, ndarray [rad]
            polar angle follows the rotation of the galaxy
            zero point is at the sun.
        z_kpc : float, ndarray [kpc]
            hight above/below the galactic plane

        Returns
        -------
        x_kpc,y_kpc,z_kpc : float, ndarray [kpc]
            galactocentric cartesian coordinates
    """
    return _coord_trans.rphiz_to_xyz(r_kpc, phi_rad, z_kpc)


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
    z_kpc : float, ndarray [kpc]
        hight above/below the galactic plane
    """
    return _coord_trans.xyz_to_rphiz(x_kpc, y_kpc, z_kpc)


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
    return _coord_trans.lb_to_ad(l_deg, b_deg)


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
    return _coord_trans.ad_to_lb(a_deg, d_deg)


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
    return _coord_trans.uvw_to_vrmulb(l_deg, b_deg, dist_kpc, u_kmps, v_kmps, w_kmps)


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
    return _coord_trans.uvw_to_vrmuad(l_deg, b_deg, dist_kpc, u_kmps, v_kmps, w_kmps)
