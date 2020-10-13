# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 13:49:52 2020

@author: schaecl

Script to propagate Gaussian beams via the ABCD formalism.
"""

import numpy as np

class qbeam:
    """
    Class to describe a q parameter of a Gaussian beam."""
    z_re = 0
    zr_im = 0
    wlen = 0
    name = ""
    DBL_MAX = np.finfo(float).max

    def __init__(self, z_start: float, zr_start: float, wlen: float, name: str):
        self.z_re = z_start
        self.zr_im = zr_start
        self.wlen = wlen
        self.name = name

    def set_q(self, z_start: float, zr_start: float, wlen: float):
        """
        Sets the parameters.

        Parameters
        ----------
        z_start : float
            z position at which the beam starts.
        zr_start : float
            Rayleigh range.
        wlen : float
            Wavelength.

        Returns
        -------
        None.

        """
        self.z_re = z_start
        self.zr_im = zr_start
        self.wlen = wlen

    def get_w0(self):
        """
        Returns the beam waist.

        Returns
        -------
        flaot
            Beam waist.

        """
        return np.sqrt(self.zr_im * self.wlen / np.pi)

    def get_z(self):
        """
        Returns the position at which the waist is to be found.

        Returns
        -------
        float
            Position z0.

        """
        return self.z_re

    def get_zr(self):
        """
        Returns the Rayleigh range.

        Returns
        -------
        float
            Rayleigh range.

        """
        return self.zr_im

    def get_beamrad(self):
        """
        Returns the wavefront curvature at z.
        TODO Implement function of z.

        Returns
        -------
        float
            Wavefront curvature at z.

        """
        r = self.z_re
        if r != 0.:
            r += (self.zr_im**2 / r)
        else:
            return self.DBL_MAX
        return r

    def get_wz(self):
        """
        Return the beam radius at z.
        TODO Implement function of z.

        Returns
        -------
        wz : float
            Beam radius at z.

        """
        wz = self.zr_im
        wz += (self.z_re**2 / self.zr_im)
        wz *= (self.wlen / np.pi)
        wz = np.sqrt(wz)
        return wz

    def get_div(self):
        """
        Returns the far field divergence of the beam.

        Returns
        -------
        float
            Far field divergence (half-angle).

        """
        return np.atan(np.sqrt(self.wlen / (self.zr_im * np.pi)))

    def q_out(self):
        """
        Prints the q parameter and beam properties.

        Returns
        -------
        None.

        """
        if self.get_beamrad() == self.DBL_MAX:
            print("q.{:6.6s}: {:12s} {:12s}\n{:9s} {:12g} {:12g}\n"
                  "{:9s} {:12s} {:12s}\n{:9s} {:12s} {:12g}\n".format(
                      self.name, "z[e-3m]", "w0[e-6m]", "",
                      self.get_z(), self.get_w0() * 1e3, "",
                      "R[e-3m]", "wz[e-6m]", "", "inf", self.get_wz() * 1e3))
        else:
            print("q.{:6.6s}: {:12s} {:12s}\n{:9s} {:12g} {:12g}\n"
                  "{:9s} {:12s} {:12s}\n{:9s} {:12g} {:12g}\n".format(
                      self.name, "z[e-3m]", "w0[e-6m]", "",
                      self.get_z(), self.get_w0() * 1e3, "",
                      "R[e-3m]", "wz[e-6m]", "", self.get_beamrad(), self.get_wz() * 1e3))

def thin_lens(f):
    """
    Returns an ABCD matrix of a thin lens.

    Parameters
    ----------
    f : float
        Focal length.

    Returns
    -------
    M : array, float
        ABCD matrix.

    """
    M = np.array([[0., 0.], [0., 0.]])
    M[0][0] = 1.
    M[0][1] = 0.
    M[1][0] = -1. / f
    M[1][1] = 1.
    return M

def thick_lens(r1, r2, d, n):
    """
    Returns an ABCD matrix for a thick lens.

    Parameters
    ----------
    r1 : float
        First surface ROC.
    r2 : float
        Second surface ROC.
    d : float
        Thickness.
    n : float
        Refractive index.

    Returns
    -------
    M : array, float
        ABCD Matrix.

    """
    M = np.array([[0., 0.], [0., 0.]])
    ni = 1. / n
    dr1 = 0.
    r1i = 0.
    r2i = 0.
    if r1 != 0.:
        r1i = 1. / r1

    dr1 = d * r1i
    if r2 != 0.:
        r2i = 1. / r2

    M[0][0] = 1. + dr1 * (ni - 1.)
    M[0][1] = d * ni
    M[1][0] = (1. - n) * (r1i - r2i) + dr1 * r2i * (2. - ni - n)
    M[1][1] = 1. - d * r2i * (ni - 1.)

    return M

def thick_p_ck_lens(r2, d, n):
    """
    Returns an ABCD matrix for a plane-concave lens.

    Parameters
    ----------
    r2 : float
        Second surface ROI.
    d : float
        Thickness of the lens.
    n : float
        Refractive index.

    Returns
    -------
    M : array, float
        ABCD matrix.

    """
    M = np.array([[0., 0.], [0., 0.]])
    ni = 1. / n
    M[0][0] = 1.
    M[0][1] = d * ni
    M[1][0] = -(1. - n) / r2
    M[1][1] = 1. - d / r2 * (ni - 1.)
    return M

def propagation(d):
    """
    Returns an ABCD matrix to propagate a beam.

    Parameters
    ----------
    d : float
        Distance along the beam axis.

    Returns
    -------
    M : array, float
        ABCD matrix.

    """
    M = np.array([[0., 0.], [0., 0.]])
    M[0][0] = 1.
    M[0][1] = d
    M[1][0] = 0.
    M[1][1] = 1.
    return M

def refraction(n1, n2, r):
    """
    Returns and ABCD matrix for refraction.

    Parameters
    ----------
    n1 : float
        First index of refraction.
    n2 : float
        Second index of refraction.
    r : float
        ROC.

    Returns
    -------
    M : array, float
        ABCD matrix.

    """
    M = np.array([[0., 0.], [0., 0.]])
    M[0][0] = 1.
    M[0][1] = 0.
    if r != 0.:
        M[1][0] = (n1 - n2) / (n2 * r)
    else:
        M[1][0] = 0.
    M[1][1] = n1 / n2
    return M

def reflection(r):
    """
    Returns an ABCD matrix for reflection.

    Parameters
    ----------
    r : float
        ROC.

    Returns
    -------
    M : array, float
        ABCD matrix.

    """
    M = np.array([[0., 0.], [0., 0.]])
    M[0][0] = 1.
    M[0][1] = 0.
    M[1][0] = -2. / r
    M[1][1] = 1.
    return M

def grin(d, n0, nmin):
    """
    Returns an ABCD matrix for a GRIN lens.

    Parameters
    ----------
    d : float
        Length of the lens.
    n0 : float
        Nominal refractive index.
    nmin : float
        Minimum refractive index.

    Returns
    -------
    M : array, float
        ABCD matrix.

    """
    M = np.array([[0., 0.], [0., 0.]])
    t1 = np.sqrt(nmin / n0)
    t2 = t1 * d
    t3 = np.sin(t2)
    M[0][0] = np.cos(t2)
    M[0][1] = t3 / t1
    M[1][0] = -t1 * t3
    M[1][1] = M[0][0]
    return M

def qtrans(M, q):
    """
    Transports a beam from q to qp via an ABCD matrix.

    Parameters
    ----------
    M : float
        ABCD matrix.
    q : qbeam
        Initial q parameter.

    Returns
    -------
    qp : qbeam
        Propagated qbeam.

    """
    qp = qbeam(0., 0., q.wlen, q.name)

    t1 = M[1][0] * q.z_re + M[1][1]
    t2 = t1 * t1
    t3 = M[0][0] * q.z_re + M[0][1]
    t4 = M[1][0]**2
    t5 = q.zr_im**2

    denom = t2 + t4 * t5
    if denom == 0.:
        print("singularity in qtrans. healing.")
        denom = 1e-15

    qp.z_re = t1 * t3 + t5 * M[0][0] * M[1][0]
    qp.z_re /= denom
    qp.zr_im = q.zr_im * (M[0][0] * t1 - M[1][0] * t3)
    qp.zr_im /= denom
    return qp

def wz(w0, z, zr):
    """
    Returns the beam radius at position z.

    Parameters
    ----------
    w0 : float
        Beam waist.
    z : float
        Position along the propagation axis where the radius should be
        evaluated.
    zr : float
        Rayleigh range.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    t = z / zr
    t *= t
    return w0 * np.sqrt(1. + t)

def rayleighr(w0, lam):
    """
    Returns the Rayleigh range.

    Parameters
    ----------
    w0 : float
        Beam waist.
    lam : float
        Wavelength.

    Returns
    -------
    float
        Rayleigh range.

    """
    return np.pi * w0**2 / lam

def stab_parm(M):
    """
    Return stability parameter of an ABCD matrix.

    Parameters
    ----------
    M : array, float
        ABCD matrix.

    Returns
    -------
    float
        Stability parameter.

    """
    return (M[0][0] + M[1][1]) / 2.

def eigenval_res(M):
    """
    Returns the eigenvalue(s) of an ABCD matrix.

    Parameters
    ----------
    M : array, float
        ABCD matrix.

    Returns
    -------
    complex
        First eigenvalue.
    complex
        Second eigenvalue.

    """
    c1_re = c1_im = 0.
    c2_re = c2_im = 0.
    g = stab_parm(M)
    t = g**2 - 1.
    if t < 0.:
        t = np.sqrt(-t)
        c1_re = c2_re = g
        c1_im = t
        c2_im = -t
    elif t > 0.:
        t = np.sqrt(t)
        c1_im = c2_im = 0.
        c1_re = g + t
        c2_re = g - t
    else:
        c1_im = c2_im = 0.
        c1_re = c2_re = g
    return complex(c1_re, c1_im), complex(c2_re, c2_im)

def eigenmod_res(M):
    """
    Return (if any) the eigenmodes of a cavity given by an ABCD matrix.

    Parameters
    ----------
    M : array, float
        ABCD matrix.

    Returns
    -------
    complex
        Eigenmode (value) of the matrix.

    """
    c_re = c_im = 0
    t1 = M[1][1] - M[0][0]
    if t1 != 0.:
        t3 = t1**2
    else:
        c_re = 0.
        t3 = M[0][1] / M[1][0]
        if t3 >= 0.:
            print("\nfound no eigenmode solution")
            c_im = 0.
        else:
            c_im = np.sqrt(-t3)
        return complex(c_re, c_im)

    t2 = 2. * M[1][0]
    c_re = t1 / t2
    t2 *= t2
    t3 = (t3 + 4. * M[0][1] * M[1][0]) / t2
    if t3 >= 0.:
        print("\nfound no eigenmode solution")
        c_re = c_im = 0.
    else:
        c_im = np.sqrt(-t3)
    return complex(c_re, c_im)

#%% Tests

if __name__ == "__main__":

    # Ray propagation through optical elements

    q = qbeam(0., rayleighr(1e-1, 632.8e-6), 632.8e-6, "in")

    M = np.matmul(propagation(172.1), \
        np.matmul(thin_lens(140.4), \
        np.matmul(propagation(1400.), \
        np.matmul(thin_lens(258.), \
        propagation(400.)))))

    qp = qtrans(M, q)

    q.q_out()
    qp.q_out()

    # The output should read

    # q.in    :      z[e-3m]     w0[e-6m]
    #                      0          100
    #                R[e-3m]     wz[e-6m]
    #                    inf          100
    # q.out   :      z[e-3m]     w0[e-6m]
    #             -0.0721956      40.0079
    #                R[e-3m]     wz[e-6m]
    #               -874.728      40.0095


    # A linear resonator, left mirror is plane, distance to right mirror is
    # 400 mm, curvature of right mirror is 2 m

    qt = qbeam(0., rayleighr(1e-1, 1064e-6), 1064e-6, "res")

    M = np.matmul(propagation(400.), np.matmul(reflection(2000.), propagation(400.)))

    c = eigenmod_res(M)

    qt.z_re = c.real
    qt.zr_im = c.imag

    qt.q_out()

    # The output should read

    # q.res   :      z[e-3m]     w0[e-6m]
    #                      0      520.524
    #                R[e-3m]     wz[e-6m]
    #                    inf      520.524
