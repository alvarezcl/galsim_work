import numpy as np
from scipy.optimize import newton
from scipy.special import gammainc, gamma

def lazyprop(fn):
    attr_name = '_lazy_' + fn.__name__
    @property
    def _lazyprop(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazyprop

class Ellipticity(object):
    def __init__(self, **kwargs):
        # always immediately set e1 and e2, which are normal class attributes
        # everything else is a lazy property
        if len(kwargs) == 0:
            self.e1 = 0.0
            self.e2 = 0.0
        if 'e' in kwargs:
            self._lazy_e = kwargs.pop('e')
            self._lazy_theta = kwargs.pop('theta', 0.0)
            self.e1 = self.e * np.cos(self.theta)
            self.e2 = self.e * np.sin(self.theta)
        elif 'g' in kwargs:
            self._lazy_g = kwargs.pop('g')
            self._lazy_theta = kwargs.pop('theta', 0.0)
            self._lazy_e = 2.0*self.g / (1 + self.g**2)
            self.e1 = self.e * np.cos(self.theta)
            self.e2 = self.e * np.sin(self.theta)
        elif 'eta' in kwargs:
            self._lazy_eta = kwargs.pop('eta')
            self._lazy_theta = kwargs.pop('theta', 0.0)
            self._lazy_e = np.tanh(self.eta/2.0)
            self.e1 = self.e * np.cos(self.theta)
            self.e2 = self.e * np.sin(self.theta)
        elif 'q' in kwargs:
            self._lazy_q = kwargs.pop('q')
            self._lazy_theta = kwargs.pop('theta', 0.0)
            self._lazy_e = (1.0 - self.q)/(1.0 + self.q)
            self.e1 = self.e * np.cos(self.theta)
            self.e2 = self.e * np.sin(self.theta)
        elif 'e1' in kwargs:
            self.e1 = kwargs.pop('e1')
            self.e2 = kwargs.pop('e2', 0.0)
            self._lazy_e = np.sqrt(self.e1**2 + self.e2**2)
            self._lazy_theta = np.arctan2(self.e2, self.e1)
        elif 'g1' in kwargs:
            self._lazy_g1 = kwargs.pop('g1')
            self._lazy_g2 = kwargs.pop('g2', 0.0)
            self._lazy_g = np.sqrt(self.g1**2 + self.g2**2)
            self._lazy_e = 2.0*self.g / (1 + self.g**2)
            self._lazy_theta = np.arctan2(self.g2, self.g1)
            self.e1 = self.e * np.cos(self.theta)
            self.e2 = self.e * np.sin(self.theta)
        elif 'eta1' in kwargs:
            self._lazy_eta1 = kwargs.pop('eta1')
            self._lazy_eta2 = kwargs.pop('eta2', 0.0)
            self._lazy_eta = np.sqrt(self.eta1**2 + self.eta2**2)
            self._lazy_e = np.tanh(self.eta/2.0)
            self._lazy_theta = np.arctan2(self.eta2, self.eta1)
            self.e1 = self.e * np.cos(self.theta)
            self.e2 = self.e * np.sin(self.theta)
        if 'theta' in kwargs: # bare theta passed
            self._lazy_theta = kwargs.pop('theta')
            self.e1 = 0.0
            self.e2 = 0.0
        if kwargs:
            raise ValueError('Too many arguments to Ellipticity')

    @lazyprop
    def e(self):
        return np.sqrt(self.e1**2 + self.e2**2)

    @lazyprop
    def g(self):
        return self.e / (1.0 + np.sqrt(1.0 - self.e**2))

    @lazyprop
    def eta(self):
        return 2.0 * np.arctanh(self.e)

    @lazyprop
    def q(self):
        return (1.0 - self.e)/(1.0 + self.e)

    @lazyprop
    def theta(self):
        return np.arctan2(self.e2, self.e1)

    @lazyprop
    def g1(self):
        return self.e1 / (1.0 + np.sqrt(1.0 - self.e**2))

    @lazyprop
    def g2(self):
        return self.e2 / (1.0 + np.sqrt(1.0 - self.e**2))

    @lazyprop
    def eta1(self):
        return self.eta * np.cos(self.theta)

    @lazyprop
    def eta2(self):
        return self.eta * np.sin(self.theta)


class Profile(object):
    def __init__(self, **kwargs):
        # centroid attributes
        self.x0 = kwargs.pop('x0', 0.0)
        self.y0 = kwargs.pop('y0', 0.0)

        # both flux and peak are lazy properties.
        if 'peak' in kwargs:
            self._lazy_peak = kwargs.pop('peak')
        else:
            self._lazy_flux = kwargs.pop('flux', 1.0)

        # Here's the tricky part.  There are many ways to parameterize the ellipticity+size
        # of the profile.  We'll make a, b, phi attributes and everything else lazy properties.
        if 'a' in kwargs and 'b' in kwargs:
            self.a = kwargs.pop('a')
            self.b = kwargs.pop('b')
            self.phi = kwargs.pop('phi', 0.0)
        elif 'C' in kwargs: # inverse covariance matrix analogue
            self._lazy_C = kwargs.pop('C')
            one_over_a_squared = 0.5 * (self.C[0,0] + self.C[1,1]
                                        + np.sqrt((self.C[0,0] - self.C[1,1])**2
                                                  + 4.0 * self.C[0,1]**2))
            one_over_b_squared = self.C[0,0] + self.C[1,1] - one_over_a_squared
            # there's degeneracy between a, b and phi at this point so enforce a > b
            if one_over_a_squared > one_over_b_squared:
                one_over_a_squared, one_over_b_squared = one_over_b_squared, one_over_a_squared
            self.a = np.sqrt(1.0 / one_over_a_squared)
            self.b = np.sqrt(1.0 / one_over_b_squared)
            if self.a == self.b:
                self.phi = 0.0
            else:
                self.phi = 0.5 * np.arctan2(
                    2.0 * self.C[0,1] / (one_over_a_squared - one_over_b_squared),
                    (self.C[0,0] - self.C[1,1]) / (one_over_a_squared - one_over_b_squared))
        elif 'I' in kwargs: # second moment matrix
            self._lazy_I = kwargs.pop('I')
            self.phi = 0.5 * np.arctan2(self.I[0,0] - self.I[1,1], 2.0 * self.I[0,1])
            self._lazy_r2 = (4.0*np.linalg.det(self.I))**(0.25)
            e = (np.sqrt((self.I[0,0] - self.I[1,1])**2 + 4.0*self.I[0,1]**2)
                 /(self.I[0,0] + self.I[1,1]))
            q = np.sqrt((1.0-e)/(1.0+e))
            self.a = self.half_light_radius / np.sqrt(q)
            self.b = self.half_light_radius * np.sqrt(q)
        else:
            if 'half_light_radius' in kwargs:
                self._lazy_half_light_radius = kwargs.pop('half_light_radius')
            elif 'FWHM' in kwargs:
                self._lazy_FWHM = kwargs.pop('FWHM')
            elif 'r2' in kwargs:
                self._lazy_r2 = kwargs.pop('r2')
            else:
                raise ValueError('Size not specified in Sersic')
            if 'phi' in kwargs:
                self.phi = kwargs.pop('phi')
                theta = 2.0*self.phi
                self._lazy_ellip = Ellipticity(theta=theta, **kwargs)
            else:
                self._lazy_ellip = Ellipticity(**kwargs)
                self.phi = 0.5*self.ellip.theta
            self.a = self.half_light_radius / np.sqrt(self.ellip.q)
            self.b = self.half_light_radius / np.sqrt(self.ellip.q)

    @lazyprop
    def half_light_radius(self):
        return np.sqrt(self.a * self.b)

    @lazyprop
    def FWHM(self):
        raise NotImplementedError('Derived class must define FWHM property getter')

    @lazyprop
    def r2(self):
        raise NotImplementedError('Derived class must define r2 property getter')

    @lazyprop
    def ellip(self):
        return Ellipticity(q=self.b/self.a, theta=2.0*self.phi)

    @lazyprop
    def peak(self):
        raise NotImplementedError('Derived class must define peak property getter')

    @lazyprop
    def flux(self):
        raise NotImplementedError('Derived class must define flux property getter')

    @lazyprop
    def C(self):
        C = np.matrix(np.identity(2, dtype=float))
        cph = np.cos(self.phi)
        sph = np.sin(self.phi)
        C[0,0] = (cph/self.a)**2 + (sph/self.b)**2
        C[0,1] = C[1,0] = 0.5 * (1.0/self.a**2 - 1.0/self.b**2) * np.sin(2.0 * self.phi)
        C[1,1] = (sph/self.a)**2 + (cph/self.b)**2
        return C

    @lazyprop
    def I(self):
        unrotI = 0.5 * self.r2**2 * np.matrix([[1./self.ellip.q, 0.0],
                                               [0.0, self.ellip.q]], dtype=float)
        cph = np.cos(self.phi)
        sph = np.sin(self.phi)
        R = np.matrix([[cph, -sph],
                       [sph, cph]], dtype=float)
        return R * unrotI * R.T

class Sersic(Profile):
    def __init__(self, **kwargs):
        # Sersic index
        if 'n' not in kwargs:
            raise ValueError('Sersic initialization requires `n` index')
        self.n = kwargs.pop('n')
        super(Sersic, self).__init__(**kwargs)

    @lazyprop
    def half_light_radius(self):
        if hasattr(self, 'a') and hasattr(self, 'b'):
            return super(Sersic, self).half_light_radius
        elif hasattr(self, '_lazy_r2'):
            return self.r2 / self.HLR_to_r2(1.0, self.n)
        elif hasattr(self, '_lazy_FWHM'):
            return self.FWHM / self.HLR_to_FWHM(1.0, self.n, self.kappa)

    @lazyprop
    def kappa(self):
        return self.n_to_kappa(self.n)

    @lazyprop
    def FWHM(self):
        return self.HLR_to_FWHM(self.half_light_radius, self.n, self.kappa)

    @lazyprop
    def r2(self):
        return self.HLR_to_r2(self.half_light_radius, self.n)

    @lazyprop
    def peak(self):
        return self.flux_to_peak(self.flux, self.n, self.half_light_radius, self.kappa)

    @lazyprop
    def flux(self):
        return self.peak_to_flux(self.peak, self.n, self.half_light_radius, self.kappa)

    def __call__(self, x, y):
        x = np.array(x)
        y = np.array(y)
        xvec = np.array([x-self.x0, y-self.y0]).T
        if len(xvec.shape) == 1:
            exponent = np.einsum('j,jk,k', xvec, self.C, xvec)
        else:
            exponent = np.einsum('ij,jk,ij->i', xvec, self.C, xvec)
        exponent **= 0.5 / self.n
        exponent *= -self.kappa
        return self.peak * np.exp(exponent)

    @staticmethod
    def n_to_kappa(n):
        '''Compute Sersic exponent factor kappa from the Sersic index'''
        kguess = 1.9992 * n - 0.3271
        return newton(lambda k: gammainc(2.0 * n, k) - 0.5, kguess)

    @classmethod
    def HLR_to_FWHM(cls, half_light_radius, n, kappa=None):
        '''Compute the full-width at half maximum given input parameters.'''
        if kappa is None:
            kappa = cls.n_to_kappa(n)
        return 2.0 * half_light_radius * (np.log(2.0) / kappa)**n

    @classmethod
    def FWHM_to_HLR(cls, FWHM, n, kappa=None):
        if kappa is None:
            kappa = cls.n_to_kappa(n)
        return 0.5 * FWHM * (kappa / np.log(2.0))**n

    @staticmethod
    def HLR_to_r2(half_light_radius, n):
        """ Using Mathematica-generated fitting function, compute the circularized second moment
        square radius sqrt(r^2) = sqrt(Ixx + Iyy); where Ixx = Iyy and Ixy = 0.0
        """
        return half_light_radius * (0.985444 + n * (0.391016 + n * (0.0739602
                                    + n * (0.00698719 + n * (0.00212432
                                    + n * (-0.000154052 + n * 0.0000219632))))))

    @staticmethod
    def r2_to_HLR(cls, r2, n):
        return r2 / self.HLR_to_r2(1.0, n)

    @classmethod
    def flux_to_peak(cls, flux, n, half_light_radius, kappa=None):
        '''Compute peak surface brightness from integrated flux and Sersic params.'''
        if kappa is None:
            kappa = cls.n_to_kappa(n)
        fluxnorm = cls.peak_to_flux(1.0, n, half_light_radius, kappa=kappa)
        return flux/fluxnorm

    @classmethod
    def peak_to_flux(cls, peak, n, half_light_radius, kappa=None):
        '''Compute flux integrated over all space given input parameters.'''
        if kappa is None:
            kappa = cls.n_to_kappa(n)
        return (2.0 * np.pi * peak * n * (kappa**(-2.0 * n)) * (half_light_radius**2.0)
                * gamma(2.0 * n))

class Moffat(Profile):
    def __init__(self, **kwargs):
        # Sersic index
        if 'beta' not in kwargs:
            raise ValueError('Moffat initialization requires `beta` index')
        self.beta = kwargs.pop('beta')
        super(Moffat, self).__init__(**kwargs)

    @lazyprop
    def half_light_radius(self):
        if hasattr(self, 'a') and hasattr(self, 'b'):
            return super(Moffat, self).half_light_radius
        elif hasattr(self, '_lazy_r2'):
            return self.r2 / self.HLR_to_r2(1.0, self.beta)
        elif hasattr(self, '_lazy_FWHM'):
            return self.FWHM / self.HLR_to_FWHM(1.0, self.beta)

    @lazyprop
    def FWHM(self):
        return self.HLR_to_FWHM(self.half_light_radius, self.beta)

    @lazyprop
    def r2(self):
        return self.HLR_to_r2(self.half_light_radius, self.beta)

    @lazyprop
    def peak(self):
        return self.flux_to_peak(self.flux, self.beta, self.C)

    @lazyprop
    def flux(self):
        return self.peak_to_flux(self.peak, self.beta, self.C)

    def __call__(self, x, y):
        x = np.array(x)
        y = np.array(y)
        xvec = np.array([x-self.x0, y-self.y0]).T
        if len(xvec.shape) == 1:
            base = np.einsum('j,jk,k', xvec, self.C, xvec)
        else:
            base = np.einsum('ij,jk,ij->i', xvec, self.C, xvec)
        return self.peak * base**(-self.beta)

    @staticmethod
    def HLR_to_FWHM(half_light_radius, beta):
        return 2.0 * half_light_radius * np.sqrt(
            (2.0**(1.0/beta) - 1.0) * (2.0**(1.0/(beta-1.0)) - 1.0))

    @staticmethod
    def FWHM_to_HLR(FWHM, beta):
        return 0.5 * FWHM / np.sqrt(
            (2.0**(1.0/beta) - 1.0) * (2.0**(1.0/(beta-1.0)) - 1.0))

    @classmethod
    def HLR_to_r2(cls, half_light_radius, beta):
        return 0.5 * cls.HLR_to_FWHM(half_light_radius, beta) / np.sqrt(
            (2.0**(1.0/beta) - 1.0) * (beta - 2.0))

    @classmethod
    def r2_to_HLR(cls, r2, beta):
        return cls.FWHM_to_HLR(2.0 * r2 * np.sqrt((2.0**(1.0/beta) - 1.0)(beta - 2.0)), beta)

    @staticmethod
    def flux_to_peak(flux, beta, C):
        return flux * (beta - 1.0) / (np.pi / np.sqrt(np.abs(np.linalg.det(C))))

    @staticmethod
    def peak_to_flux(peak, beta, C):
        return peak / (beta - 1.0) * (np.pi / np.sqrt(np.abs(np.linalg.det(C))))

def compare_Ellipticities(a, b):
    np.testing.assert_almost_equal(a.e, b.e)
    np.testing.assert_almost_equal(a.e1, b.e1)
    np.testing.assert_almost_equal(a.e2, b.e2)
    np.testing.assert_almost_equal(a.g, b.g)
    np.testing.assert_almost_equal(a.g1, b.g1)
    np.testing.assert_almost_equal(a.g2, b.g2)
    np.testing.assert_almost_equal(a.eta, b.eta)
    np.testing.assert_almost_equal(a.eta1, b.eta1)
    np.testing.assert_almost_equal(a.eta2, b.eta2)
    np.testing.assert_almost_equal(a.theta, b.theta)
    np.testing.assert_almost_equal(a.q, b.q)

def test_Ellipticity():
    a = Ellipticity(e1=0.11, e2=-0.023)

    b = Ellipticity(e=a.e, theta=a.theta)
    compare_Ellipticities(a, b)

    b = Ellipticity(g=a.g, theta=a.theta)
    compare_Ellipticities(a, b)

    b = Ellipticity(eta=a.eta, theta=a.theta)
    compare_Ellipticities(a, b)

    b = Ellipticity(q=a.q, theta=a.theta)
    compare_Ellipticities(a, b)

    b = Ellipticity(e1=a.e1, e2=a.e2)
    compare_Ellipticities(a, b)

    b = Ellipticity(g1=a.g1, g2=a.g2)
    compare_Ellipticities(a, b)

    b = Ellipticity(eta1=a.eta1, eta2=a.eta2)
    compare_Ellipticities(a, b)

def compare_Sersics(a, b):
    np.testing.assert_almost_equal(a.x0, b.x0)
    np.testing.assert_almost_equal(a.y0, b.y0)
    np.testing.assert_almost_equal(a.n, b.n)
    np.testing.assert_almost_equal(a.kappa, b.kappa)
    np.testing.assert_almost_equal(a.a, b.a)
    np.testing.assert_almost_equal(a.b, b.b)
    np.testing.assert_array_almost_equal(a.I, b.I)
    np.testing.assert_array_almost_equal(a.C, b.C)
    np.testing.assert_almost_equal(a.ellip.e, b.ellip.e)
    np.testing.assert_almost_equal(a.ellip.e1, b.ellip.e1)
    np.testing.assert_almost_equal(a.ellip.e2, b.ellip.e2)
    np.testing.assert_almost_equal(a.ellip.g, b.ellip.g)
    np.testing.assert_almost_equal(a.ellip.g1, b.ellip.g1)
    np.testing.assert_almost_equal(a.ellip.g2, b.ellip.g2)
    np.testing.assert_almost_equal(a.ellip.eta, b.ellip.eta)
    np.testing.assert_almost_equal(a.ellip.eta1, b.ellip.eta1)
    np.testing.assert_almost_equal(a.ellip.eta2, b.ellip.eta2)
    np.testing.assert_almost_equal(a.ellip.theta, b.ellip.theta)
    np.testing.assert_almost_equal(a.ellip.q, b.ellip.q)
    np.testing.assert_almost_equal(a.FWHM, b.FWHM)
    np.testing.assert_almost_equal(a.r2, b.r2)
    np.testing.assert_almost_equal(a.half_light_radius, b.half_light_radius)
    np.testing.assert_almost_equal(a.flux, b.flux)
    np.testing.assert_almost_equal(a.peak, b.peak)

def test_Sersic():
    # test radial part
    a = Sersic(n=2.2, FWHM=1.1)

    b = Sersic(n=2.2, r2=a.r2)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=a.FWHM)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, half_light_radius=a.half_light_radius)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, I=a.I)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, C=a.C)
    compare_Sersics(a, b)

    # test flux part
    a = Sersic(n=2.2, FWHM=1.0, flux=1.1)
    b = Sersic(n=2.2, FWHM=1.0, peak=a.peak)
    compare_Sersics(a, b)

    # test ellipticity part
    a = Sersic(n=2.2, FWHM=1.0, e1=0.21, e2=-0.034)
    b = Sersic(n=2.2, FWHM=1.0, eta1=a.ellip.eta1, eta2=a.ellip.eta2)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, g1=a.ellip.g1, g2=a.ellip.g2)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, e1=a.ellip.e1, e2=a.ellip.e2)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, e=a.ellip.e, theta=a.ellip.theta)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, g=a.ellip.g, theta=a.ellip.theta)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, eta=a.ellip.eta, theta=a.ellip.theta)
    compare_Sersics(a, b)

    b = Sersic(n=2.2, FWHM=1.0, q=a.ellip.q, theta=a.ellip.theta)
    compare_Sersics(a, b)

def compare_Moffats(a, b):
    np.testing.assert_almost_equal(a.x0, b.x0)
    np.testing.assert_almost_equal(a.y0, b.y0)
    np.testing.assert_almost_equal(a.beta, b.beta)
    np.testing.assert_almost_equal(a.a, b.a)
    np.testing.assert_almost_equal(a.b, b.b)
    np.testing.assert_array_almost_equal(a.I, b.I)
    np.testing.assert_array_almost_equal(a.C, b.C)
    np.testing.assert_almost_equal(a.ellip.e, b.ellip.e)
    np.testing.assert_almost_equal(a.ellip.e1, b.ellip.e1)
    np.testing.assert_almost_equal(a.ellip.e2, b.ellip.e2)
    np.testing.assert_almost_equal(a.ellip.g, b.ellip.g)
    np.testing.assert_almost_equal(a.ellip.g1, b.ellip.g1)
    np.testing.assert_almost_equal(a.ellip.g2, b.ellip.g2)
    np.testing.assert_almost_equal(a.ellip.eta, b.ellip.eta)
    np.testing.assert_almost_equal(a.ellip.eta1, b.ellip.eta1)
    np.testing.assert_almost_equal(a.ellip.eta2, b.ellip.eta2)
    np.testing.assert_almost_equal(a.ellip.theta, b.ellip.theta)
    np.testing.assert_almost_equal(a.ellip.q, b.ellip.q)
    np.testing.assert_almost_equal(a.FWHM, b.FWHM)
    np.testing.assert_almost_equal(a.r2, b.r2)
    np.testing.assert_almost_equal(a.half_light_radius, b.half_light_radius)
    np.testing.assert_almost_equal(a.flux, b.flux)
    np.testing.assert_almost_equal(a.peak, b.peak)

def test_Moffat():
    # test radial part
    a = Moffat(beta=3.0, FWHM=1.1)

    b = Moffat(beta=3.0, r2=a.r2)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=a.FWHM)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, half_light_radius=a.half_light_radius)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, I=a.I)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, C=a.C)
    compare_Moffats(a, b)

    # test flux part
    a = Moffat(beta=3.0, FWHM=1.0, flux=1.1)
    b = Moffat(beta=3.0, FWHM=1.0, peak=a.peak)
    compare_Moffats(a, b)

    # test ellipticity part
    a = Moffat(beta=3.0, FWHM=1.0, e1=0.21, e2=-0.034)
    b = Moffat(beta=3.0, FWHM=1.0, eta1=a.ellip.eta1, eta2=a.ellip.eta2)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, g1=a.ellip.g1, g2=a.ellip.g2)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, e1=a.ellip.e1, e2=a.ellip.e2)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, e=a.ellip.e, theta=a.ellip.theta)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, g=a.ellip.g, theta=a.ellip.theta)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, eta=a.ellip.eta, theta=a.ellip.theta)
    compare_Moffats(a, b)

    b = Moffat(beta=3.0, FWHM=1.0, q=a.ellip.q, theta=a.ellip.theta)
    compare_Moffats(a, b)

if __name__ == '__main__':
    test_Ellipticity()
    test_Sersic()
    test_Moffat()
