from ARTcore.ModuleDependencyConstructor import DepSolver
import numpy as np

class OAP_calculator(DepSolver):
    def __init__(self):
        self.nodes = ['fp', 'fs', 'theta', 'Rc', 'OAD', 'p', 'more_than_90']
        self.functions = {
            ('fs',): [
                (['fp', 'theta'], self.fs_fptheta)
            ],
            ('fp',): [
                (['Rc'], self.fp_Rc),
                (['fs', 'theta'], self.fp_fstheta)
            ],
            ('Rc',): [
                (['fp'], self.Rc_fp),
            ],
            ('OAD',): [
                (['fs', 'theta'], self.OAD_fstheta),
            ],
            ('theta',):[
                (['OAD', 'fp'], self.theta_OADfp),
                (['OAD', 'fs', 'more_than_90'], self.theta_OADfsm90)
            ],
            ('p',):[
                (['fs', 'theta'], lambda fs, theta: fs*(1+np.cos(np.deg2rad(theta))))
            ],
            ('more_than_90',):[
                (['theta'], self.more_than_90)
            ]
        }
        super().__init__(self.nodes, self.functions)

    def fs_fptheta(self, fp, theta):
        return 2 * fp / (1 + np.cos(np.deg2rad(theta)))

    def fp_Rc(self, Rc):
        return Rc / 2

    def Rc_fp(self, fp):
        return 2 * fp

    def OAD_fstheta(self, fs, theta):
        return fs * np.sin(np.deg2rad(theta))

    def theta_OADfp(self, OAD, fp):
        return np.rad2deg(np.arctan(OAD / (2 * fp)) * 2)

    def theta_OADfsm90(self, OAD, fs, more_than_90):
        if more_than_90:
            return 180-np.rad2deg(np.arcsin(OAD / fs))
        return np.rad2deg(np.arcsin(OAD / fs))

    def fp_fstheta(self, fs, theta):
        return (1 + np.cos(np.deg2rad(theta))) * fs / 2
    
    def more_than_90(self, theta):
        return theta > 90

class EllipsoidalMirrorCalculator(DepSolver):
    def __init__(self):
        self.nodes = ['m',
                      'f1',
                      'f2',
                      'offset_angle',
                      'incidence_angle',
                      'a',
                      'b',
                      'l',
                      'e',
                      'c']
        self.functions = {
            ('m',): [
                (['f1', 'f2'], self.m_f1f2)
            ],
            ('f1',): [
                (['f2', 'm'], self.f1_f2m),
                (['a', 'f2'], self.f1_af2),
            ],
            ('f2',): [
                (['f1', 'm'], self.f2_f1m),
                (['a', 'f1'], self.f2_af1)
            ],
            ('a',): [
                (['f1', 'f2'], self.a_f1f2),
                (['b', 'c'], self.a_bc)
            ],
            ('b',): [
                (['a', 'c'], self.b_ac),
                (['a', 'e'], self.b_ae),
            ],
            ('incidence_angle',): [
                (['f2','offset_angle', 'c'], self.incidenceangle_f2offsetc)
            ],
            ('e',): [
                (['a', 'b'], self.e_ab),
                (['a', 'c'], self.e_ac)
            ],
            ('c',): [
                (['a', 'e'], self.c_ae),
                (['a', 'b'], self.c_ab),
                (['f1', 'f2', 'offset_angle'],self.c_f1f2offset)
            ],
            ('l',):[
                (['b','a'], lambda a,b:b**2/a),
            ],
            ('offset_angle',): [
                (['c', 'f1', 'f2'], self.offset_cf1f2),
            ],
            ('f1','f2'): [
                (['a', 'c', 'offset_angle', 'm'], self.f1f2_acoffsetm),
            ]
        }
        super().__init__(self.nodes, self.functions)

    def c_ae(self, a, e):
        return a*e

    def e_ac(self, a, c):
        return c/a

    def e_ab(self, a, b):
        return np.sqrt(1-b**2/a**2)

    def b_ae(self, a, e):
        return a*np.sqrt(1 - e**2)

    def c_ab(self, a,b):
        return np.sqrt(a**2-b**2)

    def c_f1f2offset(self,f1,f2, offset_angle):
        return np.sqrt((f1**2 + f2**2 - 2*f1*f2*np.cos(np.deg2rad(offset_angle))))/2

    def m_f1f2(self, f1, f2):
        return f2 / f1

    def f1_f2m(self, f2, m):
        return f2 / m

    def f2_f1m(self, f1, m):
        return f1 * m

    def a_f1f2(self, f1, f2):
        return (f1 + f2)/2

    def f1_af2(self, a, f2):
        return 2*a-f2

    def f2_af1(self, a, f1):
        return 2*a-f1

    def b_ac(self, a, c):
        return np.sqrt(a**2-c**2)

    def a_bc(self, b, c):
        return np.sqrt(b**2 + c**2)

    def incidenceangle_f2offsetc(self, f2, offset_angle, c):
        return 180-offset_angle-np.rad2deg(np.arcsin(f2/(2*c)*np.sin(np.deg2rad(offset_angle))))

    def offset_cf1f2(self, c, f1, f2):
        return np.rad2deg(np.arccos((-(2*c)**2 + f1**2 + f2**2)/(2*f1*f2)))

    def f1f2_acoffsetm(self, a, c, offset_angle, magnification):
        delta = (2*a)**2 - 2*((2*a)**2 - (2*c)**2)/(1+np.cos(np.deg2rad(offset_angle)))
        f1 = (2*a+np.sqrt(delta))/2
        f2 = (2*a-np.sqrt(delta))/2
        if magnification>=1:
            return f2,f1
        return f1, f2

class UniformSpectraCalculator(DepSolver):
    def __init__(self):
        self.nodes = ['lambda_min',
                      'lambda_max',
                      'lambda_center',
                      'lambda_width',
                      'eV_min',
                      'eV_max',
                      'eV_center',
                      'eV_width']
        self.functions = {
            ('lambda_min',): [
                (['eV_max'], self.lambda2eV)
            ],
            ('lambda_max',): [
                (['eV_min'], self.lambda2eV)
            ],
            ('lambda_center',): [
                (['lambda_min', 'lambda_max'], self.lambda_minmaxcenter)
            ],
            ('lambda_width',): [
                (['lambda_min', 'lambda_max'], self.lambda_minmaxwidth)
            ],
            ('eV_min',): [
                (['lambda_max'], self.eV2lambda)
            ],
            ('eV_max',): [
                (['lambda_min'], self.eV2lambda)
            ],
            ('eV_center',): [
                (['eV_min', 'eV_max'], self.eV_minmaxcenter)
            ],
            ('eV_width',): [
                (['eV_min', 'eV_max'], self.eV_minmaxwidth)
            ],
            ('lambda_min', 'lambda_width'): [
                (['lambda_max', 'lambda_center'], self.lambda_maxcenterminwidth)
            ],
            ('lambda_max', 'lambda_width'): [
                (['lambda_min', 'lambda_center'], self.lambda_mincentermaxwidth)
            ],
            ('lambda_center',): [
                (['lambda_max', 'lambda_width'], lambda x1,x2: x1 - x2 / 2)
            ]
        }
        super().__init__(self.nodes, self.functions)
    def lambda_minmaxcenter(self, lambda_min, lambda_max):
        return (lambda_min + lambda_max) / 2
    def lambda_minmaxwidth(self, lambda_min, lambda_max):
        return lambda_max - lambda_min
    def lambda_centerwidth(self, lambda_center, lambda_width):
        return lambda_center - lambda_width / 2, lambda_center + lambda_width / 2
    def eV_minmaxcenter(self, eV_min, eV_max):
        return (eV_min + eV_max) / 2
    def eV_minmaxwidth(self, eV_min, eV_max):
        return eV_max - eV_min
    def eV_centerwidth(self, eV_center, eV_width):
        return eV_center - eV_width / 2, eV_center + eV_width / 2
    def eV2lambda(self, eV):
        return 1239.84193 / eV
    def lambda2eV(self, lambda_):
        return 1239.84193 / lambda_
    def lambda_maxcenterminwidth(self, lambda_max, lambda_center):
        return 2 * lambda_center - lambda_max, 2*(lambda_max - lambda_center)
    def lambda_mincentermaxwidth(self, lambda_min, lambda_center):
        return 2 * lambda_center - lambda_min, 2*(lambda_center - lambda_min)