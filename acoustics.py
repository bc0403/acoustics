"""
This is a Python module for analysis of acoustic fields and waves
in solids, using sympy and numpy for symbolic and numeric calculation,
respectively.

References:
[1] B.A. Auld, Acoustic fields and waves in solids, Vol. I,
John Wiley & Sons, New York, 1973.
======================================================================
Created by Hao JIN on 2019/05/22
"""

from sympy import Matrix, cos, sin, Rational, pprint


class ElasticMaterial:
    """
    Elastic materials used in acoustic waves and fields.

    Rotational information
    1. Rotate axes (passive transformation), not a vector (active
    transformation). A vector in the new axes v' can be represented
    by the rotational maxtix (R) and its components in old axes (v): v' = R*v.
    2. Right-hand rule: the rotation angle is defined to be positive for a
    rotation that is counter-clockwise when viewed by an observer looking
    along the rotation axis towards the origin.
    3. refer to
    https://www.mathworks.com/help/phased/ref/roty.html
    https://en.wikipedia.org/wiki/Euler_angles
    https://en.wikipedia.org/wiki/Active_and_passive_transformation
    """

    def __init__(self, density, stiffness, epsilon):
        """Initialize attributes to describe a elastic materials"""
        self.density = density
        self.stiffness = stiffness
        self.epsilon = epsilon

    def get_description(self):
        """Return the description of material"""
        print("This is a material with density, stiffness, and epsilon as:")
        pprint(self.density)
        pprint(self.stiffness)
        pprint(self.epsilon)

    def rotx_R(self, theta):
        """
        rotate of theta about the x-axis
        return the rotational matrix R
        """
        # refer to pages 19-23 of Auld's book
        return Matrix([
            [1, 0, 0],
            [0, cos(theta), sin(theta)],
            [0, -sin(theta), cos(theta)]
        ])

    def roty_R(self, theta):
        """
        rotate of theta about the y-axis
        return the rotational matrix R
        """
        return Matrix([
            [cos(theta), 0, -sin(theta)],
            [0, 1, 0],
            [sin(theta), 0, cos(theta)]
        ])

    def rotz_R(self, theta):
        """
        rotate of theta about the z-axis
        return the rotational matrix R
        """
        return Matrix([
            [cos(theta), sin(theta), 0],
            [-sin(theta), cos(theta), 0],
            [0, 0, 1]
        ])

    def rot_M(self, R):
        """
        return the transformation matrix M
        according to the rotational matrix R
        """
        # refer to page 74 of Auld's book
        R00 = R[0, 0]
        R01 = R[0, 1]
        R02 = R[0, 2]
        R10 = R[1, 0]
        R11 = R[1, 1]
        R12 = R[1, 2]
        R20 = R[2, 0]
        R21 = R[2, 1]
        R22 = R[2, 2]

        M = Matrix([
            [R00**2, R01**2, R02**2, 2*R01*R02, 2*R02*R00, 2*R00*R01],
            [R10**2, R11**2, R12**2, 2*R11*R12, 2*R12*R10, 2*R10*R11],
            [R20**2, R21**2, R22**2, 2*R21*R22, 2*R22*R20, 2*R20*R21],
            [R10*R20, R11*R21, R12*R22, R11*R22 + R12*R21, R10*R22 + R12*R20,
             R10*R21 + R11*R20],
            [R20*R00, R21*R01, R22*R02, R01*R22 + R02*R21, R00*R22 + R02*R20,
             R00*R21 + R01*R20],
            [R00*R10, R01*R11, R02*R12, R01*R12 + R02*R11, R00*R12 + R02*R10,
             R00*R11 + R01*R10]
        ])
        return M

    def rot_euler_RM(self, alpha, beta, gamma):
        """
        rotate by Z(alpha)-X(beta)-Z(gamma) euler angles,
        return the rotational matrix R and transformation matrix M
        (R, M)
        """
        R = self.rotz_R(gamma)*self.rotx_R(beta)*self.rotz_R(alpha)
        M = self.rot_M(R)
        return (R, M)

    def rot_update(self, R, M):
        """
        update matrices of stiffness and epsilon by the rotational matrix R and
        transformation matrix M, refer to p76 and p117 of Auld's book
        """
        self.stiffness = M*self.stiffness*M.T
        self.epsilon = R*self.epsilon*R.T

    def rot_euler_update(self, alpha, beta, gamma):
        """
        rotate by Z(alpha)-X(beta)-Z(gamma) euler angles,
        update matrices of stiffness and epsilon
        """
        R, M = self.rot_euler_RM(alpha, beta, gamma)
        self.rot_update(R, M)

    def cal_Gamma(self, li, lj, liK, lLj):
        """
        calculate the Christoffel matrix,
        li, the direction of wave propagation, [lx, ly, lz]
        lj, transpose of li
        liK = Matrix([
            [lx, 0, 0, 0, lz, ly],
            [0, ly, 0, lz, 0, lx],
            [0, 0, lz, ly, lx, 0]
        ])
        lLj = transpose of liK
        """
        # refer to pages 164-165, Auld's book
        return liK*self.stiffness*lLj


class PiezoMaterial(ElasticMaterial):
    """Class for piezoelectric materials."""

    def __init__(self, density, stiffness, epsilon, piezoelec):
        """Initialize attributes of the parent class."""
        super().__init__(density, stiffness, epsilon)
        self.piezoelec = piezoelec

    def get_description(self):
        """Return the description of material"""
        super().get_description()
        print("And with piezoelectric stress constants as:")
        pprint(self.piezoelec)

    def rot_update(self, R, M):
        """
        update matrices of stiffness, epsilon and piezoelectric stress
        constants by the rotational matrix R and transformation matrix M,
        refer to p76, p117, and p275 of Auld's book
        """
        super().rot_update(R, M)
        self.piezoelec = R*self.piezoelec*M.T
    
    def rot_euler_update(self, alpha, beta, gamma):
        """
        rotate by Z(alpha)-X(beta)-Z(gamma) euler angles,
        update matrices of stiffness and epsilon
        """
        R, M = super().rot_euler_RM(alpha, beta, gamma)
        self.rot_update(R, M)

    def cal_cD(self, li, lj):
        """
        calculate cD, the stiffness constants at zero electric displacement,
        li, the direction of wave propagation, [lx, ly, lz]
        lj, transpose of li
        """
        # refer to page 300, Auld's book
        cD = self.stiffness + \
            (self.piezoelec.T*lj)*(li*self.piezoelec) / \
            (li*self.epsilon*lj)[0, 0]
        return cD

    def cal_Gamma(self, li, lj, liK, lLj):
        """
        calculate the Christoffel matrix,
        li, the direction of wave propagation, [lx, ly, lz]
        lj, transpose of li
        liK = Matrix([
            [lx, 0, 0, 0, lz, ly],
            [0, ly, 0, lz, 0, lx],
            [0, 0, lz, ly, lx, 0]
        ])
        lLj = transpose of liK
        """
        # refer to pages 164-165, 300, Auld's book
        return liK*self.cal_cD(li, lj)*lLj


class Isotropic():
    """
    acoustic properties for isotropic materials
    """

    # refer to page s 363, 379, Auld's book
    # another names: c44 ~ nu; c12 ~ lambda
    # relation: c12 = c11 - 2*c44
    def __init__(self, c11, c44, eSxx):
        self.c = Matrix([
            [c11, c11 - 2*c44, c11 - 2*c44, 0, 0, 0],
            [c11 - 2*c44, c11, c11 - 2*c44, 0, 0, 0],
            [c11 - 2*c44, c11 - 2*c44, c11, 0, 0, 0],
            [0, 0, 0, c44, 0, 0],
            [0, 0, 0, 0, c44, 0],
            [0, 0, 0, 0, 0, c44]
        ])
        self.eS = Matrix([
            [eSxx, 0, 0],
            [0, eSxx, 0],
            [0, 0, eSxx]
        ])


class Cubic():
    """
    acoustic properties for Cubic materials
    Al, Au, Ag, Ni, W
    """

    # refer to pages 362, 374, 379, Auld's book
    def __init__(self, c11, c12, c44, ex4, eSxx):
        self.c = Matrix([
            [c11, c12, c12, 0, 0, 0],
            [c12, c11, c12, 0, 0, 0],
            [c12, c12, c11, 0, 0, 0],
            [0, 0, 0, c44, 0, 0],
            [0, 0, 0, 0, c44, 0],
            [0, 0, 0, 0, 0, c44]
        ])
        # for Cubic 23 and -43m
        self.e = Matrix([
            [0, 0, 0, ex4, 0, 0],
            [0, 0, 0, 0, ex4, 0],
            [0, 0, 0, 0, 0, ex4]
        ])
        self.eS = Matrix([
            [eSxx, 0, 0],
            [0, eSxx, 0],
            [0, 0, eSxx]
        ])


class Trig3m():
    """
    acoustic properties for Trig. 3m materials
    LiNbO_3, LiTaO_3
    """

    # refer to pages 362, 373, 379, Auld's book
    def __init__(self, c11, c12, c13, c14, c33, c44,
                 ex5, ey2, ez1, ez3, eSxx, eSzz):
        self.c = Matrix([
            [c11, c12, c13, c14, 0, 0],
            [c12, c11, c13, -c14, 0, 0],
            [c13, c13, c33, 0, 0, 0],
            [c14, -c14, 0, c44, 0, 0],
            [0, 0, 0, 0, c44, c14],
            [0, 0, 0, 0, c14, Rational(1, 2)*(c11 - c12)]
        ])
        self.e = Matrix([
            [0, 0, 0, 0, ex5, -ey2],
            [-ey2, ey2, 0, ex5, 0, 0],
            [ez1, ez1, ez3, 0, 0, 0]
        ])
        self.eS = Matrix([
            [eSxx, 0, 0],
            [0, eSxx, 0],
            [0, 0, eSzz]
        ])


class Trig32():
    """
    acoustic properties for Trig. 32 materials
    Quartz
    """

    # refer to pages 362, 373, 379, Auld's book
    def __init__(self, c11, c12, c13, c14, c33, c44,
                 ex1, ex4, eSxx, eSzz):
        self.c = Matrix([
            [c11, c12, c13, c14, 0, 0],
            [c12, c11, c13, -c14, 0, 0],
            [c13, c13, c33, 0, 0, 0],
            [c14, -c14, 0, c44, 0, 0],
            [0, 0, 0, 0, c44, c14],
            [0, 0, 0, 0, c14, Rational(1, 2)*(c11 - c12)]
        ])
        self.e = Matrix([
            [ex1, -ex1, 0, ex4, 0, 0],
            [0, 0, 0, 0, -ex4, -ex1],
            [0, 0, 0, 0, 0, 0]
        ])
        self.eS = Matrix([
            [eSxx, 0, 0],
            [0, eSxx, 0],
            [0, 0, eSzz]
        ])


class Hex6mm():
    """
    acoustic properties for Hex. 6mm materials
    AlN, ZnO
    """

    # refer to pages 362, 373, 379, Auld's book
    def __init__(self, c11, c12, c13, c33, c44, ex5, ez1, ez3, eSxx, eSzz):
        self.c = Matrix([
            [c11, c12, c13, 0, 0, 0],
            [c12, c11, c13, 0, 0, 0],
            [c13, c13, c33, 0, 0, 0],
            [0, 0, 0, c44, 0, 0],
            [0, 0, 0, 0, c44, 0],
            [0, 0, 0, 0, 0, Rational(1, 2)*(c11 - c12)]
        ])
        self.e = Matrix([
            [0, 0, 0, 0, ex5, 0],
            [0, 0, 0, ex5, 0, 0],
            [ez1, ez1, ez3, 0, 0, 0]
        ])
        self.eS = Matrix([
            [eSxx, 0, 0],
            [0, eSxx, 0],
            [0, 0, eSzz]
        ])


# Define function for evaluation of kt2
def eval_kt2(va, va0):
    """
    evaluate kt2 from phase velocity.
    va: phase velocity with piezoelectric effect
    va0: phase velocity without piezoelectric effect
    """
    K2 = (va/va0)**2 - 1
    kt2 = K2/(1 + K2)
    return kt2


# Define physical constants and material properties
epsilon_0 = 8.854e-12  # permittivity of free-space, F/m

# === piezoelectric materials ===
# data from B.A. Auld's book, Appendix 2
LN_auld = {  # Trig. 3m
    "rho": 4700e3,  # g/m^3
    "c11": 203e9,  # Pa
    "c12": 53e9,
    "c13": 75e9,
    "c14": 9e9,
    "c33": 245e9,
    "c44": 60e9,
    "ex5": 3.7,  # C/m^2
    "ey2": 2.5,
    "ez1": 0.2,
    "ez3": 1.3,
    "eSxx": 44*epsilon_0,  # F/m
    "eSzz": 29*epsilon_0
}

# data form Comsol Multiphysics 5.4
LN_comsol = {  # Trig. 3m
    "rho": 4700e3,  # g/m^3
    "c11": 202.897e9,  # Pa
    "c12": 52.9177e9,
    "c13": 74.9098e9,
    "c14": 8.99874e9,
    "c33": 243.075e9,
    "c44": 59.9034e9,
    "ex5": 3.69594,  # C/m^2
    "ey2": 2.53764,
    "ez1": 0.193364,
    "ez3": 1.30863,
    "eSxx": 43.6*epsilon_0,  # F/m
    "eSzz": 29.16*epsilon_0
}

# data from B.A. Auld's book, Appendix 2
Quartz_auld = {  # Trig. 32, eSxx, eSzz
    "rho": 2651e3,  # g/m^3
    "c11": 86.74e9,  # Pa
    "c12": 6.99e9,
    "c13": 11.91e9,
    "c14": -17.91e9,
    "c33": 107.2e9,
    "c44": 57.94e9,
    "ex1": 0.171,  # C/m^2
    "ex4": -0.0436,
    "eSxx": 4.5*epsilon_0,
    "eSzz": 4.6*epsilon_0
}

# data from Comsol Multiphysics 5.6
Quartz_LH_1949 = {  # same as Quartz_auld
    "rho": 2561e3,  # g/m^3
    "c11": 8.67362e10,  # Pa
    "c12": 6.98527e9,
    "c13": 1.19104e10,
    "c14": -1.79081e10,
    "c33": 1.07194e11,
    "c44": 5.79428e10,
    "ex1": 0.171,  # C/m^2
    "ex4": -0.0406,
    "eSxx": 4.428*epsilon_0,
    "eSzz": 4.634*epsilon_0
}

Quartz_RH_1949 = {
    "rho": 2561e3,  # g/m^3
    "c11": 8.67362e10,  # Pa
    "c12": 6.98527e9,
    "c13": 1.19104e10,
    "c14": -1.79081e10,
    "c33": 1.07194e11,
    "c44": 5.79428e10,
    "ex1": -0.171,  # C/m^2
    "ex4": 0.0406,
    "eSxx": 4.428*epsilon_0,
    "eSzz": 4.634*epsilon_0
}

Quartz_LH_1978 = {
    "rho": 2561e3,  # g/m^3
    "c11": 8.67362e10,  # Pa
    "c12": 6.98527e9,
    "c13": 1.19104e10,
    "c14": 1.79081e10,
    "c33": 1.07194e11,
    "c44": 5.79428e10,
    "ex1": -0.171,  # C/m^2
    "ex4": -0.0406,
    "eSxx": 4.428*epsilon_0,
    "eSzz": 4.634*epsilon_0
}

Quartz_RH_1978 = {
    "rho": 2561e3,  # g/m^3
    "c11": 8.67362e10,  # Pa
    "c12": 6.98527e9,
    "c13": 1.19104e10,
    "c14": 1.79081e10,
    "c33": 1.07194e11,
    "c44": 5.79428e10,
    "ex1": 0.171,  # C/m^2
    "ex4": 0.0406,
    "eSxx": 4.428*epsilon_0,
    "eSzz": 4.634*epsilon_0
}

# data from B.A. Auld's book, Appendix 2
ZnO_auld = {  # Hex. 6mm
    "rho": 5680e3,   # g/m^3
    "c11": 209.7e9,  # Pa
    "c12": 121.1e9,
    "c13": 105.1e9,
    "c33": 210.9e9,
    "c44": 42.47e9,
    "ex5": -0.48,  # C/m^2
    "ez1": -0.573,
    "ez3": 1.32,
    "eSxx": 8.55*epsilon_0,
    "eSzz": 10.2*epsilon_0
}

# data form Comsol Multiphysics 5.4
ZnO_comsol = {  # Hex. 6mm
    "rho": 5680e3,   # g/m^3
    "c11": 209.714e9,  # Pa
    "c12": 121.14e9,
    "c13": 105.359e9,
    "c33": 211.194e9,
    "c44": 42.3729e9,
    "ex5": -0.480508,  # C/m^2
    "ez1": -0.567005,
    "ez3": 1.32044,
    "eSxx": 8.5446*epsilon_0,
    "eSzz": 10.204*epsilon_0
}

# data from Comsol Multipysics 5.4
AlN_comsol = {  # Hex. 6mm
    "rho": 3300e3,   # g/m^3
    "c11": 410e9,  # Pa
    "c12": 149e9,
    "c13": 99e9,
    "c33": 389e9,
    "c44": 125e9,
    "ex5": -0.48,  # C/m^2
    "ez1": -0.58,
    "ez3": 1.55,
    "eSxx": 9*epsilon_0,
    "eSzz": 9*epsilon_0
}

# data from Comsol Multiphysics 5.6
PVDF_comsol = {
    "rho": 1780e3,  # g/m^3
    "c11": 3.8e9,  # Pa
    "c12": 1.9e9,
    "c13": 0.9e9,
    "c22": 3.8e9,
    "c23": 0.9e9,
    "c33": 1.2e9,
    "c44": 0.7e9,
    "c55": 0.9e9,
    "c66": 0.9e9,
    "ez1": 0.024,  # C/m^2
    "ez2": 0.001,
    "ez3": -0.027,
    "eSxx": 7.4*epsilon_0,
    "eSyy": 9.3*epsilon_0,
    "eSzz": 7.6*epsilon_0
}


# === elastic materials ===
# data from B.A. Auld's book, Appendix 2
Al_auld = {  # Cubic m3m, crystal
    "rho": 2695e3,  # g/m^3
    "c11": 108e9,  # Pa
    "c12": 61.3e9,
    "c44": 28.5e9
}

# data from B.A. Auld's book, Appendix 2
Al_auld_poly = {  # Cubic m3m, crystal
    "rho": 2695e3,  # g/m^3
    "c11": 111e9,  # Pa
    "c44": 25e9
}

# data from B.A. Auld's book, Appendix 2
Au_auld = {  # Cubic m3m, crystal
    "rho": 19300e3,  # g/m^3
    "c11": 186e9,  # Pa
    "c12": 157e9,
    "c44": 42e9
}

# data from B.A. Auld's book, Appendix 2
Au_auld_poly = {  # Cubic m3m, crystal
    "rho": 19300e3,  # g/m^3
    "c11": 207e9,  # Pa
    "c44": 28.5e9
}

# data from B.A. Auld's book, Appendix 2
W_auld = {  # Cubic m3m, crystal
    "rho": 19200e3,  # g/m^3
    "c11": 502e9,  # Pa
    "c12": 199e9,
    "c44": 152e9
}

# data from B.A. Auld's book, Appendix 2
W_auld_poly = {  # Cubic m3m, crystal
    "rho": 19200e3,  # g/m^3
    "c11": 581e9,  # Pa
    "c44": 134e9
}
