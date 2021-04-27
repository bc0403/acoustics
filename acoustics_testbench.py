"""
This is a testbench for 'acoustics.py', mainly verified by B.A. Auld's book
Hao JIN, May 22, 2019
"""

from sympy import *
import numpy as np

from acoustics import ElasticMaterial, PiezoMaterial
from acoustics import Isotropic, Cubic, Hex6mm, Trig3m, Trig32
from acoustics import LN_comsol

init_printing(use_unicode=True)

# DEFINE SYMBOLS
rho, c, epsilon, e = symbols("rho c epsilon e")
alpha, beta, theta, gamma = symbols("alpha beta theta gamma")

c11, c12, c13, c14, c33, c44 = symbols('c11 c12 c13 c14 c33 c44')  # stiffness
ex1, ex4, ex5, ey2, ez1, ez3 = symbols(
    'ex1 ex4 ex5, ey2, ez1, ez3')  # piezoelectric stress constants
eSxx, eSzz = symbols('eSxx eSzz')  # dielectric constants


# WAVE PROPAGATION
"""
Here we definde the wave propagation matrices
"""
# wave propagates along [lx ly lz]
# lx^2 + ly^2 + lz^2 == 1
lx, ly, lz = symbols('lx ly lz')

li = Matrix([
    [lx, ly, lz]
])
# pprint(li)

lj = li.T
# pprint(lj)

liK = Matrix([
    [lx, 0, 0, 0, lz, ly],
    [0, ly, 0, lz, 0, lx],
    [0, 0, lz, ly, lx, 0]
])
# pprint(liK)

lLj = liK.T
# pprint(lLj)

trig3m = Trig3m(c11, c12, c13, c14, c33, c44, ex5, ey2, ez1, ez3, eSxx, eSzz)
trig32 = Trig32(c11, c12, c13, c14, c33, c44, ex1, ex4, eSxx, eSzz)
cubic = Cubic(c11, c12, c44, ex4, eSxx)
isotropic = Isotropic(c11, c44, eSxx)
hex6mm = Hex6mm(c11, c12, c13, c33, c44, ex5, ez1, ez3, eSxx, eSzz)

# DEMONSTRATIONS AND VERIFICATIONS
# === demonstrate the ElasticMaterial class ===
print("\n=== demonstrate the ElasticMaterial class ===")
my_elastic_material = ElasticMaterial(rho, c, epsilon)
my_elastic_material.get_description()

print("\nR matrix, rotation of theta by x-axis")
pprint(my_elastic_material.rotx_R(theta))

print("\nR matrix, rotation of theta by y-axis")
pprint(my_elastic_material.roty_R(theta))

print("\nR matrix, rotation of theta by z-axis")
pprint(my_elastic_material.rotz_R(theta))

# refer to page 74, Auld's book
Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz = \
    symbols("Rxx, Rxy, Rxz, Ryx, Ryy, Ryz, Rzx, Rzy, Rzz")

R = Matrix([
    [Rxx, Rxy, Rxz],
    [Ryx, Ryy, Ryz],
    [Rzx, Rzy, Rzz]
])

print("\nM matrix, according to rotational matrix R")
pprint(my_elastic_material.rot_M(R))

# https://en.wikipedia.org/wiki/Euler_angles
print("\nR matrix, rotation by Z-X-Z euler angles")
pprint(my_elastic_material.rot_euler_RM(alpha, beta, gamma)[0])

# === demonstrate the PiezoMaterial class ===
print("\n=== demonstrate the ElasticMaterial class ===")
my_piezo_material = PiezoMaterial(rho, c, epsilon, e)
my_piezo_material.get_description()

# === verify example 1.6 of auld's book, p.21-22 ===
print("\n=== example 1.6 ===")
C = symbols("C")
S = Matrix([
    [0, C, 0],
    [C, 0, 0],
    [0, 0, 0]
])
material1 = ElasticMaterial(rho, c, epsilon)
R = material1.rotz_R(theta)
pprint(R*S*R.T)

# === verify example 3.6 of Auld's book, p. 76-77 ===
print("\n=== example 3.6 ===")
cubic1 = Cubic(c11, c12, c44, ex4, eSxx)
material1 = ElasticMaterial(rho, cubic1.c, cubic1.eS)
print("Rotation of theta angle about the z-axis:")
R1 = material1.rotz_R(theta)
M1 = material1.rot_M(R1)
print("R matrix:")
pprint(R1)
print("M matrix:")
pprint(M1)
print("stiffness before rotation:")
pprint(material1.stiffness)
material1.rot_update(R1, M1)
print("stiffness after rotation:")
pprint(material1.stiffness)
pprint(trigsimp(material1.stiffness[0, 0]))
pprint(expand_trig(material1.stiffness[0, 0]))

# === verify example 3.6 of Auld's book, p. 77 ===
print("\n=== example 3.7 ===")
pprint(material1.rot_M(material1.roty_R(theta)))
print("\n")
pprint(material1.rot_M(material1.rotx_R(theta)))

# === verify example 3.8 of Auld's book, p. 82 ===
print("\n=== example 3.8 ===")
material1 = ElasticMaterial(rho, cubic1.c, cubic1.eS)
# using original axes
Gammaij1 = liK*material1.stiffness*lLj
Gammaij1_110 = Gammaij1.subs([(lx, 1), (ly, 1), (lz, 0)])
print("stiffness before rotation:")
pprint(material1.stiffness)
print("Gamma matrix using original axes:")
pprint(Gammaij1_110)
# using transformed axes
R1 = material1.rotz_R(pi/4)
M1 = material1.rot_M(R1)
material1.rot_update(R1, M1)
print("\nstiffness after rotation:")
pprint(material1.stiffness)

Gammaij1_100 = Gammaij1.subs([(lx, 1), (ly, 0), (lz, 0)])
print("Gamma matrix usring transformed axes:")
pprint(Gammaij1_100)

# === verify example 4.2 of Auld's book, p. 117 ===
print("\n=== example 4.2 ===")
hex6mm1 = Hex6mm(c11, c12, c13, c33, c44, ex5, ez1, ez3, eSxx, eSzz)
material1 = ElasticMaterial(rho, hex6mm1.c, hex6mm1.eS)
R1 = material1.roty_R(theta)
print("\nrotational matrix R:")
pprint(R1)
print("dielectric constants before rotation:")
pprint(material1.epsilon)
print("stiffness constants before rotation:")
pprint(material1.stiffness)

M1 = material1.rot_M(R1)
material1.rot_update(R1, M1)
print("dielectric constants after rotation:")
pprint(material1.epsilon)
print("stiffness constants after rotation:")
pprint(material1.stiffness)

Gammaij1 = liK*material1.stiffness*lLj
Gammaij1_010 = Gammaij1.subs([(lx, 0), (ly, 1), (lz, 0)])
print("Gamma matrix usring transformed axes:")
pprint(Gammaij1_010)


# === verify example 8.5 of Auld's book, p. 286 ===
print("\n=== example 8.5 ===")
material = PiezoMaterial(rho, cubic.c, cubic.eS, cubic.e)
# wave propagate along [1, 1, 0]
# using original axes
Gammaij = material.cal_Gamma(li, lj, liK, lLj)
Gammaij_110 = Gammaij.subs([(lx, 1/np.sqrt(2)), (ly, 1/np.sqrt(2)), (lz, 0)])
pprint(Gammaij_110)
# using rotated axes
material.rot_euler_update(pi/4, 0, 0)
Gammaij = material.cal_Gamma(li, lj, liK, lLj)
Gammaij_100 = Gammaij.subs([(lx, 1), (ly, 0), (lz, 0)])
pprint(Gammaij_100)


# === verify example 8.7 of Auld's book, p. 296 ===
print("\n=== example 8.7 ===")
material = PiezoMaterial(rho, hex6mm.c, hex6mm.eS, hex6mm.e)
material.get_description()
Gammaij = material.cal_Gamma(li, lj, liK, lLj)
Gammaij_100 = Gammaij.subs([(lx, 1), (ly, 0), (lz, 0)])
pprint(Gammaij_100)

# === verify example 8.9 of Auld's book, p. 300 ===
print("\n=== example 8.9 ===")
material = PiezoMaterial(rho, hex6mm.c, hex6mm.eS, hex6mm.e)
Gammaij = material.cal_Gamma(li, lj, liK, lLj)
Gammaij_x0z = Gammaij.subs([(ly, 0)])
pprint(Gammaij_x0z)


# === verify table 8.1 of Auld's book, p. 308 ===
print("\n=== table 8.1 ===")
material = PiezoMaterial(rho, cubic.c, cubic.eS, cubic.e)
print("\ncubic, along [110]:")
pprint(material.cal_Gamma(li, lj, liK, lLj).subs(
    [(lx, 1/np.sqrt(2)), (ly, 1/np.sqrt(2)), (lz, 0)]))

material = PiezoMaterial(rho, hex6mm.c, hex6mm.eS, hex6mm.e)
print("\nHex. 6mm, along [001]:")
pprint(material.cal_Gamma(li, lj, liK, lLj).subs(
    [(lx, 0), (ly, 0), (lz, 1)]))
print("\nHex. 6mm, along [100]:")
pprint(material.cal_Gamma(li, lj, liK, lLj).subs(
    [(lx, 1), (ly, 0), (lz, 0)]))

material = PiezoMaterial(rho, trig3m.c, trig3m.eS, trig3m.e)
print("\nTrig. 3m along [001]:")
print("\nGamma matrix using cD:")
pprint(material.cal_Gamma(li, lj, liK, lLj).subs(
    [(lx, 0), (ly, 0), (lz, 1)]))
print("\nGamma matrix using cE:")
pprint(super(PiezoMaterial, material).cal_Gamma(li, lj, liK, lLj).subs(
    [(lx, 0), (ly, 0), (lz, 1)]))

material = PiezoMaterial(rho, trig32.c, trig32.eS, trig32.e)
print("\nTrig. 32 along [001]:")
pprint(material.cal_Gamma(li, lj, liK, lLj).subs(
    [(lx, 0), (ly, 0), (lz, 1)]))


# wave propagates along [0, 0, 1], for isotropic material
material = ElasticMaterial(rho, isotropic.c, isotropic.eS)
Gamma = material.cal_Gamma(li, lj, liK, lLj).subs(
    [(lx, 0), (ly, 0), (lz, 1)]
)
print("\nisotropic along [001], Gamma matrix: ")
pprint(Gamma)