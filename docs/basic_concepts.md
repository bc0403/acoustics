[toc]

# Basic Concepts of Acoustics

## Particle Displacement, Strain, and Stress

Acoustics is the study of time-varying deformations in material media. Particle displacement, strain, stress are the basic concepts.

### Particle Displacement

<img src="basic_concepts.assets/image-20210423091241355.png" alt="image-20210423091241355" style="zoom:33%;" />

**Particle displacement field ($\mathbf u$ )** describe how particles are displaced from their equilibrium positions.

### Strain

**Strain field ($\mathbf S$)** can be derived from particle displacement as
$$
S_{ij}(\mathbf r, t)=\frac{1}{2}(\frac{\partial u_i}{\partial r_j} + \frac{\partial u_j}{\partial r_i})
$$
where $\mathbf r$ is the position vector, $\mathbf r = \hat xx + \hat yy + \hat zz$

## Wave Propagating and Polarization

Polarization is the time-varying direction of $\mathbf u$. For example, when the acoustic wave propagates along $y$ direction, there are three basic polarizations (two pure **shear** and one pure **longitudinal**):

<img src="basic_concepts.assets/image-20210423093901548.png" alt="image-20210423093901548" style="zoom:33%;" />

## Transformation



## Symbolic Notation

* The shape of $\mathbf u$ is (3, 1);

* The shape of $\mathbf S$  is (3, 3)
  $$
  \mathbf S = \begin{bmatrix}
  S_{xx} & S_{xy} & S_{xz} \\
  S_{xy} & S_{yy} & S_{yz} \\
  S_{xz} & S_{yz} & S_{zz} \\ 
  \end{bmatrix}
  = \begin{bmatrix}
  S_1 & \frac{1}{2}S_6 & \frac{1}{2}S_5 \\
  \frac{1}{2}S_6 & S_2 & \frac{1}{2}S_4 \\
  \frac{1}{2}S_5 & \frac{1}{2}S_4 & S_3 \\
  \end{bmatrix}
  $$
  and (6, 1) in **abbreviated system**
  $$
  \mathbf S = \begin{bmatrix}
  S_1 \\ S_2 \\ S_3 \\ S_4 \\ S_5 \\ S_6 \\
  \end{bmatrix}
  $$

* The shape of $\mathbf T$ is same as $\mathbf S$

* Therefore, 
  $$
  S_I =\nabla_{Ij}u_j
  $$
  where
  $$
  \nabla_{Ij} = \begin{bmatrix}
  \frac{\partial}{\partial x} & 0 & 0 \\
  0 & \frac{\partial}{\partial y} & 0 \\
  0 & 0 & \frac{\partial}{\partial z} \\
  0 & \frac{\partial}{\partial z} & \frac{\partial}{\partial y} \\
  \frac{\partial}{\partial z} & 0 & \frac{\partial}{\partial x} \\
  \frac{\partial}{\partial y} & \frac{\partial}{\partial x} & 0 \\
  \end{bmatrix}
  $$
  

## References

* Auld, B. A. (1973). *Acoustic fields and waves in solids*. John Wiley & Sons.