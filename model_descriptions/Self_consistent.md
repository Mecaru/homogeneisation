# Self-consistent model 
[R., Hill, 1965. A Self-Consistent Mechanics of Composite Materials. J.
Mech. Phys. Solids, 13, 213–222]

Notation:
$G$ shear modulus, 
$K$ Bulk modulus, 
$\nu$ Poisson's ratio, 
$E$ Young modulus, $C$ elastic stiffness tensor,
volume fraction of fillers $f$,
index $f$ for filler, index $m$ for matrix, and $\bar{X}$ for the homogeneous medium.
## General equations for a matrix filled with identical inclusions

$$\bar{C}=C_m+f(C_f-C_m):A^{SC}$$

with $C_m$ the matrix behavior $C_f$ the filler behavior, $f$ its volume fraction, and finally $A^{SC}=\left( 1+S \bar{C}\,^{-1}(C_f-\bar{C})\right)^{-1}$
with the Eshelby tensor $S(\bar{C})$ that depends on the inclusion geometry and the homogeneous material behaviour.

## Spherical isotropic inclusions ramdomly dispersed
Interactions between inclusions are taken into account assuming that one inclusion is placed in the homogeneous medium

$$\bar{K}=K_m+\frac{f(K_f-K_m)(3\bar{K}+4\bar{G})}{3K_f+4\bar{G}}$$
$$\bar{G}=G_m+\frac{5f\bar{G}(G_f-G_m)(3\bar{K}+4\bar{G})}{3\bar{K}(3\bar{G}+2G_f)+4\bar{G}(3G_f+2\bar{G})}$$

These equations are solved with a mere fixed-point method. Initialization is $\bar{G}=G_m$ and $\bar{K}=K_m$ for $f=0$ and 
$\bar{G}=\bar{G}(f_n)$ and $\bar{K}=\bar{K}(f_n)$ for $f_{n+1}=(n+1)\delta f$.

Note that this model may diverge for conditions such as incompressible matrix and compressible inclusion at moderate to high filler concentrations depending on the constitutive phase behaviors.

## Ellipsoidal inclusions ramdomly dispersed
We consider the general case of anisotropic inclusions randomly dispersed in an isotropic matrix. The resulting equivalent homogeneous material is isotropic. 
The Eshelby tensor $S$ is given in the book [T., Mura, 1987. Micromechanics of defects in solids, second revised edition, Dordecht: Kluwer].

$n$ inclusions are randomly dispersed in space, and each inclusion $i$ is considered as a filler constitutive phase of volume fraction $f_i=f/n$. The consitutive equations writes now as:
$$\bar{C}=C_m+\sum_{i}^n f_i(C_{f_i}-C_m):A_i^{SC}$$ 
with $C_{f_i}$ and $A_i^{SC}$ the behavior of the filler and of the self-consistent localization tensor computed in the material reference. 

Note that the Eshelby tensor has been calculated with elliptic integrals that may diverge for anisotropic matrices.
## Model validation
The model has been validated on data from [V. Marcadon , E. Herve, A. Zaoui,  2007. Micromechanical modeling of packing and size effects in particulate composites]. Fig 3 p7
<img src="model_descriptions/model_validate/SC_Marcadon_G.png" alt="drawing" width="400">