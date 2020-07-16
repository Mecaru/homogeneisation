# Differential scheme

[R.M., Christensen, 1990. A critical evaluation for a class of micromechanics models. J. Mech. Phys. Solids, 38, 379-404]

Notation:
$G$ shear modulus, 
$K$ Bulk modulus, 
$\nu$ Poisson's ratio, 
$E$ Young modulus, $C$ elastic stiffness tensor,
volume fraction of fillers $f$,
index $f$ for filler, index $m$ for matrix, and $\bar{X}$ for the homogeneous medium. 

## Spherical isotropic inclusions ramdomly dispersed

The idea is to add progressively an infinitesimal amount of the reinforcing phase (the inclusion) until reaching the desired volume fraction $f$.

The differential equation to solve are:

$$\frac{d\bar{G}}{df}+\frac{15(1-\bar{\nu})(\bar{G}-G_f)}{(1-f)(7-5\bar{\nu}+2(4-5\bar{\nu})\frac{G_f}{\bar{G}})}=0$$

$$\frac{d\bar{K}}{df}+\frac{(\bar{K}-K_f)}{(1-f)(1+\frac{K_f-\bar{K}}{\bar{K}+\frac{4}{3}\bar{G}})}=0$$

with $\bar{\nu}=\frac{3\bar{K}-2\bar{G}}{2(3\bar{K}+\bar{G})}$ and limit conditions: 

for $c=0$ $\bar{G}=G_m$ and $\bar{K}=K_m$

for $c=1$ $\bar{G}=G_f$ and $\bar{K}=K_f$

## Ellipsoidal inclusions ramdomly dispersed
We consider the general case of anisotropic inclusions randomly dispersed in an isotropic matrix. The resulting equivalent homogeneous material is isotropic. 
The Eshelby tensor $S$ is given in the book [T., Mura, 1987. Micromechanics of defects in solids, second revised edition, Dordecht: Kluwer].

$n$ inclusions are randomly dispersed in space, and each inclusion $i$ is considered as a filler constitutive phase of volume fraction $f_i=f/n$. 

We start with $f=0$ and $\bar{C}_0=C_m$ and add $\delta f$ volume fraction of fillers at every step (with $\delta f$ small enough).
At each step $n+1$, the following constitutive equations are computed:
$$A_{n_i}^D=Id + S_{n}^{Esh} \bar{C}_{n}^{-1} (C_{f_i}-\bar{C}_n)$$
with $S_{n}^{Esh}$ depends on the geometry of the inclusion and of the behavior $\bar{C}_n$, 
$$\bar{C}_{n+1}=\bar{C}_n+(\delta f/n) \sum_{i}^n (C_{f_i}-\bar{C}_n):A_{n_i}^D$$ 
with $C_{f_i}$ and $A_{n_i}^{D}$ the behavior of the filler and of the self-consistent localization tensor computed in the material reference. 

Note that the Eshelby tensor has been calculated with elliptic integrals that may diverge for anisotropic matrices.
## Model validation
The model has been validated on data from : [Christensen,  R.M.,  1990.  A critical evaluation for a class of micromechanics models]. Equation (7) and (8) p5

<img src="model_descriptions/model_validate/Differentiel_Christensen_G.png" alt="drawing" width="600">
<img src="model_descriptions/model_validate/Differentiel_Christensen_K.png" alt="drawing" width="600">

and : [Jithender J.Timothy,Günther Meschke,2015. A cascade continuum micromechanics model for the effective elastic properties of porous materials]. Fig A15 p9

<img src="model_descriptions/model_validate/Differential_Timothy_K1.png" alt="drawing" width="600">