# Hashin and Shtrikman bounds
[Z., Hashin, S., Shtrikman, 1963. A Variational Approach to the Theory of
the Elastic Behaviour of Multiphase Materials. J. Mech. Phys. Solids, 11, 
127–140]

Notation:
$G$ shear modulus, 
$K$ Bulk modulus, 
$\nu$ Poisson's ratio, 
$E$ Young modulus, 
volume fraction of fillers $f$,
index $f$ for filler, index $m$ for matrix.
## Spherical isotropic inclusions ramdomly dispersed

The HS bounds are the tightest bounds without taking into consideration the geometry of the microstructure. 

Lower bound (Assuming the matrix is the softest phase):
$$K^-_{HS}=K_m+\frac{f}{\frac{1}{K_f-K_m}+\frac{3(1-f)}{3K_m+4G_m}}$$

$$G^-_{HS}=G_m+\frac{f}{\frac{1}{G_f-G_m}+\frac{6(1-f)(K_m+2G_m)}{5G_m(3K_m+4G_m})}$$

Upper bound:
$$K^+_{HS}=K_f+\frac{1-f}{\frac{1}{K_m-K_f}+\frac{3f}{3K_f+4G_f}}$$

$$G^+_{HS}=G_f+\frac{f}{\frac{1}{G_m-G_f}+\frac{6f(K_f+2G_f)}{5G_f(3K_f+4G_f})}$$


## Model validation
These models have been validated on data from [Hashin, Z., and Shtrikman, S., 1963. A Variational Approach to the Theory of the Elastic Behaviour of Multiphase Materials.J. Mech. Phys. Solids, 11, 127–140]. Fig 1 and 2 p136 : 

<img src="model_descriptions/model_validate/Hashin_Shtrickman_Ginf.png" alt="drawing" width="400">
<img src="model_descriptions/model_validate/Hashin_Shtrickman_Kinf.png" alt="drawing" width="400">
<img src="model_descriptions/model_validate/Hashin_Shtrickman_Gsup.png" alt="drawing" width="400">
<img src="model_descriptions/model_validate/Hashin_Shtrickman_Ksup.png" alt="drawing" width="400">

