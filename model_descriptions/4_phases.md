# 4-phase model
[E.H.J. Maurer, 1990. An interlayer model to describe the physical properties of particulate composites. In: Ishida H, editor, Controlled interphases in composite materials. Elsevier Science Publishing Co, 491-504]

[E. Herve and A. Zaoui, 1993. n-Layered inclusion-based micromechanical modelling. Int. J. Engng Sci., 31, 1-10]

Notation:
$G$ shear modulus, 
$K$ Bulk modulus, 
$\nu$ Poisson's ratio, 
$E$ Young modulus, 
volume fraction of fillers $f$,
index $f$ for filler, index $m$ for matrix, index $i$ for interphase and $\bar{X}$ for the homogeneous medium.
## Spherical isotropic inclusions ramdomly dispersed

Account for an interphase of volume fraction and behavior
$f_i$, $K_i$, $G_i$, and $\nu_i$ 
with $f_m+f_f+f_i=1$. 
The 3-phase self-consistent model is recovered when $f_i=0$, among other special cases. The bulk modulus $\bar{K}$ of the composite is
\begin{equation}
\bar{K}=K_m+\frac{f_f+f_i}{\frac{1}{K_e-K_m}+\frac{3f_m}{3K_m+4G_m}}
\end{equation}
where
\begin{equation}
K_e=K_i+\frac{f_f}{\frac{f_f+f_i}{K_f-K_i}+\frac{3f_i}{3K_i+4G_i}}\,.
\end{equation}
The shear modulus $\bar{\mu}$ of the composite is the positive root of
\begin{equation*}
40\,\textrm{det}[X]\,\bar{G}^2+(2\,\textrm{det}[Y]+8\,\textrm{det}[Z])G_m\,\bar{G}-5\,\textrm{det}[T]\,G_m^2=0
\end{equation*}
where the determinants of four $10\times10$ matrices $[X]$, $[Y]$, $[Z]$, and $[T]$ are involved. The only non-zero elements of $[X]$ are
\begin{gather*}      
      X(1,1)=X(2,1)=G_m/G_i\,,\quad
      X(1,2)=-f_2\,G_m/(2G_i)\,,\\
      X(1,3)=X(2,3)=X(5,7)=X(6,7)=-1\,,\quad
      X(1,4)=4/f_5\,,\\
      X(1,5)=f_2/2\,,\quad
      X(1,6)=-a_i/f_p\,,\quad
      X(2,2)=b_p\,f_2\,G_m/G_i\,,\\
      X(2,4)=-8/(3f_5)\,,\quad
      X(2,5)=-b_i\,f_2\,,\quad
      X(2,6)=-c_i/f_p\,,\\
      X(3,1)=-X(3,3)=f_1\,,\quad
      X(3,2)=-X(3,5)=f_p\,,\quad
      X(3,4)=-1/f_4\,,\\
      X(3,6)=-1/f_2\,,\quad
      X(4,1)=-X(4,3)=f_1/2\,,\quad
      X(4,2)=d_p\,f_p\,,\\
      X(4,4)=1/(3f_4)\,,\quad
      X(4,5)=-d_i\,f_p\,,\quad
      X(4,6)=-e_i/f_2\,,\\
      X(5,3)=X(6,3)=G_i/G_g\,,\quad
      X(5,4)=-4G_i/(f'_5\,G_g)\,,\\
      X(5,5)=-f'_2\,G_i/(2G_g)\,,\quad
      X(5,6)=a_i\,G_i/(f'G_g)\,,\quad
      X(5,8)=4/f'_5\,,\\
      X(5,9)=f'_2/2\,,\quad
      X(5,10)=-a_g/f'\,,\quad
      X(6,4)=8G_i/(3f'_5\,G_g)\,,\\
      X(6,5)=b_i\,f'_2\,G_i/G_g\,,\quad
      X(6,6)=c_i\,G_i/(f'G_g)\,,\quad
      X(6,8)=-8/(3f'_5)\,,\\
      X(6,9)=-b_g\,f'_2\,,\quad
      X(6,10)=-c_g/f'\,,\quad
      X(7,3)=-X(7,7)=f'_1\,,\\
      X(7,4)=-X(7,8)=1/f'_4\,,\quad
      X(7,5)=-X(7,9)=f'\,,\\
      X(7,6)=-X(7,10)=1/f'_2\,,\quad
      X(8,3)=-X(8,7)=f'_1/2\,,\\
      X(8,4)=-X(8,8)=-1/(3f'_4)\,,\quad
      X(8,5)=d_i\,f'\,,\quad
      X(8,6)=e_i/f'_2\,,\\
      X(8,9)=-d_g\,f'\,,\quad
      X(8,10)=-e_g/f'_2\,,\quad
      X(9,7)=5/2\,,\\
      X(9,9)=1+3d_g\,,\quad
      X(9,10)=1+3e_g\,,\quad
      X(10,7)=1/2\,,\\
      X(10,8)=-1/3\,,\quad
      X(10,9)=d_g\,,\quad
      X(10,10)=e_g
\end{gather*}      
where
\begin{gather*}      
      f_1=f_p^{1/3},\quad
      f_2=f_p^{2/3},\quad
      f_4=f_p^{4/3},\quad
      f_5=f_p^{5/3},\\
      f'_1=(f')^{1/3},\quad
      f'_2=(f')^{2/3},\quad
      f'_4=(f')^{4/3},\quad
      f'_5=(f')^{5/3}
\end{gather*}      
with $f'=f_p+f_i$ and
\begin{gather*}      
      a_g=-2(5-\nu_m)/(5-4\nu_m),\quad
      a_i=-2(5-\nu_i)/(5-4\nu_i),\\
      b_g=(7+2\nu_m)/(6\nu_m),\quad
      b_p=(7+2\nu_f)/(6\nu_f),\\
      b_i=(7+2\nu_i)/(6\nu_i),\quad
      c_g=2(1+\nu_m)/(5-4\nu_m),\\
      c_i=2(1+\nu_i)/(5-4\nu_i),\quad
      d_g=(7-4\nu_m)/(12\nu_m),\\
      d_p=(7-4\nu_f)/(12\nu_f),\quad
      d_i=(7-4\nu_i)/(12\nu_i),\\
      e_g=(1-2\nu_m)/(5-4\nu_m),\quad
      e_i=(1-2\nu_i)/(5-4\nu_i)\,.
\end{gather*}      
The only elements of $[Y]$ that differ from $[X]$ are
\begin{gather*}      
      Y(9,7)=Y(10,8)=0\,,\quad
      Y(9,8)=-20/3\,,\\
      Y(9,9)=-1/2-b_g\,,\quad
      Y(9,10)=a_g-c_g\,,\quad
      Y(10,7)=5/2\,,\\
      Y(10,9)=1+3d_g\,,\quad
      Y(10,10)=1+3e_g
\end{gather*}      
whereas the only elements of $[Z]$ that differ from $[Y]$ are
\begin{gather*}      
      Z(9,7)=5/2\,,\quad
      Z(9,8)=0\,,\quad
      Z(9,9)=-1/2+3b_g/2\,,\\
      Z(9,10)=a_g+3c_g/2\,,\quad
      Z(10,7)=0\,,\quad
      Z(10,8)=5/3\,,\\
      Z(10,9)=1-2d_g\,,\quad
      Z(10,10)=1-2e_g
\end{gather*}      
and the only elements of $[T]$ that differ from $[Z]$ are
\begin{gather*}      
      T(10,7)=1\,,\quad
      T(10,8)=8/3\,,\quad
      T(10,9)=b_g\,,\quad
      T(10,10)=c_g\,.
\end{gather*}      

## Model validation
This model has been poorly validated on data from [Lipsinki P., Cherkaoui M., 2007. Four Phase Model: A New Formulation to Predict the Effective Elastic Moduli of Composites. Datas : Fig 3a p5

<img src="model_descriptions/model_validate/4_Phases_Lipsinki_G1.png" alt="drawing" width="400">
<img src="model_descriptions/model_validate/4_Phases_Lipsinki_K1.png" alt="drawing" width="400">