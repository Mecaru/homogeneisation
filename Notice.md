# User Read-me

## Installation
Some sections of the code require you to have gfortran installed on your computer and available on your Path. Depending on your environnement, here are the steps you should follow when running the software for the fisrt time.

### Issue with gfortran and windows Path : 
Mac/windows : You may find a good gfortran compiler on this website http://mingw-w64.org/doku.php
Linux : apt install gfortran
If you have issue with adding gfortran with your Path on windows, this website may help https://davescience.wordpress.com/2012/03/21/gfortran-on-windows/. **Furthermore you should install mingw64 in C: problems occurs when your binaries are in a file with space in the name or too far from root.**

### Local with pip:
- Install JupyterLab
- Open a terminal on JupyterLab
- Exectute the following commands
    - pip3 install ipywidgets
    - pip3 install nodejs
    - pip3 install ipympl
    - jupyter labextension install @jupyter-widgets/jupyterlab-manager
    - jupyter labextension install jupyter-matplotlib
    - jupyter nbextension enable --py widgetsnbextension
- Restart JupyterLab

### Local with conda
- Install Anaconda
- Open a new terminal on JupyterLab
- Execute the following commands
    - conda install -c conda-forge ipywidgets
    - conda install -c conda-forge ipympl
    - conda install -c conda-forge nodejs
    - jupyter labextension install @jupyter-widgets/jupyterlab-manager
    - jupyter lab build
- Restart Anaconda and JupyterLab

### Compiling locally Fortran functions with f2py
- necessary tu use 'homogenization_main.ipynb'
- Open a new terminal on JupyterLab, and run :
    - f2py -c Fortran_tools.f -m fortran_tools
- This step must be performed only the first time you use the softaware, it has no effect after

### Running the software online
https://notebooks.ai is a website that allows you to run a JupyterLab notebook online freely.

---
## Using the software

The software enables you to compute the homogenized behaviors of simple microstructures using models implemented with Python. The three notebooks 'homogenization_main.ipynb','homogenization_main_II.ipynb', and 'homogenization_visco.ipynb' contain a user-friendly friendly interface. On the former, the user will find multiple sections:

- Computation of the homogenized behavior of defined microstructures

The section first invites the user to manually build inclusions by selecting their geometry and their behavior. Note that the inclusion is not generated as long as the 'Generate inclusion' button is not toggled. Once the inclusions generated, the user can create a microstructure by adding and removing any of the previously generated inclusions. Once again, the microstructure is generated after the 'Generate microstructure' button is toggled. The script will then display the compatible models and enable the user to chose the model he wants to use for the computation of the homogenized behavior.

The model comparison enables the user to compare multiple models results on the previously generated microstructure. 

- Automatic computation from a .txt file

Enables the computation of multiple homogenized behaviors from a text file. The input file must have a .txt extension and placed in the inputs/automatic_computation folder. example.txt gives an exemple for the expected input format. Its first line must be '*homogenisation'. Two different microstructure must be separated with a line containing the '*' sign.

When defining a microstructure, the fisrt line must contain the name of the model, the second contains the matrix behavior, and the other lines contain the inclusions info (one line per inclusion). Eah inclusion line contain the type of the inclusion (0:spheres, 1:ellipsoids), its behavior and its volume fraction. The script should be modified to also expect the ratio aspects of the ellispoids.

- Model description

Displays the description of the implemented models. The descriptions are Markdown files in the inputs/model_comparison folder. If the Latex equations are not properly displayed (known bug), the user should click on the blue bar on the left of the description twice.

- Estimates of parameters by an inverse method

Computes optimal microstructure parameters to reach a given target homogenized behavior. If the optimization problem admits multiple local minima, the program will only return one of them (not necessarily the best). Runs well for a few unknown parameters (volume fraction of an inclusion, or for a single unknown matrix behavior parameter), but can still be improved.

## Differences :

- 'homogenization_main.ipynb' and 'homogenization_main_II.ipynb'

The first needs a compilation of a fortran script with a numpy tool, f2py, whose use may be tricky on windows computer. It allows the import of a module containing fortran functions.
The second compiles first fortran functions, and then the rest of the code, avoiding the use of f2py, but resulting in a huge notebooks. 

- 'homogenization_main.ipynb' and 'homogenization_visco.ipynb'

The first handle only elastic cases, whereas the second can compute visco-elastic homogenised behavior for **single spherical isotropic inclusion**.