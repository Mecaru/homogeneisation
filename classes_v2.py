# -*- coding: utf-8 -*-
"""
Homogeneisation - classes.py

Définition des classes utilisées. Version 2 de classes, avec une définition simplifiée des modèles et prise en compte de la visco-élasticité.

Les items à modifier/améliorer/ajouter sont marqués en commentaires précédés de la mention "TODO".

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

11/06/2020
"""

#%% Useful packages
import numpy as np
from scipy import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#%% Microstructure classes
class Inclusion:
    """
    Contains info of the inclusion (type, behavior, shape).
    """
    
    def __init__(self, type_inclusion, behavior, aspect_ratio=[1.,1.], name=None, frequency=[], abscissa="frequency"):
        """
    
        type_inclusion: (int), 0 for spherical inclusion, 1 for ellipsoids
        aspect_ratio: (tuple), tuple of two floats representing the length ratio of axis 2 and 3 of the ellipsoid to the length of the axis 1 
        behavior: (dict), contains values of matrix behavior : E,K,G and nu in the isotropic case, compliance and stiffness matrix in the anisotropic
        frequency: (list), list of frequencies/temperatures associated with visco-elastic parameters
        abscissa: (str),  "frequency" or "temperature", it indicates the physical nature of the values in the frequency list
        inc_and_int: ([InclusionAndInterphase,int]), None by default, linked to the InclusionAndInterphase class instance to which the inclusion belongs if it exists, and an integer, 0 for the inclusion, 1 for the interphase
        """
        self.type_inclusion = type_inclusion
        self.aspect_ratio = aspect_ratio
        self.behavior = complete_behavior(behavior)
        self.name = name
        self.frequency = frequency
        self.abscissa = abscissa
        self.inc_and_int = None
    
    def type_to_str(self):
        """
        Transforms an integer "type_inclusion" into the corresponding string (example: 0 --> "spheres") 
        TODO: synchronize this function with the main using a dictionary to make it easier to add inclusion types
        """
        type_inclusion = self.type_inclusion
        try:
            result = dict_types[type_inclusion]
        except KeyError:
            # Le type spécifié n'est pas répertorié dans le dictionnaire
            result = None
        return result
    
    def __str__(self):
        """
        Presentation of the instance.
        """
        str_type_inclusion = self.type_to_str()
        string = "{}, {}".format(self.name, str_type_inclusion)
        if self.type_inclusion != 0:
            string += " (ratios={})".format(self.aspect_ratio)
        for parameter, value in self.behavior.items():
            if type(value) not in [list, np.ndarray]:
                string += ", {}: {:.2f}".format(parameter, value)
            else:
                string += ", {}: list".format(parameter)
        return string

    def __repr__(self):
        return str(self)

    def change_parameter(self, parameter, new_value):
        """
        Changes the value of the parameter if it exists. Updates the behavior with the function "complete_behavior".
        """
        try:
            self.behavior[parameter] = new_value
            self.behavior = complete_behavior(self.behavior)
        except:
            None

    def graph_parameter(self):
        """
       Plots the graph of the evolution of the visco-elastic parameters, if they exist
        """
        if self.frequency == []:
            None # L'inclusion ne contient pas de paramètres visco-élastiques
        else:
            plt.figure()
            # for parameter, values in self.behavior.items():
            for parameter in ['K', 'G']:
                values = self.behavior[parameter]
                # Le paramètre est visco-élastique
                if self.abscissa == "temperature":
                    plt.semilogy(self.frequency, values.real, '.', label=parameter+"'")
                    plt.semilogy(self.frequency, values.imag, '.', label=parameter+"''")
                elif self.abscissa == "frequency":
                    plt.loglog(self.frequency, values.real, '.', label=parameter+"'")
                    plt.loglog(self.frequency, values.imag, '.', label=parameter+"''")
                plt.legend()
                plt.xlabel(self.abscissa)
                plt.ylabel("Parameter value")
                plt.title("Inclusion visco-elastic behavior")
                plt.xlim(min(self.frequency), max(self.frequency))
            plt.show()

class InclusionAndInterphase:
    """
    Instance representing an inclusion and the associated interphase.
    """            
    
    def __init__(self, inclusion, interphase, name=None):
        """
        Inclusion: Inclusion Class Instance
        interphase: Inclusion class instance, represents the interphase associated with inclusion.
        The interphase and the inclusion must be of the same type (spheres or ellipsoids with the same aspect ratios).
        """
        assert inclusion.aspect_ratio==interphase.aspect_ratio
        self.inclusion = inclusion
        self.interphase = interphase
        self.name = name
        self.aspect_ratio = inclusion.aspect_ratio
        # Modification de l'attribut inc_and_int des inclusion et interphase
        inclusion.inc_and_int = [self, 0]
        interphase.inc_and_int = [self, 1]
        
    def __str__(self):
        string = "Inclusion + Interphase\n"
        string += "Inclusion: {}\n".format(str(self.inclusion))
        string += "Interphase: " + str(self.interphase)
        return string
    
    def __repr__(self):
        return str(self)

class Microstructure:
    """
    Contains information on the microstructure (matrix behaviour and inclusions)
    Contains a function that returns the Hashin-Shtrickman bounds for the microstructure in question. 
    TODO: Generalize these terminals to n phases (and not 2 as is the case here).
    """
    
    def __init__(self, behavior, dict_inclusions=dict(), frequency=[], abscissa="frequency"):
        """
       list_inclusions: (dictionnary such as {inclusion: f_i} with inclusion as an instance of class Inclusion and f_i as the volume fraction of this type of inclusion. inclusion can also be an instance of class InclusionAndInterphase, in this case, f_i is a tuple of floaters
        behavior: (dict), contains the values of the parameters of the behavior matrix
        frequency: list of frequencies associated with viscoelastic parameters
        """
        self.dict_inclusions = dict_inclusions
        self.behavior = complete_behavior(behavior)
        # Calcul de la fraction volumique de matrice f_m
        self.f_matrix = self.compute_fm()
        self.frequency = frequency
        self.abscissa = abscissa
        
    def __str__(self):
        string = "Microstructure\nf_m = {:.2f}, matrix".format(self.f_matrix, self.behavior)
        for parameter, value in self.behavior.items():
            if type(value) not in [list, np.ndarray]:
                string += ", {}: {:.2f}".format(parameter, value)
            else:
                string += ", {}: list".format(parameter)
        dict_inclusions = self.dict_inclusions
        # Display of all inclusion in microstructure
        for inclusion in dict_inclusions.keys():
            fi = dict_inclusions[inclusion]
            string += "\nf_i = {}, ".format(fi) + str(inclusion)
        return string

    def compute_fm(self):
        """
        1/ Checks if the given list of inclusions is consistent (i.e. the sum of  volume fractions inclusions is less  than 1). Else, generates an error.
        2/ If no error is generated, calculates the matrix volume fraction.
        """
        total_fi = 0 # Total of volumic fractions of inclusions
        dict_inclusions = self.dict_inclusions
        for inclusion in dict_inclusions.keys():
            fi = dict_inclusions[inclusion]
            # Case  inclusions + interphase
            if type(fi)==list:
                total_fi += fi[0] + fi[1]
            # Case simple inclusions
            else:
                total_fi += fi
        if total_fi >= 1:
            raise NameError("The total volumic fractions of the inclusions exceed 1")
        else :
            f_m = 1 - total_fi
            return f_m

    def change_fi(self, inclusion, new_f):
        """
        Updates the volume fraction of the inclusion or adds it to the dictionary if it was not present.
        Updates the volume fraction of the matrix.
        """
        self.dict_inclusions[inclusion] = new_f
        self.f_matrix = self.compute_fm()
    
    def change_parameter(self, parameter, new_value):
        """
        Change the value of the parameter if it exists. Updates the behavior with the function "complete_behavior".
        """
        try:
            self.behavior[parameter] = new_value
            self.behavior = complete_behavior(self.behavior)
        except:
            None

    def graph_parameter(self):
        """
        Plots the graph of the evolution of the visco-elastic parameters if they exist.
        """
        if self.frequency == []:
            None # Inclusion is not visco-elastic
        else:
            plt.figure()
            # for parameter, values in self.behavior.items():
            for parameter in ['K', 'G']:
                values = self.behavior[parameter]
                # Le paramètre est visco-élastique
                if self.abscissa == "temperature":
                    plt.semilogy(self.frequency, values.real, '.', label=parameter+"'")
                    plt.semilogy(self.frequency, values.imag, '.', label=parameter+"''")
                elif self.abscissa == "frequency":
                    plt.loglog(self.frequency, values.real, '.', label=parameter+"'")
                    plt.loglog(self.frequency, values.imag, '.', label=parameter+"''")
                plt.legend()
                plt.xlabel(self.abscissa)
                plt.ylabel("Parameter value")
                plt.title("Matrix visco-elastic behavior")
                plt.xlim(min(self.frequency), max(self.frequency))
            plt.show()

    def draw(self):
        """
         Method for drawing the microstructure.
        """
        inclusions = list(self.dict_inclusions.keys())
        n_fig = len(inclusions)
        if n_fig==0:
            # Microstructure sans inclusion
            return None
        fig = plt.figure(figsize=(n_fig*3 ,3))
        for index, instance in enumerate(inclusions):
            fi = self.dict_inclusions[instance]
            if type(instance)==Inclusion:
                inclusion = instance
                f_inc = fi
                interphase = None
            else:
                inclusion = instance.inclusion
                interphase = instance.interphase
                f_inc = fi[0]
                f_int = fi[1]
                
            ### Draw inclusion
            ax = fig.add_subplot(1, n_fig, index+1, projection='3d')
            # compute radius for a 10X10X10-sized VER
            c1, c2 = inclusion.aspect_ratio
            a = (1000*f_inc/(4/3*pi*c1*c2))**(1/3)
            b = c1*a
            c = c2*a

            # Radii:
            rx, ry, rz = np.array([a, b, c])

            # Set of all spherical angles:
            u = np.linspace(0, 2 * np.pi, 100)
            v = np.linspace(0, np.pi, 100)

            # Cartesian coordinates that correspond to the spherical angles:
            # (this is the equation of an ellipsoid):
            x = rx * np.outer(np.cos(u), np.sin(v))
            y = ry * np.outer(np.sin(u), np.sin(v))
            z = rz * np.outer(np.ones_like(u), np.cos(v))

            # Plot:
            ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')
            
            ### Draw interphase
            if interphase!=None:
                a = (1000*(f_inc+f_int)/(4/3*pi*c1*c2))**(1/3)
                b = c1*a
                c = c2*a
    
                # Radii:
                rx, ry, rz = np.array([a, b, c])
    
                # Set of all spherical angles:
                u = np.linspace(0, np.pi, 100)
                v = np.linspace(0, np.pi, 100)
    
                # Cartesian coordinates that correspond to the spherical angles:
                # (this is the equation of an ellipsoid):
                x = rx * np.outer(np.cos(u), np.sin(v))
                y = ry * np.outer(np.sin(u), np.sin(v))
                z = rz * np.outer(np.ones_like(u), np.cos(v))
    
                # Plot:
                ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='r')
            
            ### Draw edges of VER
            # Adjustment of the axes, so that they all have the same span:
            max_radius = 5
            for axis in 'xyz':
                getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))

            # Cube 
            points = 5*np.array([[-1, -1, -1],
                                  [1, -1, -1 ],
                                  [1, 1, -1],
                                  [-1, 1, -1],
                                  [-1, -1, 1],
                                  [1, -1, 1 ],
                                  [1, 1, 1],
                                  [-1, 1, 1]])

            r = [-5,5]
            X, Y = np.meshgrid(r, r)
            one = 5*np.ones(4).reshape(2, 2)
            ax.plot_wireframe(X,Y,one, alpha=0.5)
            ax.plot_wireframe(X,Y,-one, alpha=0.5)
            ax.plot_wireframe(X,-one,Y, alpha=0.5)
            ax.plot_wireframe(X,one,Y, alpha=0.5)
            ax.plot_wireframe(one,X,Y, alpha=0.5)
            ax.plot_wireframe(-one,X,Y, alpha=0.5)
            ax.scatter3D(points[:, 0], points[:, 1], points[:, 2])
        plt.show()
            
     ## Compute HASHIN-SHTRICKMAN bounds ##########  
    
    def khs(k1, g1, c1, k2, g2, c2):
        numerator = c2*(k2-k1)
        denominator = 1+3*c1*(k2-k1)/(4*g1+3*k1)
        return k1+numerator/denominator
    
    def ghs(k1, g1, c1, k2, g2, c2):
        numerator = c2*(g2-g1)
        denominator = 1+6*c1*(g2-g1)*(k1+2*g1)/((3*k1+4*g1)*5*g1)
        return g1+numerator/denominator
        
    def Hashin_bounds(self):
        """
        Gives the Hashin-Shtrikman bounds for single phase, isotropic...
        TODO: add the case of multiple inclusions
        """
        fm = self.f_matrix
        f = 1-fm
        km,gm = self.behavior["K"], self.behavior["G"]
        
        for inclusion in self.dict_inclusions.keys():
            try:
                kf,gf=inclusion.behavior["K"],inclusion.behavior["G"]
            except:
                return None
        
        ksup=max(Microstructure.khs(km,gm,fm,kf,gf,f),Microstructure.khs(kf,gf,f,km,gm,fm))
        kinf=min(Microstructure.khs(km,gm,fm,kf,gf,f),Microstructure.khs(kf,gf,f,km,gm,fm))
        gsup=max(Microstructure.ghs(km,gm,fm,kf,gf,f),Microstructure.ghs(kf,gf,f,km,gm,fm))
        ginf=min(Microstructure.ghs(km,gm,fm,kf,gf,f),Microstructure.ghs(kf,gf,f,km,gm,fm))
            
        
        return { 'Ginf': ginf, 'Gsup': gsup, 'Kinf': kinf, 'Ksup': ksup }
    
    def Voigt_Bound(self) : 
        
        fm = self.f_matrix
        f = 1-fm
        km,gm = self.behavior["K"], self.behavior["G"]
        
        for inclusion in self.dict_inclusions.keys():
            try:
                kf,gf=inclusion.behavior["K"],inclusion.behavior["G"]
            except:
                return None
        K_voigt = Km*fm + Kf*f
        G_voigt = Gm*gm + Gf*f
        
        return complete_behavior({'G':G_voigt , 'K':K_voigt})
    
    def Reuss_Bound(self) : 
        
        fm = self.f_matrix
        f = 1-fm
        km,gm = self.behavior["K"], self.behavior["G"]
        
        for inclusion in self.dict_inclusions.keys():
            try:
                kf,gf=inclusion.behavior["K"],inclusion.behavior["G"]
            except:
                return None
        K_reuss = 1/(fm/Km + f/Kf )
        G_reuss = 1/(fm/Gm + f/Gf )
        
        return complete_behavior({'G':G_reuss , 'K':K_reuss})
        

#%% Models classes
class Model:
    """
   Generic parent class of all model classes. 
    Contains the method for verifying the model's assumptions about a microstructure, as well as the method called when calculating the homogenized behavior.
    """
    
    def __str__(self):
        """
         Textual description of the model.
        """
        return self.name + " model"
    
    def __repr__(self):
        """
         Textual description of the model.
        """
        return str(self)
    
    def check_hypothesis(self, microstructure):
        """
         Checks if the microstructure verifies the hypothesis of the model, returns a boolean.
        """
        # Recovery of microstructure inclusions
        if self.behavior_condition=='isotropic':
            behavior_condition = set(['K', 'G', 'E', 'nu'])
        elif self.behavior_condition=='anisotropic':
            behavior_condition = set(['C', 'S'])
        # Récupération des inclusions de la microstructure
        dict_inclusions = microstructure.dict_inclusions
        instances = list(dict_inclusions.keys())
        n_instances = len(instances)
        # Initialisation of result
        result = True
        # Checking the number of inclusions
        if n_instances != self.n_inclusions:
             result = False
        # Checking the presence or absence of an interphase
        for instance in instances:
            if (type(instance)==InclusionAndInterphase)!=self.interphase:
                result = False
        # Creation of a list of inclusions without interphase
        inclusions = []
        for instance in instances:
            if type(instance)==InclusionAndInterphase:
                inclusions += [instance.inclusion, instance.interphase]
            else:
                inclusions.append(instance)
        # Checking the type of inclusion
        for inclusion in inclusions:
            if inclusion.type_inclusion != self.type_inclusion:
                result = False
        # Checking the behaviour of inclusions and matrix
        for element in inclusions + [microstructure]:
            if not set(element.behavior.keys()).issubset(behavior_condition):
                result = False
        # Returns result
        return result
    
    def compute_h_behavior(self, microstructure):
        """
        Computes the homogenized behavior of the microstructure with the chosen model.
        Verifies that the model fits the input microstructure.
        If the input elements are not viscoelastic (i.e. if the frequency list of the microstructure is empty), returns a dictionary of real parameters in the form {"parameter": value(float)}. In the isotropic case, also calculates the missing parameters (mu and E, or K and G).
        Otherwise, performs a loop on the frequency values, then returns a complete behaviour dictionary (with missing parameter values) of the form {"parameter": values(complex)]}.
        """
        # Verification of the conditions of application
        compatible = self.check_hypothesis(microstructure)
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        frequency = microstructure.frequency
        # Creation of the dictionary containing only Inclusions
        inclusions = {}
        for instance, f in microstructure.dict_inclusions.items():
            if type(instance)==InclusionAndInterphase:
                inclusions[instance.inclusion] = f[0]
                inclusions[instance.interphase] = f[1]
            else:
                inclusions[instance] = f
        # Cas élastique
        if not list(frequency):
            Cm = microstructure.behavior
            # Retrieving inclusion behavior in the format : [(inclusion.behavior, f, aspect_ratio)]
            inclusion_behaviors = [(inclusion.behavior, f, inclusion.aspect_ratio) for (inclusion,f) in inclusions.items()]
            # Calculation of homogenized behaviour
            h_behavior = self.compute_behavior(Cm, inclusion_behaviors)
            h_behavior = {parameter: value.real for (parameter,value) in h_behavior.items()} # Conversion of possibly complex values into real values
            h_behavior = complete_behavior(h_behavior)
        # Visco-elastic case
        else:
            # Initialisation of result
            h_behavior = {}
            # Calculation of the behavior as a function of frequency
            for i in range(len(frequency)):
                # Recovery of matrix behavior at frequency i
                Cm = {parameter: values[i] for (parameter,values) in microstructure.behavior.items()}
                # Retrieving inclusion behavior at frequency i
                inclusion_behaviors = [] # Initialisation
                for inclusion, f in inclusions.items():
                    inclusion_behavior = {parameter: values[i] for (parameter, values) in inclusion.behavior.items()}
                    inclusion_behaviors.append((inclusion_behavior, f, inclusion.aspect_ratio))
                # Calculation of the homogenized behavior at frequency i
                h_behavior_i = self.compute_behavior(Cm, inclusion_behaviors)
                h_behavior_i = complete_behavior(h_behavior_i)
                # Adds it to the list of behaviors
                for parameter, value in h_behavior_i.items():
                    try:
                        h_behavior[parameter].append(value)
                    except KeyError:
                        # Creation of the input associated with the parameter
                        h_behavior[parameter] = [value]
        # Return of the result
        return h_behavior
                    
        
class Mori_Tanaka(Model):
    """
   Mori-Tanaka model. Contains:
    - A description function of the model
    - A function that returns the homogenized behavior of the microstructure.
    """
    
    def __init__(self):
        """
        Definition of model hypotheses.
        """
        self.type_inclusion = 0 # Spheres
        self.behavior_condition = 'isotropic'  # The model is applied to microstructures whose inclusions and matrix are isotropic.
        self.n_inclusions = 1 # Number of different types of inclusions
        self.interphase = False # True if the model works on inclusions with interphase
        self.name = "Mori-Tanaka"
    
    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Calculates the equivalent homogeneous elastic behaviour. 
        Returns a dictionnary of behavior.
        Cm: (dict), dictionary of matrix behavior
        inclusion_behaviors(list), format [(Cf, f, aspect_ratio)] with Cf the dictionaries of inclusion behaviour, and aspect_ratio (a tuple with the two shape ratio values)
        """
        # Retrieving matrix behavior
        Km = Cm['K']
        Gm = Cm['G']
        # Retrieving inclusion behavior
        Cf, f, ratio = inclusion_behaviors[0]
        Kf = Cf['K']
        Gf = Cf['G']
        # Computation of Gh
        denominator = 5*Gm*(3*Km+4*Gm)+6*(1-f)*(Gf-Gm)*(Km+2*Gm)
        numerator = 5*f*Gm*(Gf-Gm)*(3*Km+4*Gm)
        Gh = Gm + numerator/denominator
        # Computation of Kh
        denominator = 3*Km+4*Gm+3*(1-f)*(Kf-Km)
        numerator = f*(Kf-Km)*(3*Km+4*Gm)
        Kh = Km + numerator/denominator
        return {'K': Kh, 'G': Gh}    



class Differential_Scheme(Model):
    """
    Differential scheme
    """
    
    def __init__(self):
        """
        Definition of model hypotheses.
        """
        self.type_inclusion = 0
        self.behavior_condition = 'isotropic' 
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes  
        self.interphase = False
        self.name = "Differential"
    
    ## Useful functions to compute homogenized behavior
    
    def deriv(module, f):
        """
        Function that computes the derivatives of the parameters K and G in relation to the volume fraction of inclusion. Designed to be called by the odeint function during numerical integration.
        module: list, contains the real and imaginary values of the current K, G parameters as well as Kf and Gf specific to the inclusion.
        f: float, current inclusion volume fraction.
        """
        K1, K2, G1, G2, Kf1, Kf2, Gf1, Gf2 = module
        # Creation of complex parameters
        K = K1 + K2*1j
        G = G1 + G2*1j
        Kf = Kf1 + Kf2*1j
        Gf = Gf1 + Gf2*1j
        nu = (3*K-2*G)/(6*K+2*G)
        # Computation of dK
        numerator = K-Kf
        denominator = (1-f)*(1+(Kf-K)/(K+4*G/3))
        dK = -numerator/denominator
        dK1, dK2 = dK.real, dK.imag
        # Computation of dG
        numerator = 15*(1-nu)*(G-Gf)
        denominator = (1-f)*(7-5*nu+2*(4-5*nu)*Gf/G)
        dG = -numerator/denominator
        dG1, dG2 = dG.real, dG.imag
        
        return np.array([dK1, dK2 ,dG1, dG2] + 4*[0])
    
    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Computes the equivalent homogenized behavior of the microstructure. Returns a dict with the calculated parameters.
        """
        # Retrieving matrix behavior
        Km = Cm['K']
        Gm = Cm['G']
        # Retrieving inclusion behavior
        Cf, f_finale, ratio = inclusion_behaviors[0]
        Kf = Cf['K']
        Gf = Cf['G']
        # Initialisation of integration parameters
        npoints = 100 # Number of integration points
        f = np.linspace(0, f_finale, npoints) # List of dilute volumic fractions
        initial_module = []
        for parameter in [Km, Gm, Kf, Gf]:
            initial_module += [parameter.real, parameter.imag]
        initial_module = np.array(initial_module)
        # Integration of differential equation
        module = odeint(Differential_Scheme.deriv, initial_module, f)
        # Retrieving final result
        final_module = module[-1]
        Kh1, Kh2, Gh1, Gh2 = final_module[:4]  
        # Return result
        return {'K': Kh1+1j*Kh2, 'G': Gh1+1j*Gh2}

class Autocoherent_Hill(Model):
    """
    Self-Consistent Model
    """
    def __init__(self):
        """
        Definition of model hypotheses.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = 'isotropic'  
        self.n_inclusions = 1 
        self.interphase = False 
        self.name = "Self-consistent"
        self.precision = 10**-12 ## Criterium of convergence of fixed-point algorithm
        self.n_point_fixe = 100 # Number of steps to reach final volumic fraction
    
    def Reccurence(Module,f):
        K,G,Km,Gm,Kf,Gf = Module
        ## Compute Kn+1
        numerator = f*(Kf-Km)*(3*K+4*G)
        denominator = 3*Kf+4*G
        nextK = Km + numerator/denominator
        ## Compute Gn+1
        numerator = 5*f*G*(Gf-Gm)*(3*K+4*G)
        denominator = 3*K*(3*G+2*Gf)+4*G*(3*Gf+2*G)        
        nextG = Gm + numerator/denominator
        return nextK,nextG
    
  
    def compute_behavior(self, Cm, inclusion_behaviors):
        # Retrieving matrix behavior
        Km = Cm['K']
        Gm = Cm['G']
        # Retrieving inclusion behavior
        Cf, f, ratio = inclusion_behaviors[0]
        Kf = Cf['K']
        Gf = Cf['G']
        F = np.linspace(0,f,self.n_point_fixe)
        
        Kinit = Km
        Ginit = Gm
        for i in range(len(F)) : 
            fi = F[i]
            # Initialisation of fixed-point algorithm
            K = Kinit
            G = Ginit
            # Fixed-point algorithm
            precision = self.precision
            nextK,nextG=Autocoherent_Hill.Reccurence([K,G,Km,Gm,Kf,Gf],fi)
            while abs((nextK-K)/K) > precision or abs((nextG-G)/G) > precision : 
                K,G=nextK,nextG
                nextK,NextG=Autocoherent_Hill.Reccurence([K,G,Km,Gm,Kf,Gf],fi)  
            # Updating initialisation
            Kinit = nextK
            Ginit = nextG
        return {'K': nextK, 'G': nextG}
    
class Autocoherent_III(Model):
    """
    Hypotheses : 
    -isotropic
    -spherical reinforcements 
    -elastic deformations 
    Self-coherent model. Contains :
    - A function that checks if the model is applicable to a given structure.
    - A model description function (TODO: write a function that returns a description of the model as a str and that could be called in the main)
    - A function that returns the homogenized behavior of the microstructure.
    - Functions that compute a particular characteristic (volume fraction of an inclusion, radius of an inclusion, behavior of an inclusion, etc.) from a target homogenized behavior (TODO).
    """
    
    def __init__(self):
        """
        Definition of model hypotheses.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = 'isotropic'  
        self.n_inclusions = 1 
        self.interphase = False 
        self.name = "Generalised self-consistent"
  
    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Computes the equivalent homogenized behavior of the microstructure. Returns a dict with the calculated parameters.
        """
        # Retrieving matrix behavior
        Km, Gm, num = Cm['K'], Cm['G'], Cm['nu']
        # Retrieving inclusion behavior
        Cf, f, ratio = inclusion_behaviors[0]
        Kf, Gf, nuf = Cf['K'], Cf['G'], Cf['nu']

        ## Useful constant to compute G         
        dm=(Gf/Gm)-1 
        eta1=dm*(49-50*nuf*num)+35*(dm+1)*(nuf-2*num)+35*(2*nuf-num) 
        eta2=5*nuf*(dm-7)+7*(dm+5) 
        eta3=(dm+1)*(8-10*num)+(7-5*num) 
        
        A=8*f**(10/3)*eta1*dm*(4-5*num)-2*f**(7/3)*(63*dm*eta2+2*eta1*eta3)+252*dm*eta2*f**(5/3)-50*dm*(7-12*num+8*num**2)*eta2*f+4*(7-10*num)*eta2*eta3 
        B=-4*dm*(1-5*num)*eta1*f**(10/3)+4*(63*dm*eta2+2*eta1*eta3)*f**(7/3)-504*dm*eta2*f**(5/3)+150*dm*(3-num)*num*eta2*f+3*(15*num-7)*eta2*eta3 
        D=4*dm*(5*num-7)*eta1*f**(10/3)-2*f**(7/3)*(63*dm*eta2+2*eta1*eta3)+252*dm*eta2*f**(5/3)+25*dm*(num**2-7)*eta2*f-(7+5*num)*eta2*eta3 
        
        ## Computation of G
        delta=B*B-4*D*A 
        sol1=(-B - delta**(1/2))/(2*A) 
        sol2=(-B + delta**(1/2))/(2*A) 
        sol=sol1 
        if ((sol1.real)<0) : 
            sol=sol2
            
        Gh=sol*Gm
        Kh=Km+f*(Kf-Km)/(1+(1-f)*(Kf-Km)/(Km+(4/3)*Gm))
        
        return {'K': Kh, 'G': Gh}

class Autocoherent_IV(Model):
    """
    Hypotheses : 
    Assumptions: 
    -isotropic
    -Spherical reinforcements 
    -elastic deformations 

    Self-consistent model. Contains :
    - A function that checks if the model is applicable to a microstructure.
    - A model description function 
    - A function that returns the homogenized behavior of the microstructure.
    - Functions that computes a particular characteristic (volume fraction of an inclusion, radius of an inclusion, behavior of an inclusion, etc.) from a target homogenized behavior 
    """
    def __init__(self, R_inclusion=1):
        """
       Definition of model hypotheses.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = 'isotropic'  
        self.n_inclusions = 1 
        self.interphase = True 
        self.R_inclusion = R_inclusion
        self.name = "4-phases self-consistent"

    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Calculates the equivalent homogenized behavior of the microstructure. Returns a dict with the calculated parameters. Currently, only calculates the shear modulus.
        TODO: complete with the complete calculation (K and G)
        """
        # Retrieving materials parameters
        Cf, f, ratio0 = inclusion_behaviors[0]
        Cv, cf, ratio1 = inclusion_behaviors[1]
        Km,Gm,num = Cm['K'], Cm['G'], Cm['nu']
        Kf,Gf,nuf = Cf['K'], Cf['G'], Cf['nu']
        Kv,Gv,nuv = Cv['K'], Cv['G'], Cv['nu']
        
        Rf = self.R_inclusion
        Rm = Rf/(f**(1/3))
        Rv = Rm*(f+cf)**(1/3)
        
        
        a1=(Gf/Gv)*(7+5*nuf)*(7-10*nuv)-(7-10*nuf)*(7+5*nuv) 
        b1=4*(7-10*nuf)+(Gf/Gv)*(7+5*nuf) 
        c1=(7-5*nuv)+2*(Gf/Gv)*(4-5*nuv) 
        d1=(7+5*nuv)+4*(Gf/Gv)*(7-10*nuv) 
        e1=2*(4-5*nuf)+(Gf/Gv)*(7-5*nuf) 
        f1=(4-5*nuf)*(7-5*nuv)-(Gf/Gv)*(4-5*nuv)*(7-5*nuf) 
        alpha1=(Gf/Gv)-1 
        
        a2=(Gv/Gm)*(7+5*nuv)*(7-10*num)-(7-10*nuv)*(7+5*num) 
        b2=4*(7-10*nuv)+(Gv/Gm)*(7+5*nuv) 
        c2=(7-5*num)+2*(Gv/Gm)*(4-5*num) 
        d2=(7+5*num)+4*(Gv/Gm)*(7-10*num) 
        e2=2*(4-5*nuv)+(Gv/Gm)*(7-5*nuv) 
        f2=(4-5*nuv)*(7-5*num)-(Gv/Gm)*(4-5*num)*(7-5*nuv) 
        alpha2=(Gv/Gm)-1 
    
        M1=np.zeros(shape=(4,4),dtype=complex)
        M1[0,0]=(5*(1-nuv))**(-1)*c1/3 
        M1[0,1]=(5*(1-nuv))**(-1)*Rf**2*(3*b1-7*c1)/(5*(1-2*nuf)) 
        M1[0,2]=(5*(1-nuv))**(-1)*(-12*alpha1/Rf**5) 
        M1[0,3]=(5*(1-nuv))**(-1)*4*(f1-27*alpha1)/(15*(1-2*nuf)*Rf**3) 
        M1[1,0]=0 
        M1[1,1]=(5*(1-nuv))**(-1)*(1-2*nuv)*b1/(7*(1-2*nuf)) 
        M1[1,2]=(5*(1-nuv))**(-1)*(-20*(1-2*nuv)*alpha1)/(7*Rf**7) 
        M1[1,3]=(5*(1-nuv))**(-1)*(-12*alpha1*(1-2*nuv))/(7*(1-2*nuf)*Rf**5) 
        M1[2,0]=(5*(1-nuv))**(-1)*Rf**5*alpha1/2 
        M1[2,1]=(5*(1-nuv))**(-1)*(-Rf**7*(2*a1+147*alpha1))/(70*(1-2*nuf)) 
        M1[2,2]=(5*(1-nuv))**(-1)*d1/7 
        M1[2,3]=(5*(1-nuv))**(-1)*Rf**2*(105*(1-nuv)+12*alpha1*(7-10*nuv)-7*e1)/(35*(1-2*nuf)) 
        M1[3,0]=(5*(1-nuv))**(-1)*(-5/6)*(1-2*nuv)*alpha1*Rf**3 
        M1[3,1]=(5*(1-nuv))**(-1)*7*(1-2*nuv)*alpha1*Rf**5/(2*(1-2*nuf)) 
        M1[3,2]=0 
        M1[3,3]=(5*(1-nuv))**(-1)*e1*(1-2*nuv)/(3*(1-2*nuf)) 
        
        M2=np.zeros(shape=(4,4),dtype=complex)
        M2[0,0]=(5*(1-num))**(-1)*c2/3 
        M2[0,1]=(5*(1-num))**(-1)*Rv**2*(3*b2-7*c2)/(5*(1-2*nuv)) 
        M2[0,2]=(5*(1-num))**(-1)*(-12*alpha2/Rv**5) 
        M2[0,3]=(5*(1-num))**(-1)*4*(f2-27*alpha2)/(15*(1-2*nuv)*Rv**3) 
        M2[1,0]=0 
        M2[1,1]=(5*(1-num))**(-1)*(1-2*num)*b2/(7*(1-2*nuv)) 
        M2[1,2]=(5*(1-num))**(-1)*(-20*(1-2*num)*alpha2)/(7*Rv**7) 
        M2[1,3]=(5*(1-num))**(-1)*(-12*alpha2*(1-2*num))/(7*(1-2*nuv)*Rv**5) 
        M2[2,0]=(5*(1-num))**(-1)*Rv**5*alpha2/2 
        M2[2,1]=(5*(1-num))**(-1)*(-Rv**7*(2*a2+147*alpha2))/(70*(1-2*nuv)) 
        M2[2,2]=(5*(1-num))**(-1)*d2/7 
        M2[2,3]=(5*(1-num))**(-1)*Rv**2*(105*(1-num)+12*alpha2*(7-10*num)-7*e2)/(35*(1-2*nuv)) 
        M2[3,0]=(5*(1-num))**(-1)*(-5/6)*(1-2*num)*alpha2*Rv**3 
        M2[3,1]=(5*(1-num))**(-1)*7*(1-2*num)*alpha2*Rv**5/(2*(1-2*nuv)) 
        M2[3,2]=0 
        M2[3,3]=(5*(1-num))**(-1)*e2*(1-2*num)/(3*(1-2*nuv)) 
        
        P = np.dot(M2,M1) 
        
        Z12 = P[0,0]*P[1,1]-P[1,0]*P[0,1] 
        Z14 = P[0,0]*P[3,1]-P[3,0]*P[0,1] 
        Z42 = P[3,0]*P[1,1]-P[1,0]*P[3,1] 
        Z23 = P[1,0]*P[2,1]-P[2,0]*P[1,1] 
        Z43 = P[3,0]*P[2,1]-P[2,0]*P[3,1] 
        Z13 = P[0,0]*P[2,1]-P[2,0]*P[0,1] 
    
        A = 4*Rm**10*(1-2*num)*(7-10*num)*Z12+20*Rm**7*(7-12*num+8*num**2)*Z42+12*Rm**5*(1-2*num)*(Z14-7*Z23)+20*Rm**3*(1-2*num)**2*Z13+16*(4-5*num)*(1-2*num)*Z43
        B = 3*Rm**10*(1-2*num)*(15*num-7)*Z12+60*Rm**7*(num-3)*num*Z42-24*Rm**5*(1-2*num)*(Z14-7*Z23)-40*Rm**3*(1-2*num)**2*Z13-8*(1-5*num)*(1-2*num)*Z43
        C = -Rm**10*(1-2*num)*(7+5*num)*Z12+10*Rm**7*(7-num**2)*Z42+12*Rm**5*(1-2*num)*(Z14-7*Z23)+20*Rm**3*(1-2*num)**2*Z13-8*(7-5*num)*(1-2*num)*Z43
        
        delta=B*B-4*C*A 
        sol1=(-B - delta**(1/2))/(2*A) 
        sol2=(-B + delta**(1/2))/(2*A) 
        sol=sol2 
        if (sol2.real<0):
            sol=sol1 
     
        Gh=sol*Gm
        X=(3*Km+4*Gm)*(f+cf)*( (Kf-Kv)*f*(3*Km+4*Gv)+(Kv-Km)*(cf+f)*(3*Kf+4*Gv)) 
        Y=3*(Kv-Kf)*f*( (f+cf)*(3 *Km+4*Gv)+4*(Gm-Gv)) + (3*Kf+4*Gv)*(f+cf)*(3*(cf+f)*(Km-Kv)+(3*Kv+4*Gm)) 
        Kh=Km+X/Y
        return {'K': Kh, 'G': Gh}
    
     
 
#%% Useful functions
def bulk_to_young(K, G):
    """
    Transforms  K and G modulus into E and nu
    """
    E = 9*K*G/(3*K+G)
    nu = (3*K-2*G)/(2*(3*K+G))
    return E, nu
   
def young_to_bulk(E, nu):
    """
    Transforms E and nu modulus into K and G 
    """
    K = E/(3*(1-2*nu))
    G = E/(2*(1+nu))
    return K, G   

def bulk_to_shear(K, E):
    """
    Transforms  K and E modulus into G and nu
    """
    G = 3*K*E/(9*K-E)
    nu = (3*K-E)/(6*K)
    return G, nu

def complete_behavior(behavior):
    """
    If the input behaviour is isotropic, completes it with E and nu or K and G.
    Otherwise, completes it with C or S if the input matrix is invertible.
    """
    parameters = list(behavior.keys())
    result = behavior
    nu_max = 0.4999999
    # Case of null (porous media) and incompressible (nu=0.5) media
    for parameter, values in result.items():
        # Isotropic elastic behavior
        if type(values) in [float, int]:
            if values==0:
                result[parameter] = 10**(-12)
            elif values==0.5 and parameter=='nu':
                result[parameter] = nu_max
        # Isotropic visco-elastic behavior
        elif type(values)==np.ndarray and type(values[0])!=np.ndarray:
            for i, value in enumerate(values):
                if value==0:
                    result[parameter][i] = 10**(-12)
                elif value==0.5 and parameter=='nu':
                    result[parameter][i] = nu_max
    # Isotropic K and G
    if parameters[:2]==['K', 'G'] or parameters[:2]==['G', 'K']:
        K, G = behavior['K'], behavior['G']
        E, nu = bulk_to_young(K, G)
        result['E'], result['nu'] = E, nu
    # Isotropic E and nu
    elif parameters[:2]==['E', 'nu'] or parameters[:2]==['nu', 'E']:
        E, nu = behavior['E'], behavior['nu']        
        K, G = young_to_bulk(E, nu)
        result['K'], result['G'] = K, G
    # Isotropic K and E
    elif parameters[:2]==['K', 'E'] or parameters[:2]==['E', 'K']:
        K, E = behavior['K'], behavior['E']        
        G, nu = bulk_to_shear(K, E)
        result['G'], result['nu'] = G, nu
    # Anisotropic
    elif parameters[0]=='C':
        C = behavior['C']
        try:
            S = np.linalg.inv(C)
        except:
            # C non invertible
            S = None
        result['S'] = S
    elif parameters[0]=='S':
        S = behavior['S']
        try:
            C = np.linalg.inv(S)
        except:
            # S non invertible
            C = None
        result['C'] = C
    
    # Return result
    return result


#%% Definition of model, behaviors et inclusion shape 
list_models = [Mori_Tanaka, Differential_Scheme, Autocoherent_Hill, Autocoherent_III, Autocoherent_IV]
dict_behaviors_visco = {'Elastic isotropic (K & G)': ['K', 'G'],
                        'Elastic isotropic (E & nu)': ['E', 'nu'],
                        'Visco-elastic 1': ['K', "G'", "G''"],
                        'Visco-elastic 2': ["K'", "K''", "G'", "G''"],
                        'Visco-elastic 3': ["E'", "E''", 'K']}
dict_behaviors = {'Isotropic (K & G)': ['K', 'G'],
                  'Isotropic (E & nu)': ['E', 'nu'],
                  'Anisotropic (Compliance)': ['S'],
                  'Anisotropic (Stifness)': ['C']}
dict_types = {0: 'Spheres', 1: 'Ellipsoids'}
    
    
    
    
    
    
