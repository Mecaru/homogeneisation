"""
Homogeneisation - classes.py

Définition des classes utilisées.

Les items à modifier/améliorer/ajouter sont marqués en commentaires précédés de la mention "TODO".

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

27/04/2020
"""

import numpy as np
from scipy import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class Inclusion:
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, behavior, aspect_ratio=1, name=None, frequency=[]):
        """
        TODO : Prise en compte de l'orientation
        type_inclusion : (int), 0 pour des inclusions sphériques. Voir la liste list_types (en bas du fichier) pour les autres types.
        radius : (float), valeur du rayon des inclusions sphériques. TODO : À remplacer par un paramètre plus général pour des inclusions de types différents. 
        behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). TODO :  À modifier pour représenter des comportements non isotropes.
        frequency: (list), liste des fréquences/températures associées aux paramètres visco-élastiques
        """
        self.type_inclusion = type_inclusion
        self.aspect_ratio = aspect_ratio
        self.behavior = complete_behavior(behavior)
        self.name = name
        self.frequency = frequency
    
    def type_to_str(self):
        """
        Transforme un entier "type_inclusion" en la chaîne de caractères correspondante (exemple : 0 --> "spheres") 
        TODO : synchroniser cette fonction avec le main à l'aide d'un dictionnaire pour faciliter l'ajout de types d'inclusions
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
        Présentation de l'instance.
        """
        str_type_inclusion = self.type_to_str()
        string = "{}, {}".format(self.name, str_type_inclusion)
        if self.type_inclusion != 0:
            string += " (c={})".format(self.aspect_ratio)
        for parameter, value in self.behavior.items():
            if type(value) not in [list, np.ndarray]:
                string += ", {}: {:.2f}".format(parameter, value)
            else:
                string += ", {}: Visco-elastic".format(parameter)
        return string

    def __repr__(self):
        return str(self)

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
        Trace le graphe d'évolution des paramètres visqueux-élastiques si ceux-ci existent.
        """
        if self.frequency == []:
            None # L'inclusion ne contient pas de paramètres visco-élastiques
        else:
            plt.figure()
            for parameter, values in self.behavior.items():
                if type(values) == list:
                    # Le paramètre est visco-élastique
                    plt.plot(self.frequency, values, '.', label=parameter)
                    plt.legend()
                    plt.xlabel("Frequency/Temperature")
                    plt.ylabel("Parameter value")
                    plt.title("Inclusion visco-elastic behavior")
                    plt.xlim(min(self.frequency), max(self.frequency))
            plt.show()
    
class Microstructure:
    """
    Contient des informations sur la microstructure (comportement de la matrice, inclusions, etc..). TODO : à modifier pour prendre en compte la présence ou non d'une interphase, et d'autres paramètres de modèles plus avancés.
    Contient une fonction qui renvoie les bornes de Hashin-Shtrickman pour la microstructure en question 
    TODO : Généraliser ces bornes à n phases (et pas 2 comme c'est le cas ici)
    """
    
    def __init__(self, matrix_behavior, dict_inclusions=dict(), frequency=[]):
        """
        list_inclusions : (dict), sous la forme {inclusion: f_i} avec inclusion une instance de classe Inclusion et f_i la fraction volumique de ce type d'inclusion.
        matrix_behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). TODO :  À modifier pour représenter des comportements non isotropes.
        frequency: liste des fréquences associées aux paramètres visco-élastiques
        """
        self.dict_inclusions = dict_inclusions
        self.matrix_behavior = complete_behavior(matrix_behavior)
        # Calcul de la fraction volumique de matrice f_m
        self.f_matrix = self.compute_fm()
        self.frequency = frequency
        
    def __str__(self):
        string = "Microstructure\nf_m = {:.2f}, matrix".format(self.f_matrix, self.matrix_behavior)
        for parameter, value in self.matrix_behavior.items():
            if type(value) not in [list, np.ndarray]:
                string += ", {}: {:.2f}".format(parameter, value)
            else:
                string += ", {}: Visco-elastic".format(parameter)
                # TODO : transformer cette ligne en fonction print_behavior et l'appeler lors de l'affichage du comportement homogénéisé et des bornes de hashin
        dict_inclusions = self.dict_inclusions
        # Présentation de toutes les inclusions contenues dans la microstructure
        for inclusion in dict_inclusions.keys():
            fi = dict_inclusions[inclusion]
            string += "\nf_i = {}, ".format(fi) + str(inclusion)
        return string

    def compute_fm(self):
        """
        1/ Vérifie si la liste des inclusions donnée est cohérente (i.e : la somme des fractions volumiques des inclusions est inférieure à 1). Si ce n'est pas le cas, génère une erreur.
        2/ Si aucune erreur n'est générée, calcule la fraction volumique de matrice.
        """
        total_fi = 0 # Total des fractions volumiques d'inclusions
        dict_inclusions = self.dict_inclusions
        for inclusion in dict_inclusions.keys():
            fi = dict_inclusions[inclusion]
            total_fi += fi
        if total_fi >= 1:
            raise NameError("The total volumic fractions of the inclusions exceed 1")
        else :
            f_m = 1 - total_fi
            return f_m

    def change_fi(self, inclusion, new_f):
        """
        Met à jour la fraction volumique de l'inclusion ou l'ajoute au dictionnaire si celle-ci n'y était pas présente.
        Met à jour la fraction volumique de matrice.
        """
        self.dict_inclusions[inclusion] = new_f
        self.f_matrix = self.compute_fm()
    
    def change_parameter(self, parameter, new_value):
        """
        Change the value of the parameter if it exists. Updates the behavior with the function "complete_behavior".
        """
        try:
            self.matrix_behavior[parameter] = new_value
            self.matrix_behavior = complete_behavior(self.matrix_behavior)
        except:
            None

    def graph_parameter(self):
        """
        Trace le graphe d'évolution des paramètres visqueux-élastiques si ceux-ci existent.
        """
        if self.frequency == []:
            None # L'inclusion ne contient pas de paramètres visco-élastiques
        else:
            plt.figure()
            for parameter, values in self.matrix_behavior.items():
                if type(values) == list:
                    # Le paramètre est visco-élastique
                    plt.plot(self.frequency, values, '.', label=parameter)
                    plt.legend()
                    plt.xlabel("Frequency/Temperature")
                    plt.ylabel("Parameter value")
                    plt.title("Matrix visco-elastic behavior")
                    plt.xlim(min(self.frequency), max(self.frequency))
            plt.show()

    def draw(self):
        """
        Méthode qui permet de dessiner la microstructure. Pour le moment, fonctionne uniquement avec une seule inclusion, sphérique, oblate ou prolate.
        """
        inclusions = list(self.dict_inclusions.keys())
        if len(inclusions) == 1:
            inclusion = inclusions[0]
            fi = self.dict_inclusions[inclusion]
            # Calcul du rayon pour un VER de taille 10X10X10
            ratio = inclusion.aspect_ratio
            a = (1000*fi/(4/3*pi*ratio))**(1/3)
            b = ratio*a
            
            fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
            ax = fig.add_subplot(111, projection='3d')

            # Radii:
            rx, ry, rz = np.array([b, a, a])

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
            
     ## CALCUL DES BORNES DE HASHIN-SHTRICKMAN ##########  
    
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
        Donne les bornes de Hashin-Shtrikman pour 1 seule phase, isotrope
        TODO : ajouter le cas des inclusion multiples
        """
        fm=self.f_matrix
        f=1-fm
        km,gm=self.matrix_behavior["K"],self.matrix_behavior["G"]
        
        for inclusion in self.dict_inclusions.keys():
            kf,gf=inclusion.behavior["K"],inclusion.behavior["G"]
        
        ksup=max(Microstructure.khs(km,gm,fm,kf,gf,f),Microstructure.khs(kf,gf,f,km,gm,fm))
        kinf=min(Microstructure.khs(km,gm,fm,kf,gf,f),Microstructure.khs(kf,gf,f,km,gm,fm))
        gsup=max(Microstructure.ghs(km,gm,fm,kf,gf,f),Microstructure.ghs(kf,gf,f,km,gm,fm))
        ginf=min(Microstructure.ghs(km,gm,fm,kf,gf,f),Microstructure.ghs(kf,gf,f,km,gm,fm))
            
        
        return { 'Ginf' : ginf, 'Gsup' : gsup, 'Kinf' : kinf, 'Ksup' : ksup }


    
class Mori_Tanaka:
    """
    TODO : vérifier si le modèle s'applique aussi à d'autres types d'inclusion, pour le moment seules des inclusions sphériques isotropes sont prises en compte pour tester le code.
    Modèle de Mori-Tanaka. Contient :
    - Une fonction qui vérifie si le modèle est appliquable à une microstructure.
    - Une fonction de description du modèle (TODO : écrire une fonction qui renvoie une description du modèle sous forme de str et qui pourrait être appelée dans le main)
    - Un fonction qui renvoie le comportement homogénéisé de la microstructure.
    - Des fonctions qui calculent une caractéristique particulière (fraction volumique d'une inclusion, rayon d'une inclusion, comportement d'une inclusion, etc..) à partir d'un comportement homogénéisé cible (TODO)
    """
    
    def __init__(self):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = set(['K', "G'", "G''", 'G', 'E', 'nu'])  # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes, éventuellement visco-élastiques
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.name = "Mori-Tanaka"
        
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return "Modèle de Mori-Tanaka"
    
    def __repr__(self):
        """
        Description textuelle du modèle.
        """
        return str(self)
    
    def check_hypothesis(self, microstructure):
        """
        Vérifies si la microstructure vérifie les hypothèses du modèle, renvoie un boolées. 
        TODO : Éventuellement généraliser cette fonction en l'incorporant dans une classe mère Model pour qu'elle s'applique à tous les modèles.
        """
        dict_inclusions = microstructure.dict_inclusions
        inclusions = dict_inclusions.keys()
        n_inclusions = len(inclusions)
        # vérification du nombre d'inclusions
        if n_inclusions != self.n_inclusions:
            # Le modèle ne peut pas traiter de microstructures avec autant d'inclusions de natures différentes
             #raise NameError("Wrong number of inclusion")
             return False
        for inclusion in dict_inclusions.keys():
            # Vérification du type d'inclusion
            if inclusion.type_inclusion != self.type_inclusion:
                #raise NameError("Wrong type of inclusion or microstructure")
                return False
            # vérification du comportement des inclusions
            behavior = inclusion.behavior
            if set(behavior.keys()).issubset(self.behavior_condition) == False:
                #print (list(behavior.keys()) , self.behavior_condition)
                #raise NameError("Inclusion and microstructure behavior incompatible")
                return False
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()).issubset(self.behavior_condition) == False:
            raise NameError("Inclusion and microstructure behavior incompatible")
            return False
        # À ce stade, toutes les conditions ont été vérifiées
        return True
    
    def compute_h_behavior(self, microstructure, viscoelastic=False):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis(microstructure)
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        Cm = microstructure.matrix_behavior
        dict_inclusions = microstructure.dict_inclusions
        inclusion = list(dict_inclusions.keys())[0] #Inclusion unique ici
        Cf = inclusion.behavior
        Gm, Km = Cm['G'], Cm['K']
        Gf, Kf = Cf['G'], Cf['K']
        f = dict_inclusions[inclusion]
        
        denominator = 5*Gm*(3*Km+4*Gm)+6*(1-f)*(Gf-Gm)*(Km+2*Gm)
        numerator = 5*f*Gm*(Gf-Gm)*(3*Km+4*Gm)
        Gh = Gm + numerator/denominator
        
        denominator = 3*Kf+4*Gm+3*(1-f)*(Kf-Km)
        numerator = f*(Kf-Km)*(3*Km+4*Gm)
        Kh = Km + numerator/denominator
        return complete_behavior({'K' : Kh, 'G' : Gh}) 
    

    
    def check_bounds(self,microstructure):
        Behavior_h=self.compute_h_behavior(microstructure)
        print(Behavior_h)
        Gh = Behavior_h['G']
        Kh = Behavior_h['K']
        Bounds=microstructure.Hashin_bounds()
        Gsup = Bounds['Gsup']
        Ginf = Bounds['Ginf']
        Ksup = Bounds['Ksup']
        Kinf = Bounds['Kinf']
        if Gh < Ginf or Gh > Gsup : 
            raise NameError("G out of Hashin-Shtrikman bounds")
            return False
        if Kh < Kinf or Kh > Ksup :
            raise NameError("K out of Hashin-Shtrikman bounds")
            return False
        return True
    
class Eshelby_Approximation:
    """
    TODO : vérifier si le modèle s'applique aussi à d'autres types d'inclusion, pour le moment seules des inclusions sphériques isotropes sont prises en compte pour tester le code.
    Modèle de Mori-Tanaka. Contient :
    - Une fonction qui vérifie si le modèle est appliquable à une microstructure.
    - Une fonction de description du modèle (TODO : écrire une fonction qui renvoie une description du modèle sous forme de str et qui pourrait être appelée dans le main)
    - Un fonction qui renvoie le comportement homogénéisé de la microstructure.
    - Des fonctions qui calculent une caractéristique particulière (fraction volumique d'une inclusion, rayon d'une inclusion, comportement d'une inclusion, etc..) à partir d'un comportement homogénéisé cible (TODO)
    """
    
    def __init__(self):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0
        self.behavior_condition = set(['K', "G'", "G''", 'G', 'E', 'nu']) # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.name = "Eshelby"
        
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return "Modèle d'Eshelby'"
    
    def __repr__(self):
        """
        Description textuelle du modèle.
        """
        return str(self)
    
    def check_hypothesis(self, microstructure):
        """
        Vérifies si la microstructure vérifie les hypothèses du modèle, renvoie un boolées. 
        TODO : Éventuellement généraliser cette fonction en l'incorporant dans une classe mère Model pour qu'elle s'applique à tous les modèles.
        """
        dict_inclusions = microstructure.dict_inclusions
        inclusions = dict_inclusions.keys()
        n_inclusions = len(inclusions)
        # vérification du nombre d'inclusions
        if n_inclusions != self.n_inclusions:
            # Le modèle ne peut pas traiter de microstructures avec autant d'inclusions de natures différentes
             #raise NameError("Wrong number of inclusion")
             return False
        for inclusion in dict_inclusions.keys():
            # Vérification du type d'inclusion
            if inclusion.type_inclusion != self.type_inclusion:
                #raise NameError("Wrong type of inclusion or microstructure")
                return False
            # vérification du comportement des inclusions
            behavior = inclusion.behavior
            if set(behavior.keys()).issubset(self.behavior_condition) == False:
                #print (list(behavior.keys()) , self.behavior_condition)
                #raise NameError("Inclusion and microstructure behavior incompatible")
                return False
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()).issubset(self.behavior_condition) == False:
            raise NameError("Inclusion and microstructure behavior incompatible")
            return False
        # À ce stade, toutes les conditions ont été vérifiées
        return True
    
    def compute_h_behavior(self, microstructure):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis(microstructure)
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        Cm = microstructure.matrix_behavior
        dict_inclusions = microstructure.dict_inclusions
        inclusion = list(dict_inclusions.keys())[0] #Inclusion unique ici
        Cf = inclusion.behavior
        Gm, Km = Cm['G'], Cm['K']
        Gf, Kf = Cf['G'], Cf['K']
        f = dict_inclusions[inclusion]
        
        denominator = 3*Km*(3*Gm+2*Gf) + 4*Gm*(2*Gm+3*Gf)
        numerator = 5*f*Gm*(Gf-Gm)*(3*Km+4*Gm)
        Gh = Gm + numerator/denominator

        
        denominator = 3*Kf+4*Gm
        numerator = f*(Kf-Km)*(3*Km+4*Gm)
        Kh = Km + numerator/denominator
        
        return complete_behavior({'K' : Kh, 'G' : Gh}) 
    

    
    def check_bounds(self,microstructure):
        Behavior_h=self.compute_h_behavior(microstructure)
        print(Behavior_h)
        Gh = Behavior_h['G']
        Kh = Behavior_h['K']
        Bounds=microstructure.Hashin_bounds()
        Gsup = Bounds['Gsup']
        Ginf = Bounds['Ginf']
        Ksup = Bounds['Ksup']
        Kinf = Bounds['Kinf']
        if Gh < Ginf or Gh > Gsup : 
            raise NameError("G out of Hashin-Shtrikman bounds")
            return False
        if Kh < Kinf or Kh > Ksup :
            raise NameError("K out of Hashin-Shtrikman bounds")
            return False
        return True
    

class Differential_Scheme:
    """
    TODO : 
    vérifier si le modèle s'applique aussi à d'autres types d'inclusion, pour le moment seules des inclusions sphériques isotropes sont prises en compte pour tester le code. 
    Modèle différentiel. Contient :
    - Une fonction qui vérifie si le modèle est appliquable à une microstructure.
    - Une fonction de description du modèle (TODO : écrire une fonction qui renvoie une description du modèle sous forme de str et qui pourrait être appelée dans le main)
    - Un fonction qui renvoie le comportement homogénéisé de la microstructure.
    - Une fonction qui vérifie que le comportement homogénéisé se trouve dans les bornes de Hashin-Shtrickman.
    """
    
    def __init__(self):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0
        self.behavior_condition = set(['K', "G'", "G''", 'G', 'E', 'nu']) # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.name = "Differential"
        
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return "Modèle différentiel"
    
    def __repr__(self):
        """
        Description textuelle du modèle.
        """
        return str(self)
    
    def check_hypothesis(self, microstructure):
        """
        Vérifies si la microstructure vérifie les hypothèses du modèle, renvoie un boolées. 
        TODO : Éventuellement généraliser cette fonction en l'incorporant dans une classe mère Model pour qu'elle s'applique à tous les modèles.
        """
        dict_inclusions = microstructure.dict_inclusions
        inclusions = dict_inclusions.keys()
        n_inclusions = len(inclusions)
        # vérification du nombre d'inclusions
        if n_inclusions != self.n_inclusions:
            # Le modèle ne peut pas traiter de microstructures avec autant d'inclusions de natures différentes
             #raise NameError("Wrong number of inclusion")
             return False
        for inclusion in dict_inclusions.keys():
            # Vérification du type d'inclusion
            if inclusion.type_inclusion != self.type_inclusion:
                #raise NameError("Wrong type of inclusion or microstructure")
                return False
            # vérification du comportement des inclusions
            behavior = inclusion.behavior
            if set(behavior.keys()).issubset(self.behavior_condition) == False:
                #print (list(behavior.keys()) , self.behavior_condition)
                #raise NameError("Inclusion and microstructure behavior incompatible")
                return False
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()).issubset(self.behavior_condition) == False:
            raise NameError("Inclusion and microstructure behavior incompatible")
            return False
        # À ce stade, toutes les conditions ont été vérifiées
        return True
    
    ## Fonctions utiles au calcul du comportement homogénéisé
    
    def deriv(Module,f):
        K,G,Kf,Gf=Module[0],Module[1],Module[2],Module[3]
        mu=(3*K-2*G)/(6*K+2*G)
        
        numerator=K-Kf
        denominator=(1-f)*(1+(Kf-K)/(K+4*G/3))
        dK=-numerator/denominator
        
        numerator=15*(1-mu)*(G-Gf)
        denominator=(1-f)*(7-5*mu+2*(4-5*mu)*Gf/G)
        dG=-numerator/denominator
        
        return np.array([dK,dG,0,0])
    
    def khs(k1, g1, c1, k2, g2, c2):
        numerator = c2*(k2-k1)
        denominator = 1+3*c1*(k2-k1)/(4*g1+3*k1)
        return k1+numerator/denominator
    
    def ghs(k1, g1, c1, k2, g2, c2):
        numerator = c2*(g2-g1)
        denominator = 1+6*c1*(g2-g1)*(k1+2*g1)/((3*k1+4*g1)*5*g1)
        return g1+numerator/denominator
    
    def compute_h_behavior(self, microstructure):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis(microstructure)
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        Cm = microstructure.matrix_behavior
        dict_inclusions = microstructure.dict_inclusions
        inclusion = list(dict_inclusions.keys())[0] #Inclusion unique ici
        Cf = inclusion.behavior
        Gm, Km = Cm['G'], Cm['K']
        Gf, Kf = Cf['G'], Cf['K']
        f_finale = dict_inclusions[inclusion]
        
        npoints=100
        f=np.linspace(0,f_finale,npoints)
        Module_Initial=np.array([Km,Gm,Kf,Gf])
        Module=odeint(Differential_Scheme.deriv,Module_Initial,f)
        
        Module_final=Module[npoints-1]
        Kh,Gh,Kf,Gf=Module_final   
        
        ## ajout des bornes de Hashin
        Khs=np.vectorize(Differential_Scheme.khs)
        Ghs=np.vectorize(Differential_Scheme.ghs)
        KINF=Khs(Km,Gm,1-f,Kf,Gf,f)
        GINF=Ghs(Km,Gm,1-f,Kf,Gf,f)
        KSUP=Khs(Kf,Gf,f,Km,Gm,1-f)
        GSUP=Ghs(Kf,Gf,f,Km,Gm,1-f)
        ## affichage des graphes pour K et G
        
        #plt.subplot(211)
        #plt.plot(f,Module[:,0],label="Kh")        
        #plt.plot(f,Module[:,2],label="Kf")
        #plt.plot(f,KSUP,label="Ksup")
        #plt.plot(f,KINF,label="Kinf")
        #plt.title("Kh en fonction de f ")
        #plt.legend()
        #plt.xlabel("fraction volumique")
        #plt.ylabel("Modules")
        
        #plt.subplot(212)
        #plt.title("Gh en fonction de f ")
        #plt.plot(f,Module[:,1],label="Gh")
        #plt.plot(f,Module[:,3],label="Gf")
        #plt.plot(f,GSUP,label="Gsup")
        #plt.plot(f,GINF,label="Ginf")
        #plt.legend()
        #plt.xlabel("fraction volumique")
        #plt.ylabel("Modules")
        #plt.show()
        return complete_behavior({'K' : Kh, 'G' : Gh}) 
    

    
    def check_bounds(self,microstructure):
        Behavior_h=self.compute_h_behavior(microstructure)
        Gh = Behavior_h['G']
        Kh = Behavior_h['K']
        Bounds=microstructure.Hashin_bounds()
        Gsup = Bounds['Gsup']
        Ginf = Bounds['Ginf']
        Ksup = Bounds['Ksup']
        Kinf = Bounds['Kinf']
        if Gh < Ginf or Gh > Gsup : 
            raise NameError("G out of Hashin-Shtrikman bounds")
            return False
        if Kh < Kinf or Kh > Ksup :
            raise NameError("K out of Hashin-Shtrikman bounds")
            return False
        return True

################ Tests 

#inclusion1 = Inclusion(0, {"K":30, "G":150}, 1)
#inclusion2 = Inclusion(0, {"K":90, "G":150}, 1)
#inclusion3 = Inclusion(0, {"K":150, "G":150}, 1)
#inclusion4 = Inclusion(0, {"K":230, "G":150}, 1)
#inclusion5 = Inclusion(0, {"K":300, "G":150}, 1)
#Incl=[inclusion1,inclusion2,inclusion3,inclusion4,inclusion5]
#for i in range(1):
#    inclusion=inclusion5
#    f=0.999
#    microstructure = Microstructure({"K":30, "G":15}, {inclusion:f})
#    #print("Hashin Bounds : ", microstructure.Hashin_bounds())
#    model = Differential_Scheme()
#    Ch=model.compute_h_behavior(microstructure)
    #print("Comportement homogénéisé : ", model.compute_h_behavior(microstructure))
    #print ("Dans les bornes de Hashin : ", model.check_bounds(microstructure))
#    print(inclusion.behavior['K'],inclusion.behavior['G'],inclusion.radius,f,microstructure.matrix_behavior['K'],microstructure.matrix_behavior['G'],microstructure.Hashin_bounds()['Kinf'],microstructure.Hashin_bounds()['Ksup'],microstructure.Hashin_bounds()['Ginf'],microstructure.Hashin_bounds()['Gsup'],Ch['K'],1,Ch['G'])


def bulk_to_young(K, G):
    """
    Transforme des modules K et G en modules E et nu.
    """
    E = 9*K*G/(3*K+G)
    nu = (3*K-2*G)/(2*(3*K+G))
    return E, nu
   
def young_to_bulk(E, nu):
    """
    Transforme des modules E et nu en modules K et G
    """
    K = E/(3*(1-2*nu))
    G = E/(2*(1+nu))
    return K, G
    
def complete_behavior(behavior):
    """
    Si le comportement en entrée est isotrope, le complète avec E et nu ou K et G. Sinon, le renvoie tel quel.
    """
    parameters = list(behavior.keys())
    result = behavior
    if parameters[:2] == ['K', 'G']:
        K, G = behavior['K'], behavior['G']
        E, nu = bulk_to_young(K, G)
        result['E'], result['nu'] = E, nu
    elif parameters[:2] == ['E', 'nu']:
        E, nu = behavior['E'], behavior['nu']
        K, G = young_to_bulk(E, nu)
        result['K'], result['G'] = K, G
    return result
    
list_models = [Mori_Tanaka, Eshelby_Approximation, Differential_Scheme] # Liste des modèles implémentés, à incrémenter à chaque ajout d'un nouveau modèle
dict_behaviors = {'Elastic isotropic (K & G)': ['K', 'G'], 'Elastic isotropic (E & nu)': ['E', 'nu'], 'Visco-elastic': ['K', "G'", "G''"]}
dict_types = {0: 'Spheres', 1: 'Oblate', 2: 'Prolate'} # Types de géométries admissibles et leur identifiant

# Tests
# inclusion1 = Inclusion(1, {"E":300, "nu":0.3})
#print(inclusion1)
#inclusion1 = Inclusion(0, {"K":300, "G":0.3})
#print(inclusion1)
#inclusion2 = Inclusion(0, {"K":300, "G":150})
# microstructure = Microstructure({"E":10, "nu":0.1}, {inclusion1:0.6})
#model = Mori_Tanaka()
# print(microstructure)
#print(model.check_hypothesis(microstructure))
#print(model.compute_h_behavior(microstructure))
# microstructure.change_fi(inclusion1, 0.3)
# print(microstructure)
#microstructure.draw()
