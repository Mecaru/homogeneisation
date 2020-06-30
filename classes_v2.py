# -*- coding: utf-8 -*-
"""
Homogeneisation - classes.py

Définition des classes utilisées. Version 2 de classes, avec une définition simplifiée des modèles et prise en compte de la visco-élasticité.

Les items à modifier/améliorer/ajouter sont marqués en commentaires précédés de la mention "TODO".

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

11/06/2020
"""

#%% Importation des modules
import numpy as np
from scipy import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#%% Classes microstructure
class Inclusion:
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, behavior, aspect_ratio=[1.,1.], name=None, frequency=[], abscissa="frequency"):
        """
        TODO: Prise en compte de l'orientation
        type_inclusion: (int), 0 pour des inclusions sphériques, 1 pour des inclusions ellipsoïdales
        aspect_ratio: (tuple), tuple de deux flottants représentant les rapports des longueurs des axes 2 et 3 de l'ellipsoïde sur la longueur de l'axe 1
        behavior: (dict), contient les valeurs des paramètres de la matrice de comportement, voir dict_behavior dans la dernière section du script
        frequency: (list), liste des fréquences/températures associées aux paramètres visco-élastiques
        abscissa: (str), vaut "frequency" ou "temperature", indique la nature physique des valeurs de la liste frequency
        inc_and_int: ([InclusionAndInterphase,int]), None par défaut, renvoie vers l'instance de classe InclusionAndInterphase à laquelle appartient l'inclusion, si celle-ci existe, et un entier, 0 pour l'inclusion, 1 pour l'interphase
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
        Change the value of the parameter if it exists. Updates the behavior with the function "complete_behavior".
        """
        try:
            self.behavior[parameter] = new_value
            self.behavior = complete_behavior(self.behavior)
        except:
            None

    def graph_parameter(self):
        """
        Trace le graphe d'évolution des paramètres visco-élastiques si ceux-ci existent.
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
    Instance représentant une inclusion et l'interphase associée.
    """            
    
    def __init__(self, inclusion, interphase, name=None):
        """
        inclusion: Instance de classe inclusion
        interphase: Instance de classe inclusion, représente l'interphase associée à l'inclusion
        L'interphase et l'inclusion doivent être du même type (sphères ou ellipsoïdes avec les mêmes rapports d'aspect)
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
    Contient des informations sur la microstructure (comportement de la matrice, inclusions, etc..).
    Contient une fonction qui renvoie les bornes de Hashin-Shtrickman pour la microstructure en question 
    TODO: Généraliser ces bornes à n phases (et pas 2 comme c'est le cas ici)
    """
    
    def __init__(self, behavior, dict_inclusions=dict(), frequency=[], abscissa="frequency"):
        """
        list_inclusions: (dict), sous la forme {inclusion: f_i} avec inclusion une instance de classe Inclusion et f_i la fraction volumique de ce type d'inclusion. inclusion peut aussi être une instance de classe InclusionAndInterphase, dans ce cas, f_i est un tuple de flottants
        behavior: (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment
        frequency: liste des fréquences associées aux paramètres visco-élastiques
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
            # Cas des inclusions + interphase
            if type(fi)==list:
                total_fi += fi[0] + fi[1]
            # Inclusions simples
            else:
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
            self.behavior[parameter] = new_value
            self.behavior = complete_behavior(self.behavior)
        except:
            None

    def graph_parameter(self):
        """
        Trace le graphe d'évolution des paramètres visco-élastiques si ceux-ci existent.
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
                plt.title("Matrix visco-elastic behavior")
                plt.xlim(min(self.frequency), max(self.frequency))
            plt.show()

    # def draw(self):
    #     """
    #     Méthode qui permet de dessiner la microstructure.
    #     """
    #     inclusions = list(self.dict_inclusions.keys())
    #     if len(inclusions) == 1:
    #         inclusion = inclusions[0]
    #         fi = self.dict_inclusions[inclusion]
    #         # Calcul du rayon pour un VER de taille 10X10X10
    #         c1, c2 = inclusion.aspect_ratio
    #         a = (1000*fi/(4/3*pi*c1*c2))**(1/3)
    #         b = c1*a
    #         c = c2*a
            
    #         fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
    #         ax = fig.add_subplot(111, projection='3d')

    #         # Radii:
    #         rx, ry, rz = np.array([a, b, c])

    #         # Set of all spherical angles:
    #         u = np.linspace(0, 2 * np.pi, 100)
    #         v = np.linspace(0, np.pi, 100)

    #         # Cartesian coordinates that correspond to the spherical angles:
    #         # (this is the equation of an ellipsoid):
    #         x = rx * np.outer(np.cos(u), np.sin(v))
    #         y = ry * np.outer(np.sin(u), np.sin(v))
    #         z = rz * np.outer(np.ones_like(u), np.cos(v))

    #         # Plot:
    #         ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')

    #         # Adjustment of the axes, so that they all have the same span:
    #         max_radius = 5
    #         for axis in 'xyz':
    #             getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))

    #         # Cube 
    #         points = 5*np.array([[-1, -1, -1],
    #                               [1, -1, -1 ],
    #                               [1, 1, -1],
    #                               [-1, 1, -1],
    #                               [-1, -1, 1],
    #                               [1, -1, 1 ],
    #                               [1, 1, 1],
    #                               [-1, 1, 1]])

    #         r = [-5,5]
    #         X, Y = np.meshgrid(r, r)
    #         one = 5*np.ones(4).reshape(2, 2)
    #         ax.plot_wireframe(X,Y,one, alpha=0.5)
    #         ax.plot_wireframe(X,Y,-one, alpha=0.5)
    #         ax.plot_wireframe(X,-one,Y, alpha=0.5)
    #         ax.plot_wireframe(X,one,Y, alpha=0.5)
    #         ax.plot_wireframe(one,X,Y, alpha=0.5)
    #         ax.plot_wireframe(-one,X,Y, alpha=0.5)
    #         ax.scatter3D(points[:, 0], points[:, 1], points[:, 2])

    #         plt.show()
            
    def draw(self):
        """
        Méthode qui permet de dessiner la microstructure.
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
                
            ### Tracé inclusion
            ax = fig.add_subplot(1, n_fig, index+1, projection='3d')
            # Calcul du rayon pour un VER de taille 10X10X10
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
            
            ### Tracé interphase
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
            
            ### Tracé bords VER
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

#%% Classes modèles
class Model:
    """
    Classe générique mère de toutes les classes modèles. 
    Contient la méthode permettant de vérifier les hypothèses du modèle sur une microstructure, ainsi que la méthode appelée lors du calcul du comportement homogénéisé.
    """
    
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return self.name + " model"
    
    def __repr__(self):
        """
        Description textuelle du modèle.
        """
        return str(self)
    
    def check_hypothesis(self, microstructure):
        """
        Vérifies si la microstructure vérifie les hypothèses du modèle, renvoie un booléen. 
        """
        # Récupération des inclusions de la microstructure
        dict_inclusions = microstructure.dict_inclusions
        instances = list(dict_inclusions.keys())
        n_instances = len(instances)
        # Initialisation du résultat
        result = True
        # Vérification du nombre d'inclusions
        if n_instances != self.n_inclusions:
             result = False
        # Vérification de la présence ou de l'absence d'interphase
        for instance in instances:
            if (type(instance)==InclusionAndInterphase)!=self.interphase:
                result = False
        # Construction d'une liste d'inclusions sans interphase
        inclusions = []
        for instance in instances:
            if type(instance)==InclusionAndInterphase:
                inclusions += [instance.inclusion, instance.interphase]
            else:
                inclusions.append(instance)
        # Vérification du type d'inclusion
        for inclusion in inclusions:
            if inclusion.type_inclusion != self.type_inclusion:
                result = False
        # Vérification du comportement des inclusions et de la matrice
        for element in inclusions + [microstructure]:
            if not set(element.behavior.keys()).issubset(self.behavior_condition):
                result = False
        # Renvoi du résultat
        return result
    
    def compute_h_behavior(self, microstructure):
        """
        Calcule le comportement homogénéisé de la microstructure avec le modèle.
        Vérifies que le modèle s'applique bien sur la microstructure entrée.
        Si les éléments en entrée ne sont pas visco-élastiques (i.e: si la liste frequency de la microstructure est vide), renvoie un dictionnaire de paramètres réels sous la forme {"parameter": value(float)}. Dans le cas isotrope, calcule aussi les paramètres manquants (mu et E, ou K et G).
        Sinon, réalise une boucle sur les valeurs de fréquence, puis renvoie un dictionnaire de comportement complet (avec valeurs des paramètres manquants) de la forme {"parameter": [values(complex)]}.
        """
        # Vérification des conditions d'application
        compatible = self.check_hypothesis(microstructure)
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        frequency = microstructure.frequency
        # Création du dictionnaire contenant uniquement des instances de classe Inclusion
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
            # Récupération du comportement des inclusions, format [(inclusion.behavior, f, aspect_ratio)]
            inclusion_behaviors = [(inclusion.behavior, f, inclusion.aspect_ratio) for (inclusion,f) in inclusions.items()]
            # Calcul du comportement homogénéisé
            h_behavior = self.compute_behavior(Cm, inclusion_behaviors)
            h_behavior = {parameter: value.real for (parameter,value) in h_behavior.items()} # Conversion des valeurs éventuellement complexes en valeurs réelles
            h_behavior = complete_behavior(h_behavior)
        # Cas visco-élastique
        else:
            # Initialisation du résultat
            h_behavior = {}
            # Calcul du comportement en fonction de la fréquence
            for i in range(len(frequency)):
                # Récupération du comportement de la matrice à la fréquence i
                Cm = {parameter: values[i] for (parameter,values) in microstructure.behavior.items()}
                # Récupération des comportements des inclusions à la fréquence i
                inclusion_behaviors = [] # Initialisation
                for inclusion, f in inclusions.items():
                    inclusion_behavior = {parameter: values[i] for (parameter, values) in inclusion.behavior.items()}
                    inclusion_behaviors.append((inclusion_behavior, f, inclusion.aspect_ratio))
                # Calcul du comportement homogénéisé à la fréquence i
                h_behavior_i = self.compute_behavior(Cm, inclusion_behaviors)
                h_behavior_i = complete_behavior(h_behavior_i)
                # Ajout à la liste des comportements
                for parameter, value in h_behavior_i.items():
                    try:
                        h_behavior[parameter].append(value)
                    except KeyError:
                        # Création de l'entrée associée au paramètre
                        h_behavior[parameter] = [value]
        # Renvoi du résultat
        return h_behavior
                    
        
class Mori_Tanaka(Model):
    """
    Modèle de Mori-Tanaka. Contient:
    - Une fonction de description du modèle
    - Un fonction qui renvoie le comportement homogénéisé de la microstructure.
    """
    
    def __init__(self):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = set(['K', 'G', 'E', 'nu'])  # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes
        self.interphase = False # Vrai si le modèle fonctionne sur des inclusions avec interphase
        self.name = "Mori-Tanaka"
    
    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Calcule le comportement élastique homogène équivalent. 
        Renvoie un dict de comportement.
        Cm: (dict), dictionnaire du comportement de la matrice
        inclusion_behaviors(list), format [(Cf, f, aspect_ratio)] avec Cf les dictionnaires de comportement des inclusions et aspect_ratio un tuple contenant les deux valeurs de rapports de forme
        """
        # Récupération du comportement de la matrice
        Km = Cm['K']
        Gm = Cm['G']
        # Récupération du comportement de l'inclusion
        Cf, f, ratio = inclusion_behaviors[0]
        Kf = Cf['K']
        Gf = Cf['G']
        # Calcul de Gh
        denominator = 5*Gm*(3*Km+4*Gm)+6*(1-f)*(Gf-Gm)*(Km+2*Gm)
        numerator = 5*f*Gm*(Gf-Gm)*(3*Km+4*Gm)
        Gh = Gm + numerator/denominator
        # Calcul de Kh
        denominator = 3*Kf+4*Gm+3*(1-f)*(Kf-Km)
        numerator = f*(Kf-Km)*(3*Km+4*Gm)
        Kh = Km + numerator/denominator
        return {'K': Kh, 'G': Gh}    

class Eshelby_Approximation(Model):
    """
    Modèle d'Eshelby.
    """
    
    def __init__(self):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0
        self.behavior_condition = set(['K', 'G', 'E', 'nu']) # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.interphase = False
        self.name = "Eshelby"
   
    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure.
        """
        # Récupération du comportement de la matrice
        Km = Cm['K']
        Gm = Cm['G']
        # Récupération du comportement de l'inclusion
        Cf, f, ratio = inclusion_behaviors[0]
        Kf = Cf['K']
        Gf = Cf['G']
        # Calcul de Gh
        denominator = 3*Km*(3*Gm+2*Gf) + 4*Gm*(2*Gm+3*Gf)
        numerator = 5*f*Gm*(Gf-Gm)*(3*Km+4*Gm)
        Gh = Gm + numerator/denominator
        # Calcul de Kh
        denominator = 3*Kf+4*Gm
        numerator = f*(Kf-Km)*(3*Km+4*Gm)
        Kh = Km + numerator/denominator
        # Renvoi du résultat
        return {'K': Kh, 'G': Gh}

class Differential_Scheme(Model):
    """
    Modèle différentiel.
    """
    
    def __init__(self):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0
        self.behavior_condition = set(['K', 'G', 'E', 'nu']) # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes  
        self.interphase = False
        self.name = "Differential"
    
    ## Fonctions utiles au calcul du comportement homogénéisé
    
    def deriv(module, f):
        """
        Fonction qui calcule les dérivée des paramètres K et G par rapport à la fraction volumique d'inclusion. Conçue pour être appelée par la fonction odeint lors de l'intégration numérique.
        module: list, contient les valeurs réelles et imaginaires des paramètres K, G courants ainsi que Kf et Gf propres à l'inclusion.
        f: float, fraction volumique d'inclusion courante.
        """
        K1, K2, G1, G2, Kf1, Kf2, Gf1, Gf2 = module
        # Construction des paramètres complexes
        K = K1 + K2*1j
        G = G1 + G2*1j
        Kf = Kf1 + Kf2*1j
        Gf = Gf1 + Gf2*1j
        nu = (3*K-2*G)/(6*K+2*G)
        # Calcul de dK
        numerator = K-Kf
        denominator = (1-f)*(1+(Kf-K)/(K+4*G/3))
        dK = -numerator/denominator
        dK1, dK2 = dK.real, dK.imag
        # Calcul de dG
        numerator = 15*(1-nu)*(G-Gf)
        denominator = (1-f)*(7-5*nu+2*(4-5*nu)*Gf/G)
        dG = -numerator/denominator
        dG1, dG2 = dG.real, dG.imag
        
        return np.array([dK1, dK2 ,dG1, dG2] + 4*[0])
    
    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés.
        """
        # Récupération du comportement de la matrice
        # Récupération du comportement de la matrice
        Km = Cm['K']
        Gm = Cm['G']
        # Récupération du comportement de l'inclusion
        Cf, f_finale, ratio = inclusion_behaviors[0]
        Kf = Cf['K']
        Gf = Cf['G']
        # Initialisation des paramètres d'intégration
        npoints = 100 # Nombre de points d'intégration
        f = np.linspace(0, f_finale, npoints) # Liste des fractions volumiques diluées
        initial_module = []
        for parameter in [Km, Gm, Kf, Gf]:
            initial_module += [parameter.real, parameter.imag]
        initial_module = np.array(initial_module)
        # Intégration de l'équation différentielle
        module = odeint(Differential_Scheme.deriv, initial_module, f)
        # Récupération du résultat final
        final_module = module[-1]
        Kh1, Kh2, Gh1, Gh2 = final_module[:4]  
        # Renvoi du résultat
        return {'K': Kh1+1j*Kh2, 'G': Gh1+1j*Gh2}

class Autocoherent_Hill(Model):
    """
    Modèle autocohérent de Hill.
    """
    def __init__(self):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = set(['K', 'G','E', 'nu'])  # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.interphase = False # Vrai si le modèle fonctionne sur des inclusions avec interphase
        self.name = "Self-consistent"
        self.precision = 10**-12
        self.n_point_fixe = 100
    
    def Reccurence(Module,f):
        K,G,Km,Gm,Kf,Gf = Module
        ##Calcul de Kn+1
        numerator = f*(Kf-Km)*(3*K+4*G)
        denominator = 3*Kf+4*G
        nextK = Km + numerator/denominator
        ##Calcul de Gn+1
        numerator = 5*f*G*(Gf-Gm)*(3*K+4*G)
        denominator = 3*K*(3*G+2*Gf)+4*G*(3*Gf+2*G)        
        nextG = Gm + numerator/denominator
        return nextK,nextG
    
  
    def compute_behavior(self, Cm, inclusion_behaviors):
        # Récupération du comportement de la matrice
        Km = Cm['K']
        Gm = Cm['G']
        # Récupération du comportement de l'inclusion
        Cf, f, ratio = inclusion_behaviors[0]
        Kf = Cf['K']
        Gf = Cf['G']
        F = np.linspace(0,f,self.n_point_fixe)
        
        Kinit = Km
        Ginit = Gm
        for i in range(len(F)) : 
            fi = F[i]
            # Initialisation du point fixe
            K = Kinit
            G = Ginit
            # Algorithme du point fixe
            precision = self.precision
            nextK,nextG=Autocoherent_Hill.Reccurence([K,G,Km,Gm,Kf,Gf],fi)
            while abs((nextK-K)/K) > precision or abs((nextG-G)/G) > precision : 
                K,G=nextK,nextG
                nextK,NextG=Autocoherent_Hill.Reccurence([K,G,Km,Gm,Kf,Gf],fi)  
            # Mise à jour de l'initialisation
            Kinit = nextK
            Ginit = nextG
        return {'K': nextK, 'G': nextG}
    
class Autocoherent_III(Model):
    """
    Hypothèses : 
    -isotrope
    -renforts sphériques (TO DO : PASSER AUX ELLIPSOÏDES)
    -déformations elastiques ??? (A vérifier)
    TODO : 
    -déterminer précisément les toutes les microstructures admises
    Modèle des autocohérent. Contient :
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
        self.behavior_condition = set(['K', 'G','E', 'nu'])  # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.interphase = False # Vrai si le modèle fonctionne sur des inclusions avec interphase
        self.name = "Generalised self-consistent"
  
    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés.
        """
        # Récupération du comportement de la matrice
        Km, Gm, num = Cm['K'], Cm['G'], Cm['nu']
        # Récupération du comportement de l'inclusion
        Cf, f, ratio = inclusion_behaviors[0]
        Kf, Gf, nuf = Cf['K'], Cf['G'], Cf['nu']

        ##Quelques constantes utiles au calcul de G         
        dm=(Gf/Gm)-1 
        eta1=dm*(49-50*nuf*num)+35*(dm+1)*(nuf-2*num)+35*(2*nuf-num) 
        eta2=5*nuf*(dm-7)+7*(dm+5) 
        eta3=(dm+1)*(8-10*num)+(7-5*num) 
        
        A=8*f**(10/3)*eta1*dm*(4-5*num)-2*f**(7/3)*(63*dm*eta2+2*eta1*eta3)+252*dm*eta2*f**(5/3)-50*dm*(7-12*num+8*num**2)*eta2*f+4*(7-10*num)*eta2*eta3 
        B=-4*dm*(1-5*num)*eta1*f**(10/3)+4*(63*dm*eta2+2*eta1*eta3)*f**(7/3)-504*dm*eta2*f**(5/3)+150*dm*(3-num)*num*eta2*f+3*(15*num-7)*eta2*eta3 
        D=4*dm*(5*num-7)*eta1*f**(10/3)-2*f**(7/3)*(63*dm*eta2+2*eta1*eta3)+252*dm*eta2*f**(5/3)+25*dm*(num**2-7)*eta2*f-(7+5*num)*eta2*eta3 
        
        ## Calcul de G
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
    Hypothèses : 
    -isotrope
    -renforts sphériques (TO DO : PASSER AUX ELLIPSOÏDES)
    -déformations elastiques ??? (A vérifier)
    TODO : 
    -déterminer précisément les toutes les microstructures admises
    Modèle des autocohérent. Contient :
    - Une fonction qui vérifie si le modèle est appliquable à une microstructure.
    - Une fonction de description du modèle (TODO : écrire une fonction qui renvoie une description du modèle sous forme de str et qui pourrait être appelée dans le main)
    - Un fonction qui renvoie le comportement homogénéisé de la microstructure.
    - Des fonctions qui calculent une caractéristique particulière (fraction volumique d'une inclusion, rayon d'une inclusion, comportement d'une inclusion, etc..) à partir d'un comportement homogénéisé cible (TODO)
    """
    def __init__(self, R_inclusion=1):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = set(['K', 'G','E', 'nu'])  # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.interphase = True # Vrai si le modèle fonctionne sur des inclusions avec interphase
        self.R_inclusion = R_inclusion
        self.name = "4-phases self-consistent"

    def compute_behavior(self, Cm, inclusion_behaviors):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        # Récupération des paramètres matériaux
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
    
     
 
#%% Fonctions utiles
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

def bulk_to_shear(K, E):
    """
    Transforme les modules K et E en modules G et nu
    """
    G = 3*K*E/(9*K-E)
    nu = (3*K-E)/(6*K)
    return G, nu

def complete_behavior(behavior):
    """
    Si le comportement en entrée est isotrope, le complète avec E et nu ou K et G.
    Sinon, le complète avec C ou S si la matrice entrée est inversible.
    """
    parameters = list(behavior.keys())
    result = behavior
    nu_max = 0.4999999
    # Cas de paramètres nuls (milieux poreux) et cas incompressible (nu=0.5)
    for parameter, values in result.items():
        # Comportements isotropes élastiques
        if type(values) in [float, int]:
            if values==0:
                result[parameter] = 10**(-12)
            elif values==0.5 and parameter=='nu':
                result[parameter] = nu_max
        # Comportements isotropes visco-élastiques
        elif type(values)==np.ndarray and type(values[0])!=np.ndarray:
            for i, value in enumerate(values):
                if value==0:
                    result[parameter][i] = 10**(-12)
                elif value==0.5 and parameter=='nu':
                    result[parameter][i] = nu_max
    # Isotrope K et G
    if parameters[:2]==['K', 'G'] or parameters[:2]==['G', 'K']:
        K, G = behavior['K'], behavior['G']
        E, nu = bulk_to_young(K, G)
        result['E'], result['nu'] = E, nu
    # Isotrope E et nu
    elif parameters[:2]==['E', 'nu'] or parameters[:2]==['nu', 'E']:
        E, nu = behavior['E'], behavior['nu']        
        K, G = young_to_bulk(E, nu)
        result['K'], result['G'] = K, G
    # Isotrope K et E
    elif parameters[:2]==['K', 'E'] or parameters[:2]==['E', 'K']:
        K, E = behavior['K'], behavior['E']        
        G, nu = bulk_to_shear(K, E)
        result['G'], result['nu'] = G, nu
    # Anisotrope
    elif parameters[0]=='C':
        C = behavior['C']
        try:
            S = np.linalg.inv(C)
        except:
            # C non inversible
            S = None
        result['S'] = S
    elif parameters[0]=='S':
        S = behavior['S']
        try:
            C = np.linalg.inv(S)
        except:
            # C non inversible
            C = None
        result['C'] = C
    
    # Renvoi du résultat
    return result


#%% Définition des modèles, comportements et géométries d'inclusions 
list_models = [Mori_Tanaka, Eshelby_Approximation, Differential_Scheme, Autocoherent_Hill, Autocoherent_III, Autocoherent_IV]
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
    
    
    
    
    
    
