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
    
    def __init__(self, type_inclusion, behavior, aspect_ratio=1, name=None):
        """
        TODO : Prise en compte de l'orientation
        type_inclusion : (int), 0 pour des inclusions sphériques. Voir la liste list_types (en bas du fichier) pour les autres types.
        radius : (float), valeur du rayon des inclusions sphériques. TODO : À remplacer par un paramètre plus général pour des inclusions de types différents. 
        behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). TODO :  À modifier pour représenter des comportements non isotropes.
        """
        self.type_inclusion = type_inclusion
        self.aspect_ratio = aspect_ratio
        self.behavior = complete_behavior(behavior)
        self.name = name
    
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
            string += ", {}: {:.2f}".format(parameter, value)
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

    
class Microstructure:
    """
    Contient des informations sur la microstructure (comportement de la matrice, inclusions, etc..). TODO : à modifier pour prendre en compte la présence ou non d'une interphase, et d'autres paramètres de modèles plus avancés.
    Contient une fonction qui renvoie les bornes de Hashin-Shtrickman pour la microstructure en question 
    TODO : Généraliser ces bornes à n phases (et pas 2 comme c'est le cas ici)
    """
    
    def __init__(self, matrix_behavior, dict_inclusions=dict()):
        """
        list_inclusions : (dict), sous la forme {inclusion: f_i} avec inclusion une instance de classe Inclusion et f_i la fraction volumique de ce type d'inclusion.
        matrix_behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). TODO :  À modifier pour représenter des comportements non isotropes.
        """
        self.dict_inclusions = dict_inclusions
        self.matrix_behavior = complete_behavior(matrix_behavior)
        # Calcul de la fraction volumique de matrice f_m
        self.f_matrix = self.compute_fm()
        
    def __str__(self):
        string = "Microstructure\nf_m = {:.2f}, matrix".format(self.f_matrix, self.matrix_behavior)
        for parameter, value in self.matrix_behavior.items():
            string += ", {}: {:.2f}".format(parameter, value) # TODO : transformer cette ligne en fonction print_behavior et l'appeler lors de l'affichage du comportement homogénéisé et des bornes de hashin
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
        if total_fi > 1:
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
        numerator = c2*(k2-k1)*(4*g1+3*k1)
        denominator = (4*g1+3*k1)+3*c1*(k2-k1)
        return k1+numerator/denominator
    
    def ghs(k1, g1, c1, k2, g2, c2):
        numerator = c2*(g2-g1)
        denominator = 1+6*c1*(g2-g1)*(k1+2*g1)/((3*k1+4*g1)*5*g1)
        return g1+numerator/denominator
    
    def khsporous(k,g,c):
        numerator = 4*(1-c)*k*g
        denominator = 4*g+3*c*k
        return numerator/denominator
    
    def ghsporous(k,g,c):
        numerator = (1-c)*(8*g+9*k)*g
        denominator = 4*g*(2+3*c)+3*k*(3+2*c)
        return numerator/denominator
    
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
            
            if kf<10**-5 and gf<10**-5 :
                ksup=max(Microstructure.khsporous(km,gm,f),Microstructure.khsporous(km,gm,f))
                kinf=min(Microstructure.khsporous(km,gm,f),Microstructure.khsporous(km,gm,f))
                gsup=max(Microstructure.ghsporous(km,gm,f),Microstructure.ghsporous(km,gm,f))
                ginf=min(Microstructure.ghsporous(km,gm,f),Microstructure.ghsporous(km,gm,f))
           
            else :        
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
        self.behavior_condition = set(['K', 'G','E', 'nu'])  # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
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
            if set(behavior.keys()) != self.behavior_condition:
                #print (list(behavior.keys()) , self.behavior_condition)
                #raise NameError("Inclusion and microstructure behavior incompatible")
                return False
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()) != self.behavior_condition:
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
        
        denominator = 5*Gm*(3*Km+4*Gm)+6*(1-f)*(Gf-Gm)*(Km+2*Gm)
        numerator = 5*f*Gm*(Gf-Gm)*(3*Km+4*Gm)
        Gh = Gm + numerator/denominator
        
        denominator = 3*Km+4*Gm+3*(1-f)*(Kf-Km)
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
        self.behavior_condition = set(['K', 'G','E', 'nu']) # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
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

            if set(behavior.keys()) != self.behavior_condition:
                return False
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()) != self.behavior_condition:
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
        self.behavior_condition = set(['K', 'G','E', 'nu']) # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
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
            if set(behavior.keys()) != self.behavior_condition:
                print (list(behavior.keys()) , self.behavior_condition)
                #raise NameError("Inclusion and microstructure behavior incompatible")
                return False
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()) != self.behavior_condition:
            #raise NameError("Inclusion and microstructure behavior incompatible")
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
        numerator = c2*(k2-k1)*(4*g1+3*k1)
        denominator = (4*g1+3*k1)+3*c1*(k2-k1)
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
    
    
    
class Autocoherent_Hill:
    
    
    def __init__(self):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = set(['K', 'G','E', 'nu'])  # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.name = "Autocohérent"
        self.precision = 10**-12
        self.n_point_fixe = 100
        
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return "Modèle autocohérent"
    
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
            if set(behavior.keys()) != self.behavior_condition:
                #raise NameError("Inclusion and microstructure behavior incompatible")
                return False
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()) != self.behavior_condition:
            raise NameError("Inclusion and microstructure behavior incompatible")
            return False
        # À ce stade, toutes les conditions ont été vérifiées
        return True
    
    def Reccurence(Module,f):
        K,G,Km,Gm,Kf,Gf=Module
        ##Calcul de Kn+1
        numerator = f*(Kf-Km)*(3*K+4*G)
        denominator = 3*Kf+4*G
        nextK = Km + numerator/denominator
        ##Calcul de Gn+1
        numerator = 5*f*G*(Gf-Gm)*(3*K+4*G)
        denominator = 3*K*(3*G+2*Gf)+4*G*(3*Gf+2*G)        
        nextG = Gm + numerator/denominator
        return nextK,nextG
    
  
    def compute_h_behavior(self,microstructure):
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
            while abs(nextK-K) > precision or abs(nextG-G) > precision : 
                K,G=nextK,nextG
                nextK,NextG=Autocoherent_Hill.Reccurence([K,G,Km,Gm,Kf,Gf],fi)  
            # Mise à jour de l'initialisation
            Kinit = nextK
            Ginit = nextG
        return complete_behavior({'K' : nextK, 'G' : nextG})
    
    
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


class Autocoherent_III:
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
        self.name = "Autocohérent généralisé"
        
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return "Modèle autocohérent généralisé"
    
    def __repr__(self):
        """
        Description textuelle du modèle.
        """
        return str(self)
    
    def check_hypothesis(self, microstructure):
        """
        Vérifies si la microstructure vérifie les hypothèses du modèle, renvoie un booléens. 
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
            if set(behavior.keys()) != self.behavior_condition:
                #raise NameError("Inclusion and microstructure behavior incompatible")
                return False
            # Approximation du cas poreux
            if (behavior['K'] == 0 and behavior['G'] == 0) : 
                behavior['K'] = 10**-12
                behavior['G'] = 10**-12
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()) != self.behavior_condition:
            raise NameError("Inclusion and microstructure behavior incompatible")
            return False
        
        # À ce stade, toutes les conditions ont été vérifiées
        return True
    
   
  
    def compute_h_behavior(self,microstructure):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        """
        compatible = self.check_hypothesis(microstructure)
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
            
        dict_inclusions = microstructure.dict_inclusions
        inclusion = list(dict_inclusions.keys())[0] #Inclusion unique ici
        
        Cm = microstructure.matrix_behavior
        Cf = inclusion.behavior        
        
        Km,Gm,num = Cm['K'], Cm['G'], Cm['nu']
        Kf,Gf,nuf = Cf['K'], Cf['G'], Cf['nu']
        f = dict_inclusions[inclusion]

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
        
        return complete_behavior({'K' : Kh, 'G' : Gh})
    
    

class Autocoherent_IV:
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
    def __init__(self,R_inclusion=1):
        """
        Définition des hypothèses du modèle.
        """
        self.type_inclusion = 0 # Sphères
        self.behavior_condition = set(['K', 'G','E', 'nu'])  # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 2 # Nombre d'inclusions de natures différentes 
        self.R_inclusion = R_inclusion
        self.name = "Autocohérent 4-phases"
        
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return "Modèle autocohérent 4-phases"
    
    def __repr__(self):
        """
        Description textuelle du modèle.
        """
        return str(self)
    
    def check_hypothesis(self, microstructure):
        """
        Vérifies si la microstructure vérifie les hypothèses du modèle, renvoie un booléens. 
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
            if set(behavior.keys()) != self.behavior_condition:
                #raise NameError("Inclusion and microstructure behavior incompatible")
                return False
            # Approximation du cas poreux
            if (behavior['K'] == 0 and behavior['G'] == 0) : 
                behavior['K'] = 10**-12
                behavior['G'] = 10**-12
        # Vérification su comportement de la matrice
        if set(microstructure.matrix_behavior.keys()) != self.behavior_condition:
            raise NameError("Inclusion and microstructure behavior incompatible")
            return False
        # À ce stade, toutes les conditions ont été vérifiées
        return True

  
    def compute_h_behavior(self,microstructure):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        
        compatible = self.check_hypothesis(microstructure)
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        dict_inclusions = microstructure.dict_inclusions
        inclusion = list(dict_inclusions.keys())[0]
        interphase = list(dict_inclusions.keys())[1] ## Cela suppose qu'on a bien rentré phase et interphase dans l'ordre, avec l'épaisseur adéquate
        Cf = inclusion.behavior
        Cv = interphase.behavior
        Cm = microstructure.matrix_behavior       
        
        Km,Gm,num = Cm['K'], Cm['G'], Cm['nu']
        Kf,Gf,nuf = Cf['K'], Cf['G'], Cf['nu']
        Kv,Gv,nuv = Cv['K'], Cv['G'], Cv['nu']
        
        f = dict_inclusions[inclusion]
        cf = dict_inclusions[interphase]
  
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
    
        M1=np.zeros(shape=(4,4))
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
        
        M2=np.zeros(shape=(4,4))
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
        return complete_behavior({'K' : Kh, 'G' : Gh})


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
        if (K==0 and G==0) : 
            E, nu = 0, 0.3
        else : 
            E, nu = bulk_to_young(K, G)
        result['E'], result['nu'] = E, nu
    elif parameters[:2] == ['E', 'nu']:
        E, nu = behavior['E'], behavior['nu']
        if nu >= 0.5 : 
            nu = 0.4999999999
        K, G = young_to_bulk(E, nu)
        result['K'], result['G'] = K, G
    return result
    
list_models = [Mori_Tanaka, Eshelby_Approximation, Differential_Scheme] # Liste des modèles implémentés, à incrémenter à chaque ajout d'un nouveau modèle
dict_behaviors = {'Isotropic (K & G)': ['K', 'G'], 'Isotropic (E & nu)': ['E', 'nu']}
dict_types = {0: 'Spheres', 1: 'Oblate', 2: 'Prolate'} # Types de géométries admissibles et leur identifiant


