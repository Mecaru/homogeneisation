"""
Homogeneisation - classes.py

Définition des classes utilisées.

Les items à modifier/améliorer/ajouter sont marqués en commentaires précédés de la mention "TODO".

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

27/04/2020
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt


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

    
class Microstructure:
    """
    Contient des informations sur la microstructure (comportement de la matrice, inclusions, etc..). TODO : à modifier pour prendre en compte la présence ou non d'une interphase, et d'autres paramètres de modèles plus avancés.
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
        if total_fi >= 1:
            raise NameError("The total volumic fractions of the inclusions exceed 1")
        else :
            f_m = 1 - total_fi
            return f_m
        
    def draw(self):
        """
        Méthode qui permet de dessiner la microstructure. Pour le moment, fonctionne uniquement avec une seule inclusion sphérique.
        """
        inclusions = list(self.dict_inclusions.keys())
        if len(inclusions) == 1 and inclusions[0].type_inclusion == 0:
            inclusion = inclusions[0]
            fi = self.dict_inclusions[inclusion]
            # Calcul du rayon pour un VER de taille 10X10
            r = sqrt(100*fi/pi)
            x, y = [], []
            for theta in np.linspace(0,2*pi,200):
                x.append(r*cos(theta))
                y.append(r*sin(theta))
            fig, ax = plt.subplots()
            ax.axis('equal')
            plt.plot([-5,5,5,-5,-5], [-5,-5,5,5,-5])
            plt.plot(x, y)
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
        if n_inclusions > self.n_inclusions:
            # Le modèle ne peut pas traiter de microstructures avec autant d'inclusions de natures différentes
            return False
        for inclusion in dict_inclusions.keys():
            # Vérification du type d'inclusion
            if inclusion.type_inclusion != self.type_inclusion:
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
        denominator = 1 + (1-f)*(Gf-Gm)/(Gm+Gm*(9*Km+8*Gm)/(6*(Km+2*Gm)))
        numerator = f*(Gf-Gm)
        Gh = Gm + numerator/denominator
        result = {'G' : Gh}
        # TODO : une fois le modèle complet, convertir le résultat en E et nu si les comportements en entrée sont en E et nu
        return {'G' : Gh}

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
    if parameters == ['K', 'G']:
        K, G = behavior['K'], behavior['G']
        E, nu = bulk_to_young(K, G)
        result['E'], result['nu'] = E, nu
    elif parameters == ['E', 'nu']:
        E, nu = behavior['E'], behavior['nu']
        K, G = young_to_bulk(E, nu)
        result['K'], result['G'] = K, G
    return result
    
list_models = [Mori_Tanaka] # Liste des modèles implémentés, à incrémenter à chaque ajout d'un nouveau modèle
dict_behaviors = {'Isotropic (K & G)': ['K', 'G'], 'Isotropic (E & nu)': ['E', 'nu']}
dict_types = {0: 'Spheres', 1: 'Oblate', 2: 'Prolate'} # Types de géométries admissibles et leur identifiant

# Tests
#inclusion1 = Inclusion(1, {"E":300, "nu":0.3})
#print(inclusion1)
#inclusion1 = Inclusion(0, {"K":300, "G":0.3})
#print(inclusion1)
#inclusion2 = Inclusion(0, {"K":300, "G":150})
#microstructure = Microstructure({"E":10, "nu":0.1}, {inclusion1:0.6})
#model = Mori_Tanaka()
#print(microstructure)
#print(model.check_hypothesis(microstructure))
#print(model.compute_h_behavior(microstructure))
#print(microstructure)
#microstructure.draw()
