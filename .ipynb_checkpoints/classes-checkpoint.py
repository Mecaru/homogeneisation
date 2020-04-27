"""
Homogeneisation - classes.py

Définition des classes utilisées.

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

27/04/2020
"""

import numpy as np
from math import *

list_models = [] # Liste des modèles implémentés


class Inclusion():
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, radius, behavior):
        """
        type_inclusion : (int), 0 pour des inclusions sphériques.
        radius : (float), valeur du rayon des inclusions sphériques. TODO : À remplacer par un paramètre plus général pour des inclusions de types différents. 
        behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). TODO :  À modifier pour représenter des comportements non isotropes.
        """
        self.type_inclusion = type_inclusion
        self.radius = radius
        self.behavior = behavior
        self.K = behavior['K'] # Bulk modulus
        self.G = behavior['G'] # Shear modulus
    
    def type_to_str(self):
        """
        Transforme un entier "type_inclusion" en la chaîne de caractères correspondante (exemple : 0 --> "spheres") 
        """
        type_inclusion = self.type_inclusion
        if type_inclusion == 0:
            return "spheres"
    
    def __str__(self):
        """
        Présentation de l'instance.
        """
        str_type_inclusion = self.type_to_str()
        return "Inclusion : {}, radius : {}".format(str_type_inclusion, self.radius)

    def __repr__(self):
        return str(self)

    
class Microstructure():
    """
    Contient des informations sur la microstructure (comportement de la matrice, inclusions, etc..). TODO : à modifier pour prendre en compte la présence ou non d'une interphase, et d'autres paramètres de modèles plus avancés.
    """
    
    def __init__(self, matrix_behavior, dict_inclusions=dict()):
        """
        list_inclusions : (dict), sous la forme [inclusion: f_i] avec inclusion une instance de classe Inclusion et f_i la fraction volumique de ce type d'inclusion.
        matrix_behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). TODO :  À modifier pour représenter des comportements non isotropes.
        """
        self.dict_inclusions = dict_inclusions
        self.matrix_behavior = matrix_behavior
        # Calcul de la fraction volumique de matrice f_m
        self.f_matrix = self.compute_fm()
        
    def __str__(self):
        string = "Microstructure\nf_m = {:.2f}, matrix".format(self.f_matrix)
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
            raise NameError("Inconsistent list of volumic fractions")
        else :
            f_m = 1 - total_fi
            return f_m

# Tests
inclusion1 = Inclusion(0, 1, {"K":300, "G":150})
inclusion2 = Inclusion(0, 2, {"K":300, "G":150})
microstructure = Microstructure({"K":10, "G":15}, {inclusion1:0.6, inclusion2:0.2})
print(microstructure)