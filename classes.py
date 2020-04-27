"""
Homogeneisation - classes.py

Définition des classes utilisées.

Authors : Karim AÏT AMMAR, Enguerrand LUCAS
"""

import numpy as np
from math import *

list_models = [] # Liste des modèles implémentés

class Inclusion():
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, radius, f, behavior):
        """
        type_inclusion : (int), 0 pour des inclusions sphériques.
        radius : (float), valeur du rayon des inclusions sphériques. /!\ À remplacer par un paramètre plus général pour des inclusions de types différents. 
        f : (float, 0<f<1), fraction volumique de l'inclusion.
        behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). /!\ À modifier pour représenter des comportements non isotropes.
        """
        self.type_inclusion = type_inclusion
        self.radius = radius
        self.f = f
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
    


# Tests
inclusion = Inclusion(0, 1, 0.3, {"K":300, "G":150})
print(inclusion)
    