"""
Homogeneisation - classes.py

Définition des classes utilisées.

Les items à modifier/améliorer/ajouter sont marqués en commentaires précédés de la mention "TODO".

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

27/04/2020
"""

import numpy as np
from math import *


class Inclusion_multiphase():
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, radius, behavior):
        """
        type_inclusion : (int), 10 pour des inclusions sphériques multiphases.
        n_phase : (int) nombre de phase de l'inclusion
        phases : dict qui contient les différentes phases de l'inclusion : [1 : phase_1,  2 : phase_2] 
        """
        self.type_inclusion = type_inclusion
        self.n_phase=n_phase
        self.dict_phases = dict_phases
    
    def type_to_str(self):
        """
        Transforme un entier "type_inclusion" en la chaîne de caractères correspondante (exemple : 0 --> "spheres") 
        """
        type_inclusion = self.type_inclusion
        if type_inclusion == 10:
            return "spherique multiphase"
    
    def __str__(self):
        """
        Présentation de l'instance.        
        """
        str_type_inclusion = self.type_to_str()
        s= "Inclusion : {}, Nombre de phases : {}".format(str_type_inclusion, self.n_phase)
        for phase in self.dict_phases.keys():
            s+="Phase : {}, Rayon : {}".format(phase.position, phase.radius)
        return s

    def __repr__(self):
        return str(self)

    
    
class Phase
    """
    Contient les informations propres à une phase (Geometrie, Rayon, comportement).
    """
    
    def __init__(self, type_inclusion, radius, behavior):
        """
        type_inclusion : (int), 10 pour des phases sphériques.
        position : (int), précise la position qu'occupe la phase dans l'inclusion, 1 étant la phase centrale
        radius :(int).rayon de la phase
        behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus) de chaque phase. 
        """
        self.type_phase = type_phase
        self.position = position
        self.radius = radius
        self.behavior = behavior
    
    def type_to_str(self):
        """
        Transforme un entier "type_inclusion" en la chaîne de caractères correspondante (exemple : 0 --> "spheres") 
        """
        type_inclusion = self.type_phase
        if type_inclusion == 0:
            return "spherique"
    
    def __str__(self):
        """
        Présentation de l'instance.
        """
        str_type_inclusion = self.type_to_str()
        return "Phase : {}, Rayon : {}".format(self.type_to_str(), self.radius)

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



        
class Autocohérent
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
        self.type_inclusion = 10
        self.behavior_condition = ["K", "G"] # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.n_phases=2
        
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return "Modèle autocohérent à {} phases".format(self.n_phases)
    
    def check_hypothesis(self, microstructure):
        """
        Vérifies si la microstructure vérifie les hypothèses du modèle, renvoie un booléens. 
        TODO : Éventuellement généraliser cette fonction en l'incorporant dans une classe mère Model pour qu'elle s'applique à tous les modèles.
        """
        dict_inclusions = microstructure.dict_inclusions
        inclusions = dict_inclusions.keys()
        
        # vérification du nombre d'inclusions
        n_inclusions = len(inclusions)
        if n_inclusions > self.n_inclusions:
            # Le modèle ne peut pas traiter de microstructures avec autant d'inclusions de natures différentes
            return False
        for inclusion in dict_inclusions.keys():
            # Vérification du type d'inclusion
            if inclusion.type_inclusion != self.type_inclusion:
                return False
            #vérification du nombre de phase
            if inclusion.n_phase != self.n_phase :
                return False
            
            position = 1
            radius = 0
            for phase in inclusion.dict_phases:
                # vérification du comportement des phases
                behavior = phase.behavior
                if list(behavior.keys()) != self.behavior_condition:
                    return False
                # vérification des positions et rayons des phases
                if position-phase.position != 0 and phase.radius-radius>0:
                    return False
                else :
                    position = phase.position
                    radius = phase.radius
        
        # Vérification du comportement de la matrice
        if list(microstructure.matrix_behavior.keys()) != self.behavior_condition:
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
        return 0
    
    

list_models = [Autocoherent] # Liste des modèles implémentés, penser à l'incrémenter à chaque ajout d'un nouveau modèle    
    
# Tests
phase1=Phase(0,1,5,{"K":300, "G":150})
phase2=Phase(0,2,10,{"K":300, "G":150})
phase3=Phase(1,2,10,{"K":300, "G":150})

print(phase1)
#inclusion1 = Inclusion(10, 2, [5,10], [{"K":300, "G":150},{"K":300, "G":150}])
#inclusion2 = Inclusion(0, 2, {"K":300, "G":150})
#microstructure = Microstructure({"K":10, "G":15}, {inclusion1:0.6})
#model = Mori_Tanaka()
#print(model.compute_h_behavior(microstructure))
    