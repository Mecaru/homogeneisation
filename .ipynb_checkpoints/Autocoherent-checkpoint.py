"""
Homogeneisation - classes.py

Définition des classes utilisées.

Les items à modifier/améliorer/ajouter sont marqués en commentaires précédés de la mention "TODO".

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

27/04/2020
"""

import numpy as np
from math import *

class Phase:
    """
    Contient les informations propres à une phase (Geometrie, Rayon, comportement).
    """
    
    def __init__(self, type_phase, position, radius, behavior):
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
        if type_inclusion == 10:
            return "spherique"
    
    def __str__(self):
        """
        Présentation de l'instance.
        """
        str_type_inclusion = self.type_to_str()
        return "Phase : {}, Position : {}, Rayon : {},".format(self.type_to_str(), self.position, self.radius)
    
    def __repr__(self):
        return str(self)
    

class Inclusion_multiphase:
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, n_phases, list_phases,behavior_condition):
        """
        type_inclusion : (int), 10 pour des inclusions sphériques multiphases.
        n_phases : (int) nombre de phase de l'inclusion
        phases : dict qui contient les différentes phases de l'inclusion : [1 : phase_1,  2 : phase_2] 
        """
        self.type_inclusion = type_inclusion
        self.n_phases=n_phases
        self.list_phases = list_phases
        self.behavior_condition=behavior_condition
    
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
        s= "Inclusion : {}, Nombre de phases : {} \n".format(str_type_inclusion, str(self.n_phases))
        for i in range(self.n_phases):
            phase=self.list_phases[i]
            s+="Phase : {}, Rayon : {} \n".format(phase.position, phase.radius)
        return s

    def __repr__(self):
        return str(self)
    
    def check_phases_inclusion(self):
        # vérification du nombre de phase
        if self.n_phases != len(self.list_phases):
            raise NameError("Error on phase number")
            return False
        last_position = 0
        last_radius = 0
        for i in range(self.n_phases):
            phase=self.list_phases[i]
            # vérification du type de phase (géométrie)
            if self.type_inclusion!=phase.type_phase:
                raise NameError("Phases and inclusion have different geometry")
                return False
            # vérification du comportement des phases
            behavior = phase.behavior
            if list(behavior.keys()) != self.behavior_condition:
                raise NameError("Phases and inclusion have different behavior type")
                return False
        # vérification des positions et rayons des phases
            if phase.position != last_position+1  :
                raise NameError("Error on position of phases")
                return False
            if phase.radius<last_radius :
                raise NameError("Error on radius of phases")
                return False
            else :
                last_position = phase.position
                last_radius = phase.radius
        return True


    
#class Microstructure:
#    """
#    Contient des informations sur la microstructure (comportement de la matrice, inclusions, etc..). TODO : à modifier pour prendre en compte la présence ou non d'une interphase, et d'autres paramètres de modèles plus avancés.
#    """
#    
#    def __init__(self, matrix_behavior, dict_inclusions=dict()):
#        """
#        list_inclusions : (dict), sous la forme [inclusion: f_i] avec inclusion une instance de classe Inclusion et f_i la fraction volumique de ce type d'inclusion.
#        matrix_behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). TODO :  À modifier pour représenter des comportements non isotropes.
#        """
#        self.dict_inclusions = dict_inclusions
#        self.matrix_behavior = matrix_behavior
#        # Calcul de la fraction volumique de matrice f_m
#        self.f_matrix = self.compute_fm()
#        
#    def __str__(self):
#        string = "Microstructure\nf_m = {:.2f}, matrix".format(self.f_matrix)
#        dict_inclusions = self.dict_inclusions
#        # Présentation de toutes les inclusions contenues dans la microstructure
#        for inclusion in dict_inclusions.keys():
#            fi = dict_inclusions[inclusion]
#            string += "\nf_i = {}, ".format(fi) + str(inclusion)
#        return string
#
#    def compute_fm(self):
#        """
#        1/ Vérifie si la liste des inclusions donnée est cohérente (i.e : la somme des fractions volumiques des inclusions est inférieure à 1). Si ce n'est pas le cas, génère une erreur.
#        2/ Si aucune erreur n'est générée, calcule la fraction volumique de matrice.
#        """
#        total_fi = 0 # Total des fractions volumiques d'inclusions
#        dict_inclusions = self.dict_inclusions
#        for inclusion in dict_inclusions.keys():
#            fi = dict_inclusions[inclusion]
#            total_fi += fi
#        if total_fi >= 1:
#            raise NameError("Inconsistent list of volumic fractions")
#        else :
#            f_m = 1 - total_fi
#            return f_m



        
class Autocohérent:
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
    def __init__(self,inclusion_multiphase):
        """
        Définition des hypothèses du modèle.
        """
        self.inclusion=inclusion_multiphase
        self.type_inclusion = 10
        self.behavior_condition = ["K", "G"] # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        self.n_phases=self.inclusion.n_phases
        
    def __str__(self):
        """
        Description textuelle du modèle.
        """
        return "Modèle autocohérent à {} phases".format(self.n_phases)
    
    def check_hypothesis(self):
        """
        Vérifies si la microstructure vérifie les hypothèses du modèle, renvoie un booléens. 
        TODO : Éventuellement généraliser cette fonction en l'incorporant dans une classe mère Model pour qu'elle s'applique à tous les modèles.
        """
        inclusion = self.inclusion
        # Vérification du type d'inclusion
        if inclusion.type_inclusion != self.type_inclusion:
            raise NameError("Inclusion geometry does not match model hypothesis")
            return False
        # Vérification du nombre de phase de l'inclusion
        if inclusion.n_phases != self.n_phases :
            raise NameError("Inclusion number of phases does not match model hypothesis")
            return False
        # Vérification de la bonne définition de l'inclusion
        if not inclusion.check_phases_inclusion():
            return False
        
        # A ce stade, tout est bon
        return True

######## FONCTIONS UTILISEES POUR LE CALCUL DE Kh ################

    def J(k,r,phase):
        k=phase.behavior["K"]
        g=phase.behavior["G"]
        J=np.matrix([[r,1/r**2],[3*k,-4*g/(r**4)]])
        return J
    
    def N(k,inclusion):
        list_phases=inclusion.list_phases
        return np.matmul(Autocohérent.J(k+1,list_phases[k].radius,list_phases[k+1]).I,Autocohérent.J(k,list_phases[k].radius,list_phases[k]))
    
    def Q(k,inclusion):
        Q=Autocohérent.N(0,inclusion)
        for i in range(1,k):
            Q=np.matmul(Q,Autocohérent.N(i,inclusion))
        return Q
    
    def compute_Kh(self):
        n_phases=self.n_phases
        inclusion=self.inclusion
        Qi=Autocohérent.Q(n_phases-1, inclusion)
        a=Qi[0,0]
        b=Qi[1,0]
        last_phase=inclusion.list_phases[n_phases-1]
        kn=last_phase.behavior["K"]
        gn=last_phase.behavior["G"]
        rn=last_phase.radius
        numerator = 3*kn*rn**3*a-4*gn*b
        denominator = 3*(rn**3*a+b)
        return numerator/denominator
    
    
    
    def compute_h_behavior(self):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis(microstructure)
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        Gh=compute_Gh(self)
        Kh=compute_Kh(self)
        return [Gh,Kh]

#list_models = [Autocoherent()] # Liste des modèles implémentés, penser à l'incrémenter à chaque ajout d'un nouveau modèle    
    
# Tests
phase1=Phase(10,1,5,{"K":100, "G":150})
phase2=Phase(10,2,10,{"K":100, "G":150})
phase3=Phase(10,3,15,{"K":100, "G":150}) ####TRUC TRES LOUCHE SI ON MET K : 100,100,400,100
phase4=Phase(10,4,15,{"K":100, "G":150})

inclusion1=Inclusion_multiphase(10,2,[phase1, phase2],['K', 'G'])
inclusion2=Inclusion_multiphase(10,3,[phase1, phase2, phase3],['K', 'G'])
inclusion3=Inclusion_multiphase(10,4,[phase1, phase2, phase3, phase4],['K', 'G'])

modele1=Autocohérent(inclusion1)
modele2=Autocohérent(inclusion2)
modele3=Autocohérent(inclusion3)

#print(inclusion1.check_phases_inclusion())
#print(inclusion2.check_phases_inclusion())
#print(inclusion3.check_phases_inclusion())

#print(modele1.check_hypothesis())
#print(modele2.check_hypothesis())
#print(modele3.check_hypothesis())

print(modele1.compute_Kh())
print(modele2.compute_Kh())
print(modele3.compute_Kh())
    