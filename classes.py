"""
Homogeneisation - classes.py

Définition des classes utilisées.

Les items à modifier/améliorer/ajouter sont marqués en commentaires précédés de la mention "TODO".

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

27/04/2020
"""

import numpy as np
from math import *


class Inclusion():
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, behavior, radius=0.1):
        """
        type_inclusion : (int), 0 pour des inclusions sphériques.
        radius : (float), valeur du rayon des inclusions sphériques. TODO : À remplacer par un paramètre plus général pour des inclusions de types différents. 
        behavior : (dict), contient les valeurs des paramètres de la matrice de comportement, pour le moment, K (bulk modulus) et G (shear modulus). TODO :  À modifier pour représenter des comportements non isotropes.
        """
        self.type_inclusion = type_inclusion
        self.radius = radius
        self.behavior = behavior
    
    def type_to_str(self):
        """
        Transforme un entier "type_inclusion" en la chaîne de caractères correspondante (exemple : 0 --> "spheres") 
        TODO : synchroniser cette fonction avec le main à l'aide d'un dictionnaire pour faciliter l'ajout de types d'inclusions
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
        self.type_inclusion = 0
        self.behavior_condition = ["G", "K"] # Le modèle s'applique sur des microstructures dont les inclusions et la matrice sont isotropes
        self.n_inclusions = 1 # Nombre d'inclusions de natures différentes 
        
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
             raise NameError("Wrong number of inclusion")
             return False
        for inclusion in dict_inclusions.keys():
            # Vérification du type d'inclusion
            if inclusion.type_inclusion != self.type_inclusion:
                raise NameError("Wrong type of inclusion or microstructure")
                return False
            # vérification du comportement des inclusions
            behavior = inclusion.behavior
            if list(behavior.keys()) != self.behavior_condition:
                print (list(behavior.keys()) , self.behavior_condition)
                raise NameError("Inclusion and microstructure behavior incompatible")
                return False
        # Vérification su comportement de la matrice
        if list(microstructure.matrix_behavior.keys()) != self.behavior_condition:
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
        denominator = 1 + (1-f)*(Gf-Gm)/(Gm+Gm*(9*Km+8*Gm)/(6*(Km+2*Gm)))
        numerator = f*(Gf-Gm)
        Gh = Gm + numerator/denominator
        return {'G' : Gh}
    
    def Hashin_bounds(microstructure):
        """
        Donne les bornes de Hashin-Shtrikman pour 1 seule phase, isotrope
        TODO : ajouter le cas des inclusion multiples
        """
        fm=microstructure.f_matrix
        f=1-fm
        km=microstructure.matrix_behavior["K"]
        gm=microstructure.matrix_behavior["G"]
        mum=(3*km-2*gm)/(6*km+2*gm)
        
        for inclusion in microstructure.dict_inclusions.keys():
            kf=inclusion.behavior["K"]
            gf=inclusion.behavior["G"]
            muf=(3*kf-2*gf)/(6*kf+2*gf)
            
        af =(3+4*muf)/(8*(1-muf))
        am =(3+4*mum)/(8*(1-mum))
        bf = (3-4*muf)/(4*(1-muf))
        bm = (3-4*mum)/(4*(1-mum))
        
        numerator = f*kf + fm*km/(1+af*((km-kf)/kf))
        denominator = f + fm/(1+af*((km-kf)/kf))
        ksup=numerator/denominator
        
        numerator = fm*km + f*kf/(1+am*((kf-km)/km))
        denominator = fm + f/(1+am*((kf-km)/km))
        kinf=numerator/denominator
        
        numerator = f*gf + fm*gm/(1+bf*((gm-gf)/gf))
        denominator = f + fm/(1+bf*((gm-gf)/gf))
        gsup=numerator/denominator
        
        numerator = fm*gm + f*gf/(1+bm*((gf-gm)/gm))
        denominator = fm + f/(1+am*((gf-gm)/gm))
        ginf=numerator/denominator
        
        return {'Gsup' : gsup, 'Ginf' : ginf, 'Ksup' : ksup, 'Kinf' : kinf}

    
    def check_bounds(self,microstructure):
        Behavior_h=self.compute_h_behavior(self,microstructure)
        Gh = Behavior_h['G']
        Kh = Behavior_h['K']
        Bounds=Hashin_bounds(microstructure)
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
    

        
     
#list_models = [Mori_Tanaka] # Liste des modèles implémentés, à incrémenter à chaque ajout d'un nouveau modèle
#dict_behaviors = {'Isotropic' : ['K', 'G'], 'Test' : ['Test']}

##Tests
inclusion1 = Inclusion(0, {"K":300, "G":150}, 1)
inclusion2 = Inclusion(0, {"K":300, "G":150}, 2)
microstructure = Microstructure({"K":10, "G":15}, {inclusion1:0.6})
model = Mori_Tanaka()
print(model.compute_h_behavior(microstructure), Mori_Tanaka.Hashin_bounds(microstructure))