"""
Homogeneisation - classes.py

Définition des classes utilisées.

Les items à modifier/améliorer/ajouter sont marqués en commentaires précédés de la mention "TODO".

Authors : Karim AÏT AMMAR, Enguerrand LUCAS

27/04/2020
"""

import numpy as np
import numpy.polynomial.polynomial as nppol

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
        f_phase contient la valeur de la fraction volumique de cette phase
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
        return "Phase : {}, Position : {}, Rayon : {},".format(self.type_to_str(), self.position, self.radius)
    
    def __repr__(self):
        return str(self)
    

class Inclusion_multiphase:
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, n_phases_inclusion, list_phases,behavior_condition):
        """
        type_inclusion : (int), 10 pour des inclusions sphériques multiphases.
        n_phases_inclusion : (int) nombre de phase du modèle= nombre de phase de l'inclusion +1
        list_phases : liste qui contient les différentes phases de l'inclusion de la plus interne à la plus externe, puis la phase de la matrice : [ phase_1,   phase_2, phase_matrice] 
        f_inclusion : fraction volumique de l'inclusion : égale à (R1/R2)^3
        Remarques :
            Dans le cas autocohérent à 2 phases, le rayon de la phase matrice permet de définir la fraction volumique
            Dans le cas à 3 phases, c'est le rayon de l'interphase qui joue ce rôle, mais son comportement n'a pas à être défini ??
        """
        self.type_inclusion = type_inclusion
        self.n_phases_inclusion=n_phases_inclusion
        self.list_phases = list_phases
        self.behavior_condition=behavior_condition
        self.f_inclusion = self.compute_f_inclusion()
        
    
    def compute_f_inclusion(self) : 
        if len(self.list_phases) <= 1 :
            raise NameError("There is less than 2 phase : there must be at least one inclusion phase and one matrix phase")
            return False
        else : 
            R1 = self.list_phases[0].radius
            R2 = self.list_phases[1].radius
            return (R1/R2)**3
    
    
        
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
        for i in range(self.n_phases_inclusion):
            phase=self.list_phases[i]
            s+="Phase : {}, Rayon : {} \n".format(phase.position, phase.radius)
        return s

    def __repr__(self):
        return str(self)
    
    def check_phases_inclusion(self):
        # vérification du nombre de phase
        if self.n_phases_inclusion+1 != len(self.list_phases):
            raise NameError("Error on phase number")
            return False
        if self.n_phases_inclusion < 1 : 
            raise NameError("There is less than one phase in the inclusion")
            return False
        last_position = 0
        last_radius = 0
        # vérification de la fraction volumique
        for i in range(len(self.list_phases)):
            phase=self.list_phases[i]
            # vérification du type de phase (géométrie)
            if self.type_inclusion!=phase.type_phase:
                raise NameError("Phases and inclusion have different geometry")
                return False
            # vérification du comportement des phases
            behavior = phase.behavior
            if list(behavior.keys()) != list(self.behavior_condition):
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
        self.n_phases_modèle=self.inclusion.n_phases_inclusion+1
        
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
        if inclusion.n_phases_inclusion+1 != self.n_phases_modèle :
            raise NameError("Inclusion number of phases does not match model hypothesis")
            return False
        # Vérification de la bonne définition de l'inclusion
        if not inclusion.check_phases_inclusion():
            return False
        
        # A ce stade, tout est bon
        return True

######## FONCTIONS UTILISEES POUR LE CALCUL DE Kh ################

    def J(r,phase):
        k=phase.behavior["K"]
        g=phase.behavior["G"]
        J=np.matrix([[r,1/r**2],[3*k,-4*g/(r**3)]])
        return J
    
    def N(k,inclusion):
        list_phases=inclusion.list_phases
        Rk=list_phases[k].radius
        return np.matmul(Autocohérent.J(Rk,list_phases[k+1]).I,Autocohérent.J(Rk,list_phases[k]))
    
    def Q(k,inclusion):
        Q=Autocohérent.N(0,inclusion)
        for i in range(1,k):
            Q=np.matmul(Q,Autocohérent.N(i,inclusion))
        return Q
    
    def compute_Kh(self):
        n_phases_inclusion=self.n_phases_modèle-1
        inclusion=self.inclusion
        Q=Autocohérent.Q(n_phases_inclusion-1, inclusion)
        a=Q[0,0]
        b=Q[1,0]
        last_phase=inclusion.list_phases[self.n_phases_modèle-1]
        kn=last_phase.behavior["K"]
        gn=last_phase.behavior["G"]
        rn=last_phase.radius
        numerator = 3*kn*rn**3*a-4*gn*b
        denominator = 3*(rn**3*a+b)
        return numerator/denominator


######## FONCTIONS UTILISEES POUR LE CALCUL DE Gh ################  

    def L(r,phase):
        k=phase.behavior["K"]
        g=phase.behavior["G"]
        mu=(3*k-2*g)/(6*k+2*g)
        Line1=[r , -6*mu*r**3/(1-2*mu), 3/r**4, (5-4*mu)/(r**2*(1-2*mu))]
        Line2=[r , (4*mu-7)*r**3/(1-2*mu) , -2/r**4 , 2/r**2]
        Line3=[g , 3*mu*g*r**2/(1-2*mu) , -12*g/r**5 , 2*(mu-5)*g/((1-2*mu)*r**3)]
        Line4=[g , -(7+2*mu)*g*r**2/(1-2*mu) , 8*g/r**5 , 2*(1+mu)*g/((1-2*mu)*r**3)]
        L=np.matrix([Line1,Line2,Line3,Line4])
        return L
    
    def M(k,inclusion):
        list_phases=inclusion.list_phases
        Rk=list_phases[k].radius
        return np.matmul(Autocohérent.L(Rk,list_phases[k+1]).I,Autocohérent.L(Rk,list_phases[k]))
    
    def P(k,inclusion):
        P=Autocohérent.M(0,inclusion)
        for i in range(1,k):
            P=np.matmul(P,Autocohérent.P(i,inclusion))
        return P
    
    def Z(a,b,inclusion):
        n_phases_inclusion=inclusion.n_phases_inclusion
        P1=Autocohérent.P(n_phases_inclusion-1,inclusion)
        return P1[a-1,0]*P1[b-1,1]-P1[b-1,0]*P1[a-1,1]
    
    def A(inclusion):
        n_phases_inclusion=inclusion.n_phases_inclusion
        last_phase=inclusion.list_phases[-1]
        k=last_phase.behavior["K"]
        g=last_phase.behavior["G"]
        mu=(3*k-2*g)/(6*k+2*g)
        r=last_phase.radius
        Z=Autocohérent.Z
        terme1=4*r**10*(1-2*mu)*(7-10*mu)*Z(1,2,inclusion)
        terme2=20*r**7*(7-12*mu+8*mu**2)*Z(4,2,inclusion)
        terme3=12*r**5*(1-2*mu)*(Z(1,4,inclusion)-7*Z(2,3,inclusion))
        terme4=20*r**3*(1-2*mu)**2*Z(1,3,inclusion)
        terme5=16*(4-5*mu)*(1-2*mu)*Z(4,3,inclusion)
        return terme1+terme2+terme3+terme4+terme5
    
    def B(inclusion):
        n_phases_inclusion=inclusion.n_phases_inclusion
        last_phase=inclusion.list_phases[-1]
        k=last_phase.behavior["K"]
        g=last_phase.behavior["G"]
        mu=(3*k-2*g)/(6*k+2*g)
        r=last_phase.radius
        Z=Autocohérent.Z
        terme1=3*r**10*(1-2*mu)*(15*mu-7)*Z(1,2,inclusion)
        terme2=60*r**7*mu*(mu-3)*Z(4,2,inclusion)
        terme3=-24*r**5*(1-2*mu)*(Z(1,4,inclusion)-7*Z(2,3,inclusion))
        terme4=-40*r**3*(1-2*mu)**2*Z(1,3,inclusion)
        terme5=-8*(1-5*mu)*(1-2*mu)*Z(4,3,inclusion)
        return terme1+terme2+terme3+terme4+terme5
    
    def C(inclusion):
        n_phases_inclusion=inclusion.n_phases_inclusion
        last_phase=inclusion.list_phases[-1]
        k=last_phase.behavior["K"]
        g=last_phase.behavior["G"]
        mu=(3*k-2*g)/(6*k+2*g)
        r=last_phase.radius
        Z=Autocohérent.Z
        terme1=-r**10*(1-2*mu)*(5*mu+7)*Z(1,2,inclusion)
        terme2=10*r**7*(7-mu**2)*Z(4,2,inclusion)
        terme3=12*r**5*(1-2*mu)*(Z(1,4,inclusion)-7*Z(2,3,inclusion))
        terme4=20*r**3*(1-2*mu)**2*Z(1,3,inclusion)
        terme5=-8*(7-5*mu)*(1-2*mu)*Z(4,3,inclusion)
        return terme1+terme2+terme3+terme4+terme5


    def compute_Gh(self):
        n_phases_inclusion=self.n_phases_modèle-1
        inclusion=self.inclusion
        Pi=Autocohérent.P(n_phases_inclusion-1, inclusion)
        A=Autocohérent.A(inclusion)
        B=Autocohérent.B(inclusion)
        C=Autocohérent.C(inclusion)
        last_phase=inclusion.list_phases[self.n_phases_modèle-1]
        gn=last_phase.behavior["G"]
        return gn*max(nppol.polyroots([C, B, A]))
  
    def compute_h_behavior(self):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis()
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        Gh = self.compute_Gh()
        Kh = self.compute_Kh()
        return {'K' : Kh, 'G' : Gh}

#list_models = [Autocoherent()] # Liste des modèles implémentés, penser à l'incrémenter à chaque ajout d'un nouveau modèle    

               
               
class Autocohérent_II:
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
        self.n_phases_modèle=2
        
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
        if inclusion.n_phases_inclusion + 1 != self.n_phases_modèle :
            raise NameError("Inclusion number of phases does not match model hypothesis")
            return False
        # Vérification de la bonne définition de l'inclusion
        if not inclusion.check_phases_inclusion():
            return False
        
        # A ce stade, tout est bon
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
  
    def compute_h_behavior(self,précision):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis()
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        inclusion=self.inclusion
        phase_inclusion,phase_matrix = inclusion.list_phases
        Ci, Cm = phase_inclusion.behavior, phase_matrix.behavior 
        Km,Gm = Cm['K'], Cm['G']
        Ki,Gi = Ci['K'], Ci['G']
        f = inclusion.f_inclusion
        # calcul par reccurence des modules
        K,G = Km , Gm
        nextK,nextG=Autocohérent_II.Reccurence([K,G,Km,Gm,Ki,Gi],f)
        while nextK-K > précision or nextG-G > précision : 
            K,G=nextK,nextG
            nextK,NextG=Autocohérent_II.Reccurence([K,G,Km,Gm,Ki,Gi],f)            
        return {'K' : nextK, 'G' : nextG}
    
    
class Autocohérent_III:
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
        self.n_phases_modèle=3
        
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
        if inclusion.n_phases_inclusion + 1 != self.n_phases_modèle :
            raise NameError("Inclusion number of phases does not match model hypothesis")
            return False
        # Vérification de la bonne définition de l'inclusion
        if not inclusion.check_phases_inclusion():
            return False
        
        # A ce stade, tout est bon
        return True

    
    def compute_Kh(self):
        inclusion=self.inclusion
        phase_inclusion,interphase,phase_matrix = inclusion.list_phases
        Ci, Cm = phase_inclusion.behavior, phase_matrix.behavior 
        Km,Gm = Cm['K'], Cm['G']
        Ki,Gi = Ci['K'], Ci['G']
        f = inclusion.f_inclusion
        denominator=1/(Ki-Km)+3*(1-f)/(3*Km+4*Gm)
        return Km+f/denominator
    
    def compute_Gh(self):
        inclusion=self.inclusion
        phase_inclusion,interphase,phase_matrix = inclusion.list_phases
        Ci, Cm = phase_inclusion.behavior, phase_matrix.behavior 
        Km,Gm = Cm['K'], Cm['G']
        mum=(3*Km-2*Gm)/(6*Km+2*Gm)
        Ki,Gi = Ci['K'], Ci['G']
        mui=(3*Ki-2*Gi)/(6*Ki+2*Gi)        
        f = inclusion.f_inclusion
        #constantes utiles
        g = Gi/Gm-1
        e1 = (49+35*mui-70*mum-50*mui*mum)*g+105*(mui-mum)
        e2 = (7+5*mui)*g+35*(1-mui)
        e3 = 2*(4-5*mum)*g+15*(1-mum)
        D = 2*(63*g*e2+2*e1*e3)*f**(7/3)-252*g*e2*f**(5/3)
        # termes de l'equation du second degré suivant JDhomogeneisation
        A = 8*(5*mum-4)*g*e1*f**(10/3)+D+50*(8*mum**2-12*mum+7)*g*e2*f+4*(10*mum-7)*e2*e3
        B = 4*(5*mum-1)*g*e1*f**(10/3)+2*D-150*(mum-3)*mum*g*e2*f+3*(15*mum-7)*e2*e3
        C = -4*(5*mum-7)*g*e1*f**(10/3)+D-25*(mum**2-7)*g*e2*f+(5*mum+7)*e2*e3
       
      

##Calcul séparé des termes suivant Christiensen
#
#        g = Gi/Gm
#        e1 = (g-1)*(49-50*mui*mum)+35*g*(mui-2*mum)+35*(2*mui-mum)
#        e2 = 5*mui*(g-8)+7*(Gi+Gm+4)
#        e3 = 2*(4-5*mum)*g+(7-5*mum)
#        
#        A1 = 8*(g-1)*(4-5*mum)*e1*f**(10/3)
#        A2 = -2*(63*(g-1)*e2+2*e1*e3)*f**(7/3)
#        A3 = 252*(g-1)*e2*f**(5/3)
#        A4 = -25*(g-1)*(7-12*mum+8*mum**2)*e2*f
#        A5 = 4*(7-10*mum)*e2*e3
#        A = A1+A2+A3+A4+A5
#        
#        B1 = -4*(g-1)*(1-5*mum)*e1*f**(10/3)
#        B2 =  4*(63*(g-1)*e2+2*e1*e3)*f**(7/3)
#        B3 = -504*(g-1)*e2*f**(5/3)
#        B4 = 150*(g-1)*(3-mum)*mum*e2*f
#        B5 = 3*(15*mum-7)*e2*e3
#        B = B1+B2+B3+B4+B5
#        
#        C1 = 4*(g-1)*(5*mum-7)*e1*f**(10/3)
#        C2 = -2*(63*(g-1)*e2+2*e1*e3)*f**(7/3)
#        C3 = 252*(g-1)*e2*f**(5/3)
#        C4 = 25*(g-1)*(mum**2-7)*e2*f
#        C5 = -(7-5*mum)*e2*e3
#        C = C1+C2+C3+C4+C5
#        
#        print(A,B,C,Gm*max(nppol.polyroots([C, -B, A])))
        
        return Gm*max(nppol.polyroots([C, -B, A]))
        
  
    def compute_h_behavior(self):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis()
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        Kh = self.compute_Kh()
        Gh = self.compute_Gh()
        return {'K' : Kh, 'G' : Gh}
#    
    
##### Tests
#phase1=Phase(10,1,5,{"K":200, "G":450})
#phase2=Phase(10,2,10,{"K":200, "G":450})
#matrix3=Phase(10,3,15,{"K":100, "G":150})
#
#
#
#inclusion2=Inclusion_multiphase(10,2,[phase1, phase2, matrix3],['G', 'K'])
#
#modele1=Autocohérent(inclusion1)
#modele3=Autocohérent(inclusion2)
#modele4=Autocohérent_III(inclusion2)
#
#
##print(inclusion1.check_phases_inclusion())
##print(inclusion2.check_phases_inclusion())
##print(inclusion3.check_phases_inclusion())
#
##print(modele1.check_hypothesis())
##print(modele2.check_hypothesis())
##print(modele3.check_hypothesis())
##print(modele1.compute_h_behavior())
#
#
#print(modele3.compute_h_behavior())
##print(modele4.compute_h_behavior())
#
##print(modele3.compute_h_behavior())
##print(modele4.compute_h_behavior(0.01))
#
#
phase4=Phase(10,1,5,{"K":100, "G":150})
matrix5=Phase(10,2,10,{"K":100, "G":150})
inclusion1=Inclusion_multiphase(10,1,[phase4, matrix5],['G', 'K'])
modele2=Autocohérent_II(inclusion1)
print(modele2.compute_h_behavior(0.1))

phase4=Phase(10,1,5,{"K":150, "G":150})
matrix5=Phase(10,2,10,{"K":100, "G":150})
inclusion1=Inclusion_multiphase(10,1,[phase4, matrix5],['G', 'K'])
modele2=Autocohérent_II(inclusion1)
print(modele2.compute_h_behavior(0.1))

phase4=Phase(10,1,5,{"K":200, "G":150})
matrix5=Phase(10,2,10,{"K":100, "G":150})
inclusion1=Inclusion_multiphase(10,1,[phase4, matrix5],['G', 'K'])
modele2=Autocohérent_II(inclusion1)
print(modele2.compute_h_behavior(0.1))

phase4=Phase(10,1,5,{"K":200, "G":250})
matrix5=Phase(10,2,10,{"K":100, "G":150})
inclusion1=Inclusion_multiphase(10,1,[phase4, matrix5],['G', 'K'])
modele2=Autocohérent_II(inclusion1)
print(modele2.compute_h_behavior(0.1))

phase4=Phase(10,1,5,{"K":200, "G":350})
matrix5=Phase(10,2,10,{"K":100, "G":150})
inclusion1=Inclusion_multiphase(10,1,[phase4, matrix5],['G', 'K'])
modele2=Autocohérent_II(inclusion1)
print(modele2.compute_h_behavior(0.1))    

phase4=Phase(10,1,5,{"K":200, "G":450})
matrix5=Phase(10,2,10,{"K":100, "G":150})
inclusion1=Inclusion_multiphase(10,1,[phase4, matrix5],['G', 'K'])
modele2=Autocohérent_II(inclusion1)
print(modele2.compute_h_behavior(0.1))    

phase4=Phase(10,1,5,{"K":200, "G":450})
matrix5=Phase(10,2,15,{"K":100, "G":150})
inclusion1=Inclusion_multiphase(10,1,[phase4, matrix5],['G', 'K'])
modele2=Autocohérent_II(inclusion1)
print(modele2.compute_h_behavior(0.1))    

phase4=Phase(10,1,5,{"K":200, "G":450})
matrix5=Phase(10,2,20,{"K":100, "G":150})
inclusion1=Inclusion_multiphase(10,1,[phase4, matrix5],['G', 'K'])
modele2=Autocohérent_II(inclusion1)
print(modele2.compute_h_behavior(0.1))   