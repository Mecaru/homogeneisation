"""
Homogeneisation - Autocoherent.py

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
    

class Interphase(Phase):
    
    def __init__(self, type_phase, inclusion_phase, thickness, behavior):
        Phase.__init__(self, type_phase, inclusion_phase.position+1, inclusion_phase.radius+thickness, behavior)
    
    

class Inclusion_multiphase:
    """
    Contient les informations propres à une inclusion (type, géométrie, comportement, etc...).
    """
    
    def __init__(self, type_inclusion, list_phases, matrix_behavior, f_inclusion, behavior_condition):
        """
        type_inclusion : (int), 10 pour des inclusions sphériques multiphases.
        n_phases_inclusion : (int) nombre de phase du modèle= nombre de phase de l'inclusion +1
        list_phases : liste qui contient les différentes phases de l'inclusion de la plus interne à la plus externe, puis la phase de la matrice : [ phase_1,   phase_2, phase_matrice] 
        f_inclusion : fraction volumique de l'inclusion : égale à (R1/R2)^3
        """
        self.type_inclusion = type_inclusion
        self.n_phases_inclusion=len(list_phases)
        self.list_phases = list_phases
        self.behavior_condition=behavior_condition
        self.matrix_behavior = matrix_behavior
        self.f_inclusion = f_inclusion
        
    
    def check_f_inclusion(self) : 
        if self.n_phases_inclusion <= 1 :
            return True
        else : 
            R1 = self.list_phases[0].radius
            R2 = self.list_phases[1].radius
            return (self.f_inclusion-(R1/R2)**3)<0.001
        
    def compute_f_interphase(self):
        if self.n_phases_inclusion <= 1 : 
            raise NameError ("No interphase in inclusion")
            return 0
        else : 
            R1 = self.list_phases[0].radius
            R2 = self.list_phases[1].radius
            f=self.f_inclusion
            return f*((R2/R1)**3-1)
    
    
        
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
        if self.n_phases_inclusion != len(self.list_phases):
            raise NameError("Error on phase number")
            return False
        if self.n_phases_inclusion < 1 : 
            raise NameError("There is less than one phase in the inclusion")
            return False
        last_position = 0
        last_radius = 0
        # vérification de la fraction volumique
        if not self.check_f_inclusion() :
            raise NameError ("Rayons des phases incompatibles avec les fractions volumiques")
            return True
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
        self.n_phases_modèle=2 # Nombre de phases qui doivent etre définies pour calculer le modèle
        self.précision = 0.01
        
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
  
    def compute_h_behavior(self):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis()
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        inclusion=self.inclusion
        phase_inclusion = inclusion.list_phases[0]
        Ci, Cm = phase_inclusion.behavior, inclusion.matrix_behavior 
        Km,Gm = Cm['K'], Cm['G']
        Ki,Gi = Ci['K'], Ci['G']
        f = inclusion.f_inclusion
        # calcul par reccurence des modules
        K,G = Km , Gm
        précision = self.précision
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
        self.n_phases_modèle=2 # Nombre de phases qui doivent etre définies pour calculer le modèle
        
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

  
    def compute_h_behavior(self):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        """
        compatible = self.check_hypothesis()
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        inclusion=self.inclusion
        phase_inclusion = inclusion.list_phases[0]
        Cf, Cm = phase_inclusion.behavior, inclusion.matrix_behavior 
        Km,Gm = Cm['K'], Cm['G']
        num = (3*Km-2*Gm)/(6*Km+2*Gm)
        Kf,Gf = Cf['K'], Cf['G']
        nuf = (3*Kf-2*Gf)/(6*Kf+2*Gf)
        f = inclusion.f_inclusion
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
        return Kh,Gh
    
    

class Autocohérent_IV:
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
        self.n_inclusions = 2 # Nombre d'inclusions de natures différentes 
        self.n_phases_modèle = 3 # Nombre de phases qui doivent etre définies pour calculer le modèle
        
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

  
    def compute_h_behavior(self):
        """
        Calcule le comportement homogénéisé équivalent de la microstructure. Renvoie un dict avec les paramètres calculés. Pour le moment, ne calcul que le module de cisaillement.
        TODO : compléter avec le calcul complet (K et G)
        """
        compatible = self.check_hypothesis()
        if not compatible:
            raise NameError("The microstructure does not match the model hypothesis")
        inclusion=self.inclusion
        phase_inclusion = inclusion.list_phases[0]
        interphase = inclusion.list_phases[1]
        Cf, Cv, Cm = phase_inclusion.behavior, interphase.behavior, inclusion.matrix_behavior 
        Km,Gm = Cm['K'], Cm['G']
        num = (3*Km-2*Gm)/(6*Km+2*Gm)
        Kf,Gf = Cf['K'], Cf['G']
        nuf = (3*Kf-2*Gf)/(6*Kf+2*Gf)
        Kv,Gv = Cv['K'], Cv['G']
        nuv = (3*Kv-2*Gv)/(6*Kv+2*Gv)
        f = inclusion.f_inclusion
        cf = inclusion.compute_f_interphase()

        R3 = 1/(4*np.pi)**(1/3) 
        R2 = (R3**3*(f+cf))**(1/3) 
        R1 = (R3**3*f)**(1/3)  
        
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
        M1[0,1]=(5*(1-nuv))**(-1)*R1**2*(3*b1-7*c1)/(5*(1-2*nuf)) 
        M1[0,2]=(5*(1-nuv))**(-1)*(-12*alpha1/R1**5) 
        M1[0,3]=(5*(1-nuv))**(-1)*4*(f1-27*alpha1)/(15*(1-2*nuf)*R1**3) 
        M1[1,0]=0 
        M1[1,1]=(5*(1-nuv))**(-1)*(1-2*nuv)*b1/(7*(1-2*nuf)) 
        M1[1,2]=(5*(1-nuv))**(-1)*(-20*(1-2*nuv)*alpha1)/(7*R1**7) 
        M1[1,3]=(5*(1-nuv))**(-1)*(-12*alpha1*(1-2*nuv))/(7*(1-2*nuf)*R1**5) 
        M1[2,0]=(5*(1-nuv))**(-1)*R1**5*alpha1/2 
        M1[2,1]=(5*(1-nuv))**(-1)*(-R1**7*(2*a1+147*alpha1))/(70*(1-2*nuf)) 
        M1[2,2]=(5*(1-nuv))**(-1)*d1/7 
        M1[2,3]=(5*(1-nuv))**(-1)*R1**2*(105*(1-nuv)+12*alpha1*(7-10*nuv)-7*e1)/(35*(1-2*nuf)) 
        M1[3,0]=(5*(1-nuv))**(-1)*(-5/6)*(1-2*nuv)*alpha1*R1**3 
        M1[3,1]=(5*(1-nuv))**(-1)*7*(1-2*nuv)*alpha1*R1**5/(2*(1-2*nuf)) 
        M1[3,2]=0 
        M1[3,3]=(5*(1-nuv))**(-1)*e1*(1-2*nuv)/(3*(1-2*nuf)) 
        
        M2=np.zeros(shape=(4,4))
        M2[0,0]=(5*(1-num))**(-1)*c2/3 
        M2[0,1]=(5*(1-num))**(-1)*R2**2*(3*b2-7*c2)/(5*(1-2*nuv)) 
        M2[0,2]=(5*(1-num))**(-1)*(-12*alpha2/R2**5) 
        M2[0,3]=(5*(1-num))**(-1)*4*(f2-27*alpha2)/(15*(1-2*nuv)*R2**3) 
        M2[1,0]=0 
        M2[1,1]=(5*(1-num))**(-1)*(1-2*num)*b2/(7*(1-2*nuv)) 
        M2[1,2]=(5*(1-num))**(-1)*(-20*(1-2*num)*alpha2)/(7*R2**7) 
        M2[1,3]=(5*(1-num))**(-1)*(-12*alpha2*(1-2*num))/(7*(1-2*nuv)*R2**5) 
        M2[2,0]=(5*(1-num))**(-1)*R2**5*alpha2/2 
        M2[2,1]=(5*(1-num))**(-1)*(-R2**7*(2*a2+147*alpha2))/(70*(1-2*nuv)) 
        M2[2,2]=(5*(1-num))**(-1)*d2/7 
        M2[2,3]=(5*(1-num))**(-1)*R2**2*(105*(1-num)+12*alpha2*(7-10*num)-7*e2)/(35*(1-2*nuv)) 
        M2[3,0]=(5*(1-num))**(-1)*(-5/6)*(1-2*num)*alpha2*R2**3 
        M2[3,1]=(5*(1-num))**(-1)*7*(1-2*num)*alpha2*R2**5/(2*(1-2*nuv)) 
        M2[3,2]=0 
        M2[3,3]=(5*(1-num))**(-1)*e2*(1-2*num)/(3*(1-2*nuv)) 
        
        P = np.dot(M2,M1) 
        
        Z12 = P[0,0]*P[1,1]-P[1,0]*P[0,1] 
        Z14 = P[0,0]*P[3,1]-P[3,0]*P[0,1] 
        Z42 = P[3,0]*P[1,1]-P[1,0]*P[3,1] 
        Z23 = P[1,0]*P[2,1]-P[2,0]*P[1,1] 
        Z43 = P[3,0]*P[2,1]-P[2,0]*P[3,1] 
        Z13 = P[0,0]*P[2,1]-P[2,0]*P[0,1] 
    
        A = 4*R3**10*(1-2*num)*(7-10*num)*Z12+20*R3**7*(7-12*num+8*num**2)*Z42+12*R3**5*(1-2*num)*(Z14-7*Z23)+20*R3**3*(1-2*num)**2*Z13+16*(4-5*num)*(1-2*num)*Z43
        B = 3*R3**10*(1-2*num)*(15*num-7)*Z12+60*R3**7*(num-3)*num*Z42-24*R3**5*(1-2*num)*(Z14-7*Z23)-40*R3**3*(1-2*num)**2*Z13-8*(1-5*num)*(1-2*num)*Z43
        C = -R3**10*(1-2*num)*(7+5*num)*Z12+10*R3**7*(7-num**2)*Z42+12*R3**5*(1-2*num)*(Z14-7*Z23)+20*R3**3*(1-2*num)**2*Z13-8*(7-5*num)*(1-2*num)*Z43
            
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
        return Kh,Gh
    
    

R1,R2,e = 5,11,0
f_inclusion = (R1/R2)**3
matrix_behavior={"K":100, "G":150}  
inclusion_behavior = {"K":200, "G":450}           
              
#phase1=Phase(10,1,R1,inclusion_behavior)
#phase_matrix=Phase(10,2,R2,matrix_behavior)
#inclusion1=Inclusion_multiphase(10,[phase1, phase_matrix],matrix_behavior,f_inclusion,['G', 'K'])
#modele2=Autocohérent_III(inclusion1)
#print(modele2.compute_h_behavior())

phase1=Phase(10,1,R1,inclusion_behavior)
interphase=Interphase(10,phase1,e,inclusion_behavior)
phase_matrix=Phase(10,2,R2,matrix_behavior)

inclusion1=Inclusion_multiphase(10,[phase1],matrix_behavior,f_inclusion,['K', 'G'])
inclusion2=Inclusion_multiphase(10,[phase1,interphase],matrix_behavior,f_inclusion,['K', 'G'])

modele1=Autocohérent_II(inclusion1)
modele2=Autocohérent_III(inclusion1)
modele3=Autocohérent_IV(inclusion2)


print(modele1.compute_h_behavior())
print(modele2.compute_h_behavior())
print(modele3.compute_h_behavior())