HowTo add a model

Models are defined in the 'classes_v2.py' file. A template for your models is defined below.

class YourModel(Model):
    """
    An optional description of your model.
    """
    def __init__(self):
        self.type_inclusion = 0 for spheres, 1 for ellipsoids
        self.behavior_condition = 'isotropic' or 'anisotropic'
        self.n_inclusion = maximum number of inclusions of different nature allowed
        self.interphase = True if the model needs an inclusion with an interphase
        self.name = 'YourModel name'
        
    def compute_behavior(self, Cm, inclusion_behaviors):
        # Cm is a 'behavior dictionnary' containing the matrix parameters.
        # If the matrix is isotropic, you can recover its bulk, shear, Young and Poisson modulus with:
        Km, Gm, Em, num = Cm['K'], Cm['G'], Cm['E'], Cm['nu']
        # If the matrix is anisotropic, you can revover its compliance and stifness matrices (arrays) with:
        S, C = Cm['S'], Cm['C']
        # inclusion_behaviors is a list of tuples with as many tuples as inclusions in the microstructure
        # Each tuple contains three values (Cf, f, aspect_ratio)
        # Cf is a 'behavior dictionnary' (same as Cm)
        # f is a float (volume fraction)
        # aspect_ratio is a tuple containing the two aspect ratios of the inclusion ((1,1) for a sphere)
        # The function should return a 'behavior dictionnary'
        return {'K': K_model, 'G': G_model}
        
