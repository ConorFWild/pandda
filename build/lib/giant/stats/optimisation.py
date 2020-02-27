from scitbx.array_family import flex
from scitbx import simplex


class _LeastSquaresFitter(object):

    def __init__(self, x_values, ref_values, scl_values=None, weights=None):
        """Calculate the optimum scaling of x_values -> ref_values (or scl_values -> ref_values, if given) for a target function"""

        self.x_values = flex.double(x_values)

        self.ref_values = flex.double(ref_values)

        if scl_values is not None:
            self.scl_values = flex.double(scl_values)
        else:
            self.scl_values = self.x_values

        if weights is not None:
            self.weight_array = flex.double(weights)
        else:
            self.weight_array = flex.double(self.x_values.size(), 1.0/self.x_values.size())

        assert self.x_values.size() == self.ref_values.size()
        assert self.x_values.size() == self.scl_values.size()
        assert self.x_values.size() == self.weight_array.size()

        # Initialise simplex
        self.starting_values = self.initialise_parameters()
        # Calculate scaling
        self.run(initial_simplex=self.starting_simplex)

    def run(self, initial_simplex):
        """Calculate scaling"""
        # Optimise the simplex
        self.optimised = simplex.simplex_opt(dimension = len(initial_simplex[0]),
                                             matrix    = initial_simplex,
                                             evaluator = self)
        # Extract solution
        self.optimised_values = self.optimised.get_solution()

        # Scaled input values
        self.out_values = self.transform(values=self.scl_values)
        # Calculate rmsds
        self.unscaled_rmsd = self.rmsd_to_ref(values=self.scl_values)
        self.scaled_rmsd = self.rmsd_to_ref(values=self.out_values)

        return self.optimised_values

    def rmsd_to_ref(self, values, sel=None):
        """Calculate the rmsd between a vector and the reference values"""
        assert len(values) == len(self.ref_values)
        if sel:
            ref = self.ref_values.select(sel)
            val = values.select(sel)
        else:
            ref = self.ref_values
            val = values

        return (ref-val).norm()/(ref.size()**0.5)

    def target(self, vector):
        """Target function for the simplex optimisation"""
        scaled = self.transform(values=self.scl_values, params=vector)
        diff = (scaled-self.ref_values)
        diff_sq = diff*diff
        result = flex.sum(self.weight_array*diff_sq)
        return result

    def transform(self, values, params=None):
        """Function defining how the fitting parameters are used to transform the input vectors"""
        if params is None:
            params=self.optimised_values
        values = flex.double(values)
        assert values.size() == self.x_values.size()
        return self._scale(values=values, params=params)

    def new_x_values(self, x_values):
        self.x_values = flex.double(x_values)
        self.ref_values = None
        self.scl_values = None
        self.out_values = None
        self.weight_array = None

class LinearScaling(_LeastSquaresFitter):

    def initialise_parameters(self):
        """Initialise starting simplex"""
        v0 = 0.0    # 0th order - offset
        v1 = 1.0    # 1st order - scale
        self.starting_simplex = [    flex.double([v0,    v1]),
                                     flex.double([v0,    v1+0.1]),
                                     flex.double([v0+0.1,v1])      ]
        return [v0,v1]

    def _scale(self, values, params):
        v0,v1 = params
        out = v0 + values*v1
        return out

class ExponentialScaling(_LeastSquaresFitter):

    def initialise_parameters(self):
        """Initialise starting simplex"""
        v1 = 0.0    # 1st order - scale
        self.starting_simplex = [    flex.double([v1-10]),
                                     flex.double([v1+10])
                                ]
        return [v1]

    def _scale(self, values, params):
        v1, = params
        out = flex.exp(v1*self.x_values)*values
        return out

