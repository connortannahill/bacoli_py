import numpy as np
from numpy import array
import inspect

class ProblemDefinition:
    """Defines a system of 1D parabolic PDEs"""


    def __init__(self, npde, f, bndxa, bndxb, uinit, derivf=None, 
                 difbxa=None, difbxb=None):
        """
        Parameters
        ----------
        npde : int
            The number of PDE's in this system.
        f : callback function
            System of PDE's to be solved.
        bndxa : callback function 
            Left boundary conditions for the system of PDEs.
        bndxb : callback function 
            Right boundary conditions for the system of PDEs.
        uinit : callback function
            Initial conditions for the set of PDEs.
        derivf : callback function 
            Partial derivatives of the PDE system.
        difbxa : callback function 
            Partial derivatives of the left boundary conditions.
        difbxb : callback function 
            Partial derivatives of the right boundary conditions.

        Raises
        ------
        ValueError
            If npde is the wrong type
        ValueError
            If only one of difbxa, difbxb are provided.
        """

        try:
            self.npde = array(npde, dtype=np.int)
        except ValueError:
            raise ValueError('npde must be an integer value.')

        self.f = f 

        self.bndxa = bndxa

        self.bndxb = bndxb 

        self.uinit = uinit 

        if (difbxa == None and difbxb != None) \
                or (difbxa != None and difbxb == None):
            raise ValueError('either neither or both of difbxa and '
                           + 'difbxb must be provided.')

        # Check if any of the three optional callback functions have been
        # provided. If not, set to dummy function and set appropriate flag.
        if derivf != None:
            self.derivf = derivf
            self.is_derivf = True
        else:
            self.derivf = dummy_derivf
            self.is_derivf = False

        if difbxa != None:
            self.difbxa = difbxa
            self.is_difbxa = True
        else:
            self.difbxa = dummy_difbx
            self.is_difbxa = False

        if difbxb != None:
            self.difbxb = difbxb
            self.is_difbxb = True
        else:
            self.difbxb = dummy_difbx
            self.is_difbxb = False

# Dummy callback functions to be used if not defined by user. Fortran
# interface requires that all arguments always be passed.
def dummy_derivf(t, x, u, ux, uxx, dfdu, dfdux, dfduxx):
    pass

def dummy_difbx(t, u, ux, dbdu, dbdux, dbdt):
    pass

