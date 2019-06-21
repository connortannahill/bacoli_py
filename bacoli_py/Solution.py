class Solution:

    """An object to hold information from calls to bacoli_py.

    Includes solution information as well as information about the 
    computation performed by BACOL(R)I. These objects are returned to the user
    following a successful call to bacoli_py.solve().
    """

    def __init__(self, tspan, xspan, u, ux=None):
        """Solution object.

        Parameters
        ----------
        tspan : castable to floating point ndarray
            Points along the temporal domain at which the PDE system has been computed.

        xspan : castable to floating point ndarray
            Points along x-axis at which the solution of the PDE system has been computed.
        
        u : ndarray with shape=(npde, len(tspan), len(xspan))
            Solution output by BACOL(R)I at a points in time and space.

        ux : ndarray with shape=(npde, len(tspan), len(xspan))
            Solution derivatives output by BACOL(R)I at a points in time and space.
        """

        self.xspan = xspan
        self.tspan = tspan
        self.u = u 
        self.ux = ux
