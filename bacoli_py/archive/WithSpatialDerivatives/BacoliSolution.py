"""
    Class to contain information from calls to bacoli_py's solve method.

    @author: Connor Tannahill
"""

class BacoliSolution:
    """ An object to hold information from calls to bacoli_py.

    Includes solution information as well as information about the 
    computation performed by BACOLI. These objects are returned to the user
    following a successful call to bacoli_py.solve().
    """

    def __init__(self, output_points, output_times, solution, bacoli_info, 
                 spatial_derivatives=None):
        """Construct BacoliSolution object.

        Args:
            output_points (:obj:`list` of (float)): Points along x-axis at
                which the solution of the PDE system has been computed.
            
            output_times (:obj:`list` of (float)): Points along the tempoaral
                domain at which the PDE system has been computed.
            
            solution (:obj:`list` of :obj:`list` of (float)): Solution output
                by BACOLI at a points in time and space.
            
            bacoli_info (:obj:`dict`): Dictionary containing information about
                the parameters used to compute this solution, as well as info.
                about the computation performed by BACOLI.
        """

        self.output_points = output_points
        self.output_times = output_times
        self.solution = solution
        self.bacoli_info = bacoli_info
        self.spatial_derivatives = spatial_derivatives
    
    def print_to_file(self, file_names=None):
        """ Print PDE solutions to ouput files.

        For each PDE in the system an output file is opened and its solution 
        is copied there.

        Note: 
            If file_name argument is given, list must contain as many strings
            as there are PDE's to be output.

        Args:
            file_name (:obj:`list` of :obj:`str`, optional): Names of output 
            files.
        """

        npde = len(self.solution)

        if file_names is not None and len(file_names) != npde:
            raise ValueError("Error: Not enough file_names given. Must provide "
                             + "one for each PDE in system.")

        # Print solution to file.
        for i in range(npde):

            # Open output file.
            if file_names is None:
                fname = "Points{0:d}".format(i+1)
            else:
                fname = file_names[i]
            
            output_file = open(fname, "w")

            # Print to file.
            if self.spatial_derivatives is None:
                for j in range(len(self.output_times)):
                    for k in range(len(self.output_points)):
                        out_str = '{0: <18.18f} {1: <18.18f} {2: <18.18f}\n'.format(
                            self.output_points[k], self.output_times[j], 
                            self.solution[i][j][k])

                        output_file.write(out_str)
            else:
                for j in range(len(self.output_times)):
                    for k in range(len(self.output_points)):
                        for n in range(len(self.spatial_derivatives)):
                            out_str = '{0: <18f} {1: <18f} {2: <18f} {3: <18f}\n'.format(
                                self.output_points[k], self.output_times[j],
                                self.solution[i][j][k], self.spatial_derivatives[i][j][k][n])

                            output_file.write(out_str) 