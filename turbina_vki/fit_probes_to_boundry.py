__author__ = 'wgryglas'

"""
Script calculating probes positions that would fit to the mesh boundary
 (provided set of points usually is not transofrmed)
"""

def perform(dirs, files, par):
    from RedBoundaryReader import readBoundaryNodes, write_coordinates
    from wg.tools.geom import fit_coords

    bX, bY = readBoundaryNodes(files.boundary_coords)
    pX, pY = readBoundaryNodes(files.probes_base_coords)

    X, Y, params = fit_coords(pX, pY, bX, bY)

    write_coordinates(files.probes_coords, X, Y)





import unittest
class ComapreFit(unittest.TestCase):
    def setUp(self):
        import settings
        perform(settings.dirs, settings.files, settings.par)

    def runTest(self):
        try:
            from settings import files
            import matplotlib.pyplot as plt
            from RedBoundaryReader import readBoundaryNodes

            bX, bY = readBoundaryNodes(files.boundary_coords)
            X, Y = readBoundaryNodes(files.probes_coords)
            pX, pY = readBoundaryNodes(files.probes_base_coords)

            plt.figure()
            plt.plot(bX, bY, 'r.')
            plt.plot(pX, pY, 'g.', markersize=20)
            plt.plot(X, Y, 'b.', markersize=20)

            for i in range(len(X)):
                plt.annotate(s='', xy=(X[i], Y[i]), xytext=(pX[i], pY[i]), arrowprops=dict(arrowstyle='->'))

            plt.show()
        except Exception as e:
            self.fail("Compare fit exception shold not appear, msg:\n", e)
