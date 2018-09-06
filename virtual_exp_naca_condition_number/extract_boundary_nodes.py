from RedBoundaryReader import RedBoundaryReader
from RedBoundaryReader import readBoundaryNodes
from wg.occ.curve2D import *
from wg.tools.function import empty
import matplotlib.pyplot as plt
from bearded_octo_wookie.RED.RedReader import GetReader

"""
Script which reads tecplot(?) boundary data and thus extracts boundary nodes coordinates. Those nodes are later
saved in single file in format:
N=...
X1 Y1
X2 Y2
.....
where N is number of nodes
"""
__author__ = 'wgryglas'

def perform(dirs, files, params, organizer):
    boundaryData = RedBoundaryReader(files.boundary_source).read()
    boundaryData.write(files.boundary_coords, ['X', 'Y'])


# READ SAVED DATA
# x, y = readBoundaryNodes(files.boundary_coords)
# plt.figure()
# plt.plot(x, y, '.')



# SOME TESTS WITH STEP BOUNDARY FILE
# bounds1 = Bounds().set_from_data(x, y)
# geomFile = "/home/wgryglas/AVIO/data/tunel_vki/geom/palisada.step"
# getFile = "/home/wgryglas/Desktop/kaskada_ns.get"
#
# curves = GetReader(getFile)

# curves = Curves(geomFile)
# scale = 1./curves.bounds().size()[0]
# angle = 26.7/180*np.pi
#
# curves.translate(-curves.center())
# curves.translate([])
# curves.scale(scale)
# curves.mirror(axis=[0, 0, 0, 0, 1, 0])
# curves.rotate(angle, axis=[0, 0, 0, 0, 0, 1])
# curves.mirror(axis=[0, 0, 0, 1, 0, 0]).mirror(axis=[0, 0, 0, 0, 1, 0])
# curves.rotate(-angle, [0, 0, 0, 0, 0, 1])
#
# bounds2 = curves.bounds()
# scale = bounds1.size()[0]/bounds2.size()[0]
# curves.scale(scale)
# bounds2 = curves.bounds()
# curves.translate(-(bounds2.max-bounds1.max))


#xt, yt, zt = curves.tesselatation()


# xy = np.array(curves.getXYList(np.linspace(0., 1., 100), 1))
#
#
# # plt.figure()
# # plt.plot(xt, yt, 'r.')
# plt.plot(xy[:, 0], xy[:, 1], 'r.')
# plt.show()