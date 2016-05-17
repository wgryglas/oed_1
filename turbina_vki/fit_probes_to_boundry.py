__author__ = 'wgryglas'

"""
Script calculating probes positions that would fit to the mesh boundary
 (provided set of points usually is not transofrmed)
"""

# from RedBoundaryReader import readBoundaryNodes, write_coordinates
# from wg.tools.geom import fit_coords
#
# import matplotlib.pyplot as plt
# from settings import *
#
# boundary_coords = "/home/wgryglas/AVIO/data/tunel_vki/boundaryNodesCoordinates.dat"
# probes_coords = "/home/wgryglas/AVIO/data/tunel_vki/probe_base_points"
# out_file = "/home/wgryglas/AVIO/data/tunel_vki/probe_points"
#
#
# bX, bY = readBoundaryNodes(files.boundary_coords)
#
# pX, pY = readBoundaryNodes(files.probes_coords)
#
#
# X, Y, params = fit_coords(pX, pY, bX, bY)
#
# print params
#
# # write_coordinates(out_file, *fit_coords(pX, pY, bX, bY))
#
# #TEST
# #X, Y = readBoundaryNodes(out_file)
#
# plt.figure()
# plt.plot(bX, bY, 'r.')
# plt.plot(pX, pY, 'g.', markersize=20)
# plt.plot(X, Y, 'b.', markersize=20)
#
# for i in range(len(X)):
#     plt.annotate(s='', xy=(X[i],Y[i]), xytext=(pX[i], pY[i]), arrowprops=dict(arrowstyle='->'))
#
# plt.show()