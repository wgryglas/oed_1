
filePath = "/home/wgryglas/AVIO/data/pszaltys_data_naca0012/fin/aoa2_2_ma_0_2_surf.dat"

from RedBoundaryReader import RedBoundaryReader


data = RedBoundaryReader(filePath).read()

print data.variables

import matplotlib.pyplot as plt

d = data["surface1"]

plt.plot(d[:, 0], d[:, 1], ".")
plt.show()

import numpy as np

# np.savetxt("/home/wgryglas/AVIO/data/pszaltys_data_naca0012/fin/boundary.dat", d[:,:2])

