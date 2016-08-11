

# line = 'VARIABLES = "X", "Y", "rho", "rhoE", "rhoU", "rhoV", "nut", "p", "mach", "mach_iso"'
#
# import re
#
# variables=[]
# if re.match('.*VARIABLES.*', line):
#         for var in re.findall('"([a-zA-Z,0-9,_]*)"',line):
#             variables.append(var)
#
# print ", ".join(variables)


# import re
# name = "coeff_0-65_mach_iso.npz"
#
# datas = re.findall(r'coeff_([0-9]+-[0-9]*)_([0-9,aA-zZ]+).npz', name)
#
# print datas[0][0]
# print datas[0][1]

# import re
# text = "mach_iso_8_0-52.dat"
#
# groups = re.findall(r'([0-9,aA-zZ]+[_[0-9,aA-zZ]*]*)_([0-9]+)_([0-9]+[-[0-9]+]*).dat', text)
# print groups

# import numpy as np
#
# v1 = np.array([[1,2,3,4,5],[1,4,8,16,32]]).T
# print v1
#
# v2 = np.array([-1,-2,-3,-4,-5])
#
# print v1 * v2[:, np.newaxis]
#
# l = np.array([0, 0])
#
# np.add.at(v2, l, np.array([10,12]))
#
# print v2

# s = '1 1 "2inl2et2"'
# import re
#
# print re.findall(r'([0-9]+)\s([0-9]+)\s\"(([aA-zZ]|[0-9])+)\"', s)[0][:3]

l = [1,2,3,4]
print l[-2:]