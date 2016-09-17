def getSortedCoords(bX, bY):
    import numpy as np

    bSet = np.array([bX, bY]).T

    pnts = list()

    left = [p for p in bSet]

    start = np.argmax(bX) # this will set trailing edge as 0 in parameter space

    pnts.append(bSet[start])
    left.pop(start)

    for i in range(1, len(bX)):
        p = pnts[-1]
        minD = np.inf
        res = None
        resId = None
        for i, p2 in enumerate(left):
            v = p2 - p
            dist = v.dot(v)
            if dist < minD:
                minD = dist
                res = p2
                resId = i

        if res is not None:
            pnts.append(res)
            left.pop(resId)

    pnts = np.array(pnts)

    return pnts


from settings import dirs, files, organizer

bX, bY = organizer.load_mesh_boundary_coordinates()

sortedPnts = getSortedCoords(bX, bY)[::-1] #assure right order

with open(dirs.root+"/orderedbndcoords.dat","w") as f:
    for p in sortedPnts:
        f.write("(%lf,\t%lf),\n" % (p[0], p[1]))
