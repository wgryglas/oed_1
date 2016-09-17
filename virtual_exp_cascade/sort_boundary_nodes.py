def getSortedCoords(bX, bY):
    import numpy as np

    bSet = np.array([bX, bY]).T

    pnts = list()

    left = [p for p in bSet]

    pnts.append(bSet[0])
    left.pop(0)

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
