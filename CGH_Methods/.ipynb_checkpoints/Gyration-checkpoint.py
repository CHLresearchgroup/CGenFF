import numpy as np
import MDAnalysis
import CGH_Methods.pickleMethods as pickle

def analysis(lUniverses, lGroupSelection, lTimeOffsets=[], sSaveAs="GyrationResults"):
    Rgyr = []
    u = 0
    for universe in lUniverses:
        lGroups = []
        for group in lGroupSelection:
            lGroups.append(universe.select_atoms(group))
        for ts in universe.trajectory:
            tempRgyr = []
            if u < 1:
                tempRgyr.append(universe.trajectory.time)
            else:
                tempRgyr.append(universe.trajectory.time + lTimeOffsets[u-1])
            for group in lGroups:
                tempRgyr.append(group.radius_of_gyration())
            Rgyr.append(tempRgyr)
        u += 1
    Rgyr = np.array(Rgyr).T
    dGyration = {}
    dGyration['time'] = Rgyr[0]
    for x in range(len(lGroupSelection)):
        dGyration[lGroupSelection[x]] = Rgyr[x+1]
    pickle.save(dGyration, sSaveAs)