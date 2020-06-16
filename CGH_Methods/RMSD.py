import MDAnalysis
import MDAnalysis.analysis.rms
import numpy as np
import CGH_Methods.pickleMethods as pickle

'''
l- list
b- boolean
u- universe
s- string
f- float
i- integer
np- numpy array
'''


def test():
    print("RMSD loaded")
    
    
def analysis(lUniverses, uRef,sSelection, lGroups, sSaveAs="RMSDResults", lTimeOffsets=[]):
    lRMSD = []
    for universe in lUniverses:
        R = MDAnalysis.analysis.rms.RMSD(universe, uRef, select=sSelection,
                               groupselections=lGroups,
                               filename ="rmsd.dat")
        R.run()
        if len(lRMSD) < 1:
            lRMSD = R.rmsd.T
        else:
            lRMSD = np.concatenate((lRMSD.T, R.rmsd), axis=0).T
    if len(lTimeOffsets) > 0:
        for x in range(len(lTimeOffsets)):
            fBegIndex = 0
            for i in range(x+1):
                fBegIndex += len(lUniverses[i].trajectory)
            for y in range(fBegIndex, fBegIndex + len(lUniverses[x+1].trajectory)):
                lRMSD[1][y] = lRMSD[1][y] + lTimeOffsets[x]
    dRMSD = {}
    dRMSD['frame'] = lRMSD[0]
    dRMSD['time'] = lRMSD[1]
    dRMSD[sSelection] = lRMSD[2]
    for x in range(len(lGroups)):
        dRMSD[lGroups[x]] = lRMSD[x+3]
    pickle.save(dRMSD, sSaveAs)
