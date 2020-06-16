import MDAnalysis
import numpy as np
import MDAnalysis.analysis.dihedrals as Dihedrals
import CGH_Methods.pickleMethods as pickle
from copy import deepcopy


def analysis(lUniverses, lAtomSelections, lTimeOffsets=[], sSaveAs="DihedralResults"):
    npDihedral = []
    for universe in lUniverses:
        lQuartets = []
        for x in lAtomSelections:
            lQuartets.append(universe.select_atoms(x))
        D = Dihedrals.Dihedral(lQuartets)
        D.run()
        npTempDihedral = D.angles.T
        if len(npDihedral) == 0:
            npDihedral = deepcopy(npTempDihedral)
        else:
            npDihedral = np.concatenate((npDihedral, npTempDihedral), axis=1)
            
    dDihedral = {}
    for x in range(len(lAtomSelections)):
        dDihedral[lAtomSelections[x]] = npDihedral[x]
        
    #Add time points
    dDihedral['time'] = []
    for y in range(len(lUniverses)):
        for x in range(0,len(lUniverses[y].trajectory) * 100, 100):
            if ((y-1) == 0) | ((y-1) > 0):
                dDihedral['time'].append(x + lTimeOffsets[y-1])
            else:
                dDihedral['time'].append(x)
    pickle.save(dDihedral, sSaveAs)