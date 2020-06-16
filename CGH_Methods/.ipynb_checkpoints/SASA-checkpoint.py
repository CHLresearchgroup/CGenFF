import gromacs.tools as tools
import gromacs.formats as wrap
import os
import numpy as np
import CGH_Methods.pickleMethods as pickleMethods
from copy import deepcopy
import MDAnalysis

def analysis(lUniverses, lTpr, lXtc, sSurface, lOutput, lTimeOffsets=[], sSaveAs='SASAResults'):
    makeIndex(lUniverses, sSurface, lOutput)
    lSASA = []
    npSASA = np.array(lSASA)
    lOutputInts = []
    for y in range(1,len(lOutput)+1):
        lOutputInts.append(y)
    for x in range(len(lTpr)):
        cSASA = tools.G_sas(s=lTpr[x], f=lXtc[x], n='SASAIndex.ndx', o='SASA.xvg', surface=["0"], output=lOutputInts)
        cSASA.run()
        wSASA = wrap.XVG()
        wSASA.read(filename='SASA.xvg')
        npTempSASA = wSASA.array
        if x == 0:
            npSASA = deepcopy(npTempSASA)
        else:
            npSASA = np.concatenate((npSASA, npTempSASA), axis=1)       
    dSASA = createDictionary(lUniverses, npSASA, sSurface, lOutput, lTimeOffsets)
    #save dictionary
    pickleMethods.save(dSASA, sSaveAs)
    
    
def makeIndex(lUniverses, surface, lSelections):
    os.system('rm SASAIndex.ndx')
    indexWriter = MDAnalysis.selections.gromacs.SelectionWriter('SASAIndex.ndx', mode='w', numterms=None, preamble=None)
    indexWriter.write(lUniverses[0].select_atoms(surface), name="0")
    for x in range(1,len(lSelections) +1):
        indexWriter.write(lUniverses[0].select_atoms(lSelections[x -1]), name=str(x))
    indexWriter.close()
    
def createDictionary(lUniverses, npSASA, sSurface, lOutput, lTimeOffsets=[]):
    # offset the times if multiple runs
    if len(lTimeOffsets) > 0:
        for x in range(len(lTimeOffsets)):
            fBegIndex = 0
            for i in range(x+1):
                fBegIndex += len(lUniverses[i].trajectory)
            for y in range(fBegIndex, fBegIndex + len(lUniverses[x+1].trajectory)):
                npSASA[0][y] = npSASA[0][y] + lTimeOffsets[x]
    
    #convert nparray to dictionary
    dSASA = {}
    dSASA['time'] = npSASA[0]
    dSASA[sSurface] = npSASA[1]
    for x in range(len(lOutput)):
        dSASA[lOutput[x]] = npSASA[x+2]
    return dSASA
    