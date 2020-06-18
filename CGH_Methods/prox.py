# Methods for analyzing Residues in proximity to Specific atoms
# 


import MDAnalysis
import MDAnalysis.lib.NeighborSearch as NSearch
import numpy as np
   
def save(dObj, sFilename):
    import pickle
    '''Given an object and a file name, write the object to the file using pickle.'''

    f = open(sFilename, "wb")
    p = pickle.Pickler(f)
    p.dump(dObj)
    f.close()

def load(sFilename):
    import pickle
    '''Given a file name, load and return the object stored in the file.'''
    
    f = open(sFilename, "rb")
    u = pickle.Unpickler(f)
    dObj = u.load()
    f.close()
    return dObj

def analysis(lUniverses, searchAtoms, radius, fileName):
#Inputs:
#       a list of MD Analysis Universes,
#       a search string (in MD Analysis atom selection grammar)
#       Radius (in Angstroms) around the search atoms in which to search
#       string to call the analysis result file
#Outputs:
#       a pickle file containing a list of dictionaries with the results in the following format:
#       [{information about a residue that is whithin the search radius}, {}, {}, {}]

    resList = []
    for u in lUniverses:
        SearchClass = NSearch.AtomNeighborSearch(u.atoms)
        nr = []
        for ts in u.trajectory:
            nr.append(SearchClass.search(u.select_atoms(searchAtoms), radius, 'R'))
        nr = np.array(nr)
        for a in range(len(nr)):
            for i in range(len(nr[a])):
                resList.append(nr[a][i])    
    resList = np.array(resList)
    unRes = np.unique(resList, return_counts=True)
    unRes = np.array(unRes)
    unRes = unRes.T
    unRes = unRes[np.argsort(unRes[:, 1])] #sorts into ascending order of counts
    unRes = np.flip(unRes, axis=0) # flips to descending order of counts
    lResidues = []
    for i in range(1,len(unRes)):
        dRes = {}
        dRes['segid'] = unRes[i][0].segid
        dRes['resname'] = unRes[i][0].resname
        dRes['resid'] = unRes[i][0].resid
        dRes['count'] = unRes[i][1]
        dRes['countPercentage'] = unRes[i][1]/len(lUniverses[0].trajectory)
        dRes['full_string'] = str(dRes['segid']) + " " + str(dRes['resname']) + " " + str(dRes['resid'])
        lResidues.append(dRes)
    save(lResidues, fileName)

def plot(fileName, percentage=False):
# given an analysis result file name, plot the results 
# works only for a single results file
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    
    lResidues = load(fileName)
    # Make Arrays for Plotting
    lX = []
    lY = []
    for i in range(20):
        lX.append(lResidues[i]['full_string'])
        if percentage:
            lY.append(lResidues[i]['countPercentage'])
        else:
            lY.append(lResidues[i]['count'])

    # Plot
    import matplotlib.pyplot as plt
    fig = plt.figure(1)
    plt.suptitle('Residues neighboring Search Atoms (within (radius)A)',y = 1, fontsize = 'x-large')
    gridspec.GridSpec(2,4)
    ax = plt.subplot2grid((2,4), (0,0), colspan=4, rowspan=2)
    ax.plot(lX, lY,  '.', ms = 10, color='black')
    ax.set_xlabel("Residues")
    ax.set_ylabel("Counts")
    plt.xticks(lX, rotation=45, horizontalalignment='right')
    fig.set_size_inches(w=11,h=4)
    plt.show()
    
def compResList(lSystems):
#Inputs:
#        list of results lists
#Outputs:
#        combined list of results (without counts)
# takes only the first 20 residues from each results list
    lCombined = []
    for ls in lSystems:
        if len(lCombined) > 0:
            for res in range(20):
                check = 0
                for x in range(len(lCombined)):
                    if (ls[res]['resid'] == lCombined[x]['resid']) & (ls[res]['segid'] == lCombined[x]['segid']):
                        check += 1
                if check == 0:
                    lCombined.append(ls[res])
        else:
            for res in range(20):
                lCombined.append(ls[res])
    return lCombined



def xval(lSystems, bChain=False):
#Inputs:
#        list of lists of results
#        include chain information? t/f
# Outputs:
#        x values of combined list of residues for plotting
    lCombined = compResList(lSystems)
    lX = []
    if bChain:
        for i in range(len(lCombined)):
            lX.append(str(lCombined[i]['segid']) + " " + str(lCombined[i]['resname']) + " " + str(lCombined[i]['resid']))
    else:
        for i in range(len(lCombined)):
            lX.append(str(lCombined[i]['resname']) + " " + str(lCombined[i]['resid']))
    return lX

def yval(lSystems, bPercentage=False):
#Inputs:
#        list of lists of results
#        plot as fraction of total frames? t/f
# Outputs:
#        y values of combined list of residues for plotting
    lCombined = compResList(lSystems)
    lY = []
    for ls in lSystems:
        lTempY = []
        for res in lCombined:
            lTempY.append(find(ls, res['resid'], res['segid'], bPercentage))
        lY.append(lTempY)
    return lY

        
def find(ls, resid, segid, percentage= False):
#returns t/f that the resid/segid is found in ls
    bFound = False
    for dic in ls:
        if (dic['resid']  == resid) & (dic['segid'] == segid):
            if percentage:
                return dic['countPercentage']
            else: 
                return dic['count']
            bFound = True
    if bFound == False:
        return 0
    
def merge(ls):
    from operator import itemgetter as ig
    from copy import deepcopy
    lReturn = []
    lTemp = []
    lFrames = []
    for l in ls:
        lFrames.append(round(l[0]['count']/l[0]['countPercentage']))
    for l in ls:
        if len(lTemp) < 1:
            for res in l:
                lTemp.append(deepcopy(res))
        else:
            for res in l:
                check = 0
                for x in range(len(lTemp)):
                    if res['resid'] == lTemp[x]['resid']:
                        check += 1
                        lTemp[x]['count'] += res['count']
                if check == 0:
                    lTemp.append(res)
    lReturn = deepcopy(lTemp)
    for res in lReturn:
        res['countPercentage'] = res['count'] / sum(lFrames)
    lReturn.sort(key=ig('count'), reverse=True)
    return lReturn     
    
    
