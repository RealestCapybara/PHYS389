import numpy as np
from Chain import Chain
from Pendulum import Pendulum
import matplotlib.pyplot as plt
import pickle
import tomli as toml
from operator import itemgetter
import os


def checkSettings(settings):
    try:
        sys = settings["system"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'system' table")
    try:
        filename = sys["filename"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'filename' value")
    try:
        anim = sys["anim"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'anim' value")
    try:
        plot = sys["plot"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'plot' table")
    try:
        x = plot["x"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'x' value")
    try:
        v = plot["v"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'v' value")
    try:
        pos = plot["pos"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'pos' value")
    try:
        vel = plot["vel"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'vel' value")
    try:
        pole = plot["pole"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'pole' value")
    try:
        KE = plot["KE"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'KE' value")
    try:
        PE = plot["PE"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'PE' value")
    try:
        ME = plot["ME"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'ME' value")

    try:
        totalKE = plot["totalKE"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'totalKE' value")
    try:
        totalPE = plot["totalPE"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'totalPE' value")
    try:
        totalME = plot["totalME"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'totalME' value")

    try:
        LM = plot["LM"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'LM' value")
    try:
        AM = plot["AM"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'AM' value")
    try:
        totalLM = plot["totalLM"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'totalLM' value")
    try:
        totalAM = plot["totalAM"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'totalAM' value")
    
    if not isinstance(filename, str):
        raise TypeError("'filename' must be a str")
    if not isinstance(anim, bool):
        raise TypeError("'anim' must be a boolean")
    if not isinstance(x, list):
        raise TypeError("'x' must be a list")
    if len(x) != 3:
        raise ValueError("'x' must have len=3")
    if not all(isinstance(val, bool) for val in x):
        raise TypeError("'x' must only contain boolean values")
    if not isinstance(v, list):
        raise TypeError("'v' must be a list")
    if len(v) != 3:
        raise ValueError("'v' must have len=3")
    if not all(isinstance(val, bool) for val in v):
        raise TypeError("'v' must only contain boolean values")
    if not isinstance(pos, list):
        raise TypeError("'pos' must be a list")
    if len(pos) != 2:
        raise ValueError("'pos' must have len=2")
    if not all(isinstance(val, bool) for val in pos):
        raise TypeError("'pos' must only contain boolean values")
    if not isinstance(vel, list):
        raise TypeError("'vel' must be a list")
    if len(vel) != 2:
        raise ValueError("'vel' must have len=2")
    if not all(isinstance(val, bool) for val in vel):
        raise TypeError("'vel' must only contain boolean values")
    
    if not isinstance(pole, bool):
        raise TypeError("'pole' must be a boolean value")
    if not isinstance(KE, bool):
        raise TypeError("'KE' must be a boolean value")
    if not isinstance(PE, bool):
        raise TypeError("'PE' must be a boolean value")
    if not isinstance(ME, bool):
        raise TypeError("'ME' must be a boolean value")
    if not isinstance(totalKE, bool):
        raise TypeError("'totalKE' must be a boolean value")
    if not isinstance(totalPE, bool):
        raise TypeError("'totalPE' must be a boolean value")
    if not isinstance(totalME, bool):
        raise TypeError("'totalME' must be a boolean value")

    if not isinstance(LM, list):
        raise TypeError("'LM' must be a list")
    if len(LM) != 3:
        raise ValueError("'LM' must have len=3")
    if not all(isinstance(val, bool) for val in LM):
        raise TypeError("'LM' must only contain boolean values")

    if not isinstance(AM, list):
        raise TypeError("'AM' must be a list")
    if len(AM) != 3:
        raise ValueError("'AM' must have len=3")
    if not all(isinstance(val, bool) for val in AM):
        raise TypeError("'AM' must only contain boolean values")

    if not isinstance(totalLM, list):
        raise TypeError("'totalLM' must be a list")
    if len(totalLM) != 3:
        raise ValueError("'totalLM' must have len=3")
    if not all(isinstance(val, bool) for val in totalLM):
        raise TypeError("'totalLM' must only contain boolean values")

    if not isinstance(totalAM, list):
        raise TypeError("'totalAM' must be a list")
    if len(totalAM) != 3:
        raise ValueError("'totalAM' must have len=3")
    if not all(isinstance(val, bool) for val in totalAM):
        raise TypeError("'totalAM' must only contain boolean values")

    return sys

def plotValues(data):
    """Returns a dict with the values to plot."""
    
    #each key I want to exist in the returned dictionary (depends on config file)
    #the first key list are values for each pendulum
    pkeyList = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'theta', 'phi', 'thetadot',
               'phidot', 'whichPole', 'KE', 'PE', 'ME', 'LMx', 'AMx', 'LMy', 
               'AMy', 'LMz', 'AMz', 'absLM', 'absAM']
    #the second list of keys are values of the entire system
    sysKeyList = ['totalKE', 'totalPE', 'totalME',  'totalLMx', 'totalAMx', 
                  'totalLMy', 'totalAMy', 'totalLMz', 'totalAMz', 'absTotalLM',
                  'absTotalAM', 'totalAbsLM', 'totalAbsAM']

    #creates an empty dictionary containing lists attributed to all possible keys
    #uses comprehension to add numbers to the keys for each pendula
    FrameDict = {f'{val}{i}':[] for i in range(len(data[0].pList))
                 for val in pkeyList}
    #adds on the keys for the system. These don't need numbers added to the end
    #so need to be treated differently.
    FrameDict.update({k:[] for k in sysKeyList})

    for i, v in enumerate(data):
        Pos = v.Positions()
        Vel = v.Velocities()
        KE = v.KE()
        PE = v.PE()
        ME = [i+j for i,j in zip(KE,PE)]
        LM = v.LMom()
        AM = v.AMom()
        
        for i, p in enumerate(v.pList):
            FrameDict[f'x{i}'].append(Pos[i][0])
            FrameDict[f'y{i}'].append(Pos[i][1])
            FrameDict[f'z{i}'].append(Pos[i][2])

            FrameDict[f'vx{i}'].append(Vel[i][0])
            FrameDict[f'vy{i}'].append(Vel[i][1])
            FrameDict[f'vz{i}'].append(Vel[i][2])

            FrameDict[f'vx{i}'].append(Vel[i][0])
            FrameDict[f'vy{i}'].append(Vel[i][1])
            FrameDict[f'vz{i}'].append(Vel[i][2])

            FrameDict[f'LMx{i}'].append(LM[i][0])
            FrameDict[f'LMy{i}'].append(LM[i][1])
            FrameDict[f'LMz{i}'].append(LM[i][2])
            FrameDict[f'absLM{i}'].append(np.linalg.norm(LM[i]))

            FrameDict[f'AMx{i}'].append(AM[i][0])
            FrameDict[f'AMy{i}'].append(AM[i][1])
            FrameDict[f'AMz{i}'].append(AM[i][2])
            FrameDict[f'absAM{i}'].append(np.linalg.norm(AM[i]))

            FrameDict[f'KE{i}'].append(KE[i])
            FrameDict[f'PE{i}'].append(PE[i])
            FrameDict[f'ME{i}'].append(ME[i])

            FrameDict[f'theta{i}'].append(p.pos[0])
            FrameDict[f'phi{i}'].append(p.pos[1])
            FrameDict[f'thetadot{i}'].append(p.vel[0])
            FrameDict[f'phidot{i}'].append(p.vel[1])

        FrameDict['totalKE'].append(sum(KE))
        FrameDict['totalPE'].append(sum(PE))
        FrameDict['totalME'].append(sum(ME))

        totLM = sum(LM)
        FrameDict['totalLMx'].append(totLM[0])
        FrameDict['totalLMy'].append(totLM[1])
        FrameDict['totalLMz'].append(totLM)
        FrameDict['absTotalLM'].append(np.linalg.norm(totLM))
        FrameDict['totalAbsLM'].append(sum([np.linalg.norm(p) for p in LM]))

        totAM = sum(AM)
        FrameDict['totalAMx'].append(totAM[0])
        FrameDict['totalAMy'].append(totAM[1])
        FrameDict['totalAMz'].append(totAM[2])
        FrameDict['absTotalAM'].append(np.linalg.norm(totAM))
        FrameDict['totalAbsAM'].append(sum([np.linalg.norm(p) for p in AM]))

    return FrameDict

def filterDiagnosticData(dictionary, settings):
    #creates a reference dictionary to connect settings keys with 
    #data dictionary keys
    refDict = {'x': ('x', 'y', 'z'),
               'v': ('vx', 'vy', 'vz'),
               'pos': ('theta', 'phi'),
               'vel': ('thetadot', 'phidot'),
               'pole': ('whichPole',),
               'KE': ('KE',),
               'PE': ('PE',),
               'ME': ('ME',),
               'totalKE': ('totalKE',),
               'totalPE': ('totalPE',),
               'totalME': ('totalME',),
               'LM': ('LMx', 'LMy', 'LMz'),
               'AM': ('AMx', 'AMy', 'AMz'),
               'absLM': ('absLM',),
               'absAM': ('absAM',),
               'totalLM': ('totalLMx', 'totalLMy', 'totalLMz'),
               'totalAM': ('totalAMx', 'totalAMy', 'totalAMz'),
               'absTotalLM': ('absTotalLM',),
               'absTotalAM': ('absTotalAM',),
               'totalAbsLM': ('totalAbsLM',),
               'totalAbsAM': ('totalAbsAM',)}

    #goes through settings
    for key, value in settings['plot'].items():
        #if the value is just a lone boolean, it is converted to a 1 element
        #list to simplify later filtering.
        if not isinstance(value, list):
            value = [value]
        #iterates through the settings values (the boolean values)
        for i, v in enumerate(value):
            #if v is false, then the data needs to be removed.
            if not v:
                #goes through the dictionary and gets all keys that match the
                #reference key from refDict
                dataKeys =  [k for k in dictionary.keys()
                             if k.startswith(refDict[key][i])]
                #then goes through the dictionary keys and removes each example
                for k in dataKeys:
                    dictionary.pop(k, None)
    return dictionary


if __name__ == "__main__":
    file = os.path.dirname(__file__)
    path = os.path.join(file,"plotconfig.toml")

    with open(path, 'rb') as f:
        settings = toml.load(f)

    settings = checkSettings(settings)

    filename = settings['filename']

    picklepath = os.path.join(file, f"{filename}.pickle")

    with open(picklepath, 'rb') as file:
        data = pickle.load(file)

    timedelta = data[0]
    data = data[1:]
    steps = len(data)
    t = np.arange(0, steps*timedelta, timedelta)
    dictionary = plotValues(data)
    dictionary = filterDiagnosticData(dictionary, settings)
    for key, value in dictionary.items():
        plt.plot(t, value, label=key)
    plt.legend()
    plt.show()
