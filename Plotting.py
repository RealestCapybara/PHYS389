import numpy as np
from Chain import Chain
from Pendulum import Pendulum
import matplotlib.pyplot as plt
import pickle
import tomllib as toml
from operator import itemgetter


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
        totalE = plot["totalE"]
    except KeyError:
        raise KeyError(
            "plotconfig file has misnamed or absent 'totalE' value")
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

    return sys

def plotValues(data):
    """Returns a dict with the values to plot."""

    #THIS ALL NEEDS WORK - VERY CLUNKY AND DIFFICULT TO UNDERSTAND
    #ensure unit testing is finished before coming back to improve.

    keyList = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'theta', 'phi', 'thetadot',
               'phidot', 'whichPole', 'KE', 'PE', 'totalE']

    FrameDict = {f'{val}{i}':[] for i in range(len(data[0].pList))
                 for val in keyList[:11]}
    FrameDict.update({k:[] for k in keyList[11:]})

    #iterates through each frame and collects v, the system element and i, the
    #frame
    for i, v in enumerate(data):
        #gets the list of pendula in the given frame
        pList = v.pList
        #holds the relevant generation functions and the final dictionary 
        #keys in a dict associated with the respective settings key
        #this is done to avoid repetetive code
        keyfuncs1 = {
            'x': [(lambda p: p.toCartesian()), ('x', 'y', 'z')],
            'v': [(lambda p: p.Velocity()), ('vx', 'vy', 'vz')],
            'pos': [(lambda p: p.pos), ('theta', 'phi')],
            'vel': [(lambda p: p.vel), ('thetadot', 'phidot')]}
        #stores three functions in a dictionary. This is done to avoid 
        #repetetive code.
        keyfuncs2 = {
            'KE': (lambda sys: sys.KE()),
            'PE': (lambda sys: sys.PE()),
            'totalE': (lambda sys: sys.KE()+sys.PE())}
        #iterates through the settings to check each value
        for key in ['x', 'v', 'pos', 'vel', 'pole', 'KE', 'PE', 'totalE']:
            #if the 'true' key is in the list of keys, then this code can be 
            #executed
            if key in keyfuncs1.keys():
                #creates a 2D list containing the final values for each 
                #pendulum (row)
                newVal = [keyfuncs1[key][0](p) for p in pList]
                #puts the values into a dictionary element with a key for each
                #combiantion of desired value and pendulum.
                dictElement = {f"{coord}{pi}":newVal[pi][ci]
                               for pi in range(len(newVal))
                               for ci, coord in enumerate(
                                   keyfuncs1[key][1])}
            #another value. This key acts differently to all others, so is 
            #treated seperately.
            if key == 'pole':
                #checks what coordinate basis the pendula are using, and if
                #they are in the z-basis, puts a 10, -10 otherwise.
                newVal = [(lambda b: 10 if b else -10)(p.zPolar)
                          for p in pList]
                #creates a dictionary element, with a key for each pendulum
                dictElement = {f"whichPole{pi}":v for pi, v, in
                               enumerate(newVal)}

            #checks if the key is in this second set of grouped functions
            if key in keyfuncs2.keys():
                #as these are functions applied to the entire system, only one
                #value needs to be stored.
                newVal = keyfuncs2[key](v)
                dictElement = {key: newVal}

            #new elements are added to FrameDict
            for key, value in dictElement.items():
                #print(f"{key}:{value}")
                FrameDict[key].append(value)

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
               'totalE': ('totalE',)}
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
    with open("plotconfig.toml", 'rb') as file:
        settings = toml.load(file)

    settings = checkSettings(settings)

    filename = settings['filename']

    with open(f"{filename}.pickle", 'rb') as file:
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
