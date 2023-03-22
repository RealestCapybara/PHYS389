import numpy as np
from Chain import Chain
from Pendulum import Pendulum
import matplotlib.pyplot as plt
import pickle
import tomllib as toml


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

def plotValues(data, settings):
    """Returns a dict with the values to plot."""

    #THIS ALL NEEDS WORK - VERY CLUNKY AND DIFFICULT TO UNDERSTAND
    #ensure unit testing is finished before coming back to improve.
    FrameDict = {}
    #iterates through each frame and collects v, the system element and i, the frame
    for i, v in enumerate(data):
        #gets the list of pendula in the given frame
        pList = v.pList
        #iterates through the settings to check each value
        for key, value in settings["plot"].items():
            #checks if setting have a given diagnostic value listed as true.
            if (lambda x: any(x) if isinstance(x, list) else x)(value):
                #holds the relevant generation functions and the final dictionary keys in a dict associated with the respective settings key
                #this is done to avoid repetetive code
                keyfuncs1 = {
                          'x': [(lambda p: p.toCartesian()), ('x', 'y', 'z')],
                          'v': [(lambda p: p.Velocity()), ('vx', 'vy', 'vz')],
                          'pos': [(lambda p: p.pos), ('theta', 'phi')],
                          'vel': [(lambda p: p.vel), ('thetadot', 'phidot')]}
                #if the 'true' key is in the list of keys, then this code can be executed
                if key in keyfuncs1.keys():
                    #creates a 2D list containing the final values for each pendulum (row)
                    newVal = [keyfuncs1[key][0](p) for p in pList]
                    #replaces the values with a None value if the specific coordinate is specified as false in settings
                    newVal = [[p[i] if val else None
                               for i, val in enumerate(value)] for p in newVal]
                    #puts the values into a dictionary element with a key for each combiantion of desired value and pendulum.
                    dictElement = {f"{coord}{pi}":[newVal[pi][ci]]
                                   for pi in range(len(newVal))
                                   for ci, coord in enumerate(
                                       keyfuncs1[key][1])}
                #another value. This key acts differently to all others, so is treated seperately.
                if key == 'pole':
                    #checks what coordinate basis the pendula are using, and if they are in the z-basis, puts a 10, -10 otherwise.
                    newVal = [(lambda b: 10 if b else -10)(p.zPolar)
                              for p in pList]
                    #creates a dictionary element, with a key for each pendulum
                    dictElement = {f"whichPole{pi}":[v] for pi, v, in enumerate(newVal)}
                #stores three functions in a dictionary. This is done to avoid repetetive code.
                keyfuncs2 = {'KE': (lambda sys: sys.KE()),
                             'PE': (lambda sys: sys.PE()),
                             'totalE': (lambda sys: sys.KE()+sys.PE())}
                #checks if the key is in this second set of grouped functions
                if key in keyfuncs2.keys():
                    #as these are functions applied to the entire system, only one value needs to be stored.
                    newVal = keyfuncs2[key](v)
                    dictElement = {key: [newVal]}

                #if the keys of the dictionary element are a subset of the desired dictionary for the entire simulation, then the keys have already been added.
                if not set(dictElement.keys()) <= set(FrameDict.keys()):
                    #if the keys dont exist in the main dictionary, it is simply updated with these new values.
                    FrameDict.update(dictElement)
                else:
                    #otherwise, the new elements are simply added onto the existing values, such that the values for that key can be iterated through for plotting.
                    for key, value in dictElement.items():
                        FrameDict[key].extend(value)

    return FrameDict

if __name__ == "__main__":
    with open("plotconfig.toml", 'rb') as file:
        settings = toml.load(file)

    checkSettings(settings)

    settings = settings['system']
    filename = settings['filename']

    with open(f"{filename}.pickle", 'rb') as file:
        data = pickle.load(file)

    timedelta = data[0]
    data = data[1:]
    for sys in data[:10]:
        print(sys)
    steps = len(data)
    t = np.arange(0, steps*timedelta, timedelta)
    dictionary = plotValues(data, settings)

    for key, value in dictionary.items():
        plt.plot(t, value, label=key)
    plt.legend()
    plt.show()
