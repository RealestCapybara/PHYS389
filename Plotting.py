import numpy as np
from Chain import Chain
from Pendulum import Pendulum
import matplotlib.plt as plt
import pickle
import tomllib as toml

filename="test"

with open(f"{filename}.pickle", 'rb') as file:
    data = pickle.load(file)

timestep = data[0]

data = data[1:]

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
    States = []
    for i, v in enumerate(data):
        pList = v.pList
        for key, value in settings["plot"].items():
            if (lambda x: any(x) if isinstance(x, list) else x)(value):
                if key == 'x':
                    newVal = [p.toCartesian() for p in pList]
                    newVal = [[p[i] if val else None
                               for i, val in enumerate(value)] for p in newVal]
                    newVal = np.array(newVal).transpose()

                if key == 'v':
                    newVal = [p.Velocity() for p in pList]
                    newVal = [[p[i] if val else None
                               for i, val in enumerate(value)] for p in newVal]
                    newVal = np.array(newVal).transpose()

                if key == 'pos':
                    newVal = [p.pos for p in pList]
                    newVal = [[p[i] if val else None
                               for i, val in enumerate(value)] for p in newVal]
                    newVal = np.array(newVal).transpose()

                if key == 'vel':
                    newVal = [p.vel for p in pList]
                    newVal = [[p[i] if val else None
                               for i, val in enumerate(value)] for p in newVal]
                    newVal = np.array(newVal).transpose()

                if key == 'pole':
                    newVal = [p.zPolar if value else None for p in pList]
                if key == 'KE':
                    newVal = (lambda sys: sys.KE() if value else None)(v)
                if key == 'PE':
                    newVal = (lambda sys: sys.PE() if value else None)(v)
                if key == 'totalE':
                    newVal = (lambda sys: sys.PE() + sys.KE()
                              if value else None)(v)

        if settings['plot']['KE']:
            T = v.KE()
        if settings[]
        V = v.PE()
        E = T + P
        xs = []
        ys = []
        zs = []
        vxs = []
        vys = []
        vzs = []
        thetas = []
        phis = []
        thetadots = []
        phidots = []
        poles = []
        for p in pList:
            x, y, z = p.toCartesian()
            vx, vy, vz = p.Velocity()
            theta = p.pos[0]
            phi = p.pos[1]
            thetadot = p.vel[0]
            phidot = p.vel[1]
            whichPole = (lambda x: 10 if x.zPolar==True else -10)(p)
            xs.append(x)
            ys.append(y)
            zs.append(z)
            vxs.append(vx)
            vys.append(vy)
            vzs.append(vz)
            thetas.append(theta)
            phis.append(phi)
            thetadots.append(thetadot)
            phidots.append(phidot)
            poles.append(whichPole)


