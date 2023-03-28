import numpy as np
from Pendulum import Pendulum
from Chain import Chain
import tomli as toml
import pickle
from tqdm import tqdm
from copy import deepcopy

def genPList(data):
    try:
        data = data['system']
    except KeyError:
        raise KeyError("config file has misnamed or absent 'system' table")
    try:
        PDatList = data['pendula']
    except KeyError:
        raise KeyError("config file has missing or misnamed 'pendula' array")

    if not isinstance(PDatList, list):
        raise ValueError("'pendula' item in config file must be an array")

    if len(PDatList) == 0:
        raise ValueError("config file must have at least one pendula item")

    pList = []

    for i, v in enumerate(PDatList):
        try:
            length = v['length']
        except KeyError:
            length = 10
            raise Warning(f"\
length for pendula {i} has not been specified, defaulting to length=10m")
        try:
            mass = v['mass']
        except KeyError:
            mass = 1
            raise Warning(f"\
mass for pendula {i} has not been specified, defaulting to mass=1kg")
        try:
            theta = v["theta"]
        except KeyError:
            theta = [0, 'r']
        try:
            phi = v["phi"]
        except KeyError:
            phi = [0, 'r']
        try:
            thetadot = v["thetadot"]
        except KeyError:
            thetadot = [0, 'r']
        try:
            phidot = v["phidot"]
        except KeyError:
            phidot = [0, 'r']

        vals = [theta, phi, thetadot, phidot]

        if not all(isinstance(val, list) for val in vals):
            raise TypeError("theta, phi, thetadot, and phidot must be arrays")

        if not all(len(val)==2 for val in vals):
            raise ValueError(
                "theta, phi, thetadot, and phidot must have len=2")

        if not all(isinstance(val[0], (float, int)) for val in vals):
            raise TypeError(
                "Initial values of theta, phi, thetadot, and phidot must be \
                of type int or float")
        if not all(isinstance(val[1], str) for val in vals):
            raise TypeError(
                "Second values of theta, phi, thetadot, and phidot must be \
                of type int")
        if not all((val[1] in ['r', 'p', 'd']) for val in vals):
            raise ValueError(
                "Second values of theta, phi, thetadot, and phidot must be \
                'r', 'p', or 'd' as these represent the accepted angle units.")

        vals = [(lambda v: np.pi*v[0] if v[1]=='p' else
                 (v[0]*np.pi/180 if v[1]=='d' else v[0]))(val) for val in vals]

        try:
            zPolar = v["zPole"]
        except KeyError:
            zPolar = True

        if not isinstance(zPolar, bool):
            raise TypeError("zPole must be a boolean value")

        p = Pendulum(length=length, mass=mass, pos=vals[:2], vel=vals[2:],
                     zPolar=zPolar)

        pList.append(p)

    return(pList)

def SystemValues(data):
    try:
        data = data['system']
    except KeyError:
        raise KeyError("config file has misnamed or absent 'system' table")

    try:
        g = data["g"]
    except KeyError:
        g = 9.81
        raise Warning("g not specified, defaulting to 9.81")

    try:
        timedelta = data["timedelta"]
    except KeyError:
        raise KeyError("'timedelta' value has been misnamed or is absent")

    try:
        steps = data["steps"]
    except KeyError:
        raise KeyError("'steps' value has been misnamed or is absent")

    try:
        filename = data["filename"]
    except KeyError:
        raise KeyError("'filename' value has been misnamed or is absent")

    if not isinstance(g, (float, int)):
        raise TypeError("'g' must be a float or int")

    if not isinstance(timedelta, (float, int)):
        raise TypeError("'timedelta' must be a float or int")

    if not isinstance(steps, int):
        raise TypeError("'steps' must be an int")

    if not isinstance(filename, str):
        raise TypeError("'filename' must be a str")

    return g, timedelta, steps, filename


if __name__ == '__main__':

    with open("config.toml", "rb") as f:
        data = toml.load(f)

    pList = genPList(data)

    g, timedelta, steps, filename = SystemValues(data)

    System = Chain(pList=pList, g=g)
    System.fixCoords()
    SystemEveryTick = [timedelta]

    for i in tqdm(range(steps), desc="Simulating..."):
        SystemEveryTick.append(deepcopy(System))
        System.RK4(timedelta)
        System.fixCoords()
    for sys in SystemEveryTick[:10]:
        print(sys)

    print(f"Saving to {filename}.pickle...")

    with open(f"{filename}.pickle", 'wb') as file:
        pickle.dump(SystemEveryTick, file, protocol=pickle.HIGHEST_PROTOCOL)
