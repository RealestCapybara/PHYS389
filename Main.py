import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from copy import deepcopy

class Pendulum():
    def __init__(self, pos, vel, acc, mass, zPolar=True):
        try:
            self.pos = np.array(pos)
            self.vel = np.array(vel)
            self.acc = np.array(acc)
        except TypeError:
            raise TypeError("pos, vel, and acc must be array-like")

        if len(self.pos) != 3:
            raise ValueError("pos must have len == 3")
        if len(self.vel) != 2:
            raise ValueError("vel must have len == 2")
        if len(self.acc) != 2:
            raise ValueError("acc must have len == 2")

        try:
            self.mass = float(mass)
        except TypeError:
            raise TypeError("mass must be float-like")

        if self.mass <= 0:
            raise ValueError("mass must be positive and non-zero")

        self.zPolar = bool(zPolar)

    def diffTheta(self):
        l = self.pos[0]
        theta = self.pos[1]
        phi = self.pos[2]

        if self.zPolar:
            return np.array([l*np.cos(theta)*np.cos(phi),
                             l*np.cos(theta)*np.sin(phi),
                             l*np.sin(theta)])
        else:
            return np.array([-l*np.sin(theta),
                             l*np.cos(theta)*np.cos(phi),
                             l*np.cos(theta)*np.sin(phi)])

    def diffPhi(self):
        l = self.pos[0]
        theta = self.pos[1]
        phi = self.pos[2]

        if self.zPolar:
            return np.array([-l*np.sin(theta)*np.sin(phi),
                             l*np.sin(theta)*np.cos(phi),
                             0])
        else:
            return np.array([0,
                             -l*np.sin(theta)*np.sin(phi),
                             l*np.sin(theta)*np.cos(phi)])
    def diffTheta2(self):
        l = self.pos[0]
        theta = self.pos[1]
        phi = self.pos[2]

        if self.zPolar:
            return np.array([-l*np.sin(theta)*np.cos(phi),
                             -l*np.sin(theta)*np.sin(phi),
                             l*np.cos(theta)])
        else:
            return np.array([-l*np.cos(theta),
                             -l*np.sin(theta)*np.cos(phi),
                             -l*np.sin(theta)*np.sin(phi)])

    def diffPhi2(self):
        l = self.pos[0]
        theta = self.pos[1]
        phi = self.pos[2]

        if self.zPolar:
            return np.array([-l*np.sin(theta)*np.cos(phi),
                             -l*np.sin(theta)*np.sin(phi),
                             0])
        else:
            return np.array([0,
                             -l*np.sin(theta)*np.cos(phi),
                             -l*np.sin(theta)*np.sin(phi)])

    def diffThetaPhi(self):
        l = self.pos[0]
        theta = self.pos[1]
        phi = self.pos[2]

        if self.zPolar:
            return np.array([-l*np.cos(theta)*np.sin(phi),
                             l*np.cos(theta)*np.cos(phi),
                             0])
        else:
            return np.array([0,
                             -l*np.cos(theta)*np.sin(phi),
                             l*np.cos(theta)*np.cos(phi)])

    def fixCoord(self):
        if self.pos[1] >= (np.pi/4):
            return None

        theta = deepcopy(self.pos[1])
        phi = deepcopy(self.pos[2])

        cosTheta = np.cos(theta)
        sinTheta = np.sin(theta)
        cosPhi = np.cos(phi)
        sinPhi = np.sin(phi)

        thetadot = deepcopy(self.vel[0])
        phidot = deepcopy(self.vel[1])
        v = thetadot * self.diffTheta() + phidot * self.diffPhi()

        if self.zPolar:
            phi2 = np.arctan2(-cosTheta, sinTheta * sinPhi)
            theta2 = np.arctan2(np.sqrt(sinTheta**2 * sinPhi**2 + cosTheta**2),
                                sinTheta * cosPhi)
        else:
            phi2 = np.arctan2(sinTheta*cosPhi, cosTheta)
            theta2 = np.arctan2(np.sqrt(cosTheta**2 + sinTheta**2 * cosPhi**2),
                                -sinTheta * sinPhi)

        self.pos[1] = theta2
        self.pos[2] = phi2

        self.vel[0] = (1/self.pos[0])**2 * v * self.diffTheta
        self.vel[1] = (1/(self.pos[0] * np.sin(theta2)))**2 *v*self.diffPhi()

        self.zPolar = not self.zPolar

    def toCartesian(self):
        l = self.pos[0]
        theta = self.pos[1]
        phi = self.pos[2]
        if self.zPolar():
            return np.array([l*np.sin(theta)*np.cos(phi),
                             l*np.sin(theta)*np.sin(phi),
                             -l*np.cos(theta)])
        else:
            return np.array([l*np.cos(theta),
                             l*np.sin(theta)*np.cos(phi),
                             l*np.sin(theta),np.sin(phi)])


