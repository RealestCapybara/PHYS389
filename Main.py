import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from mpl_tookkits import Axes3D

class Pendulum():
    def __init__(self, length, pos, vel, acc, mass, zPolar=True):
        try:
            self.pos = np.array(pos)
            self.vel = np.array(vel)
            self.acc = np.array(acc)
        except TypeError:
            raise TypeError("pos, vel, and acc must be array-like")

        if len(self.pos) != 2:
            raise ValueError("pos must have len == 2")
        if len(self.vel) != 2:
            raise ValueError("vel must have len == 2")
        if len(self.acc) != 2:
            raise ValueError("acc must have len == 2")

        try:
            self.mass = float(mass)
            self.length = float(length)
        except TypeError:
            raise TypeError("mass and length must be float-like")

        if self.mass <= 0 or self.length <=0:
            raise ValueError("mass and length must be positive and non-zero")

        self.zPolar = bool(zPolar)

    def __mul__(self, other):
        try:
            other = np.float32(other)
        except TypeError:
            raise TypeError("Pendulum can only be multiplied by a scalar")
        return Pendulum(pos=self.pos*other,vel=self.vel*other,
                        acc=self.acc*other,mass=self.mass,length=self.length,
                        zPolar=self.zPolar)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        if type(other) is not Pendulum:
            raise TypeError(
"Pendulum object can only be added to other Pendulum objects")
        if self.zPolar != other.zPolar:
            other.fixCoord(force=True)
        return Pendulum(pos=self.pos+other.pos,vel=self.vel+other.vel,
                        acc=self.acc+other.acc,mass=self.mass,
                        length=self.length,zPolar=self.zPolar)

    def dT(self):
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return np.array([l*np.cos(theta)*np.cos(phi),
                             l*np.cos(theta)*np.sin(phi),
                             l*np.sin(theta)])
        else:
            return np.array([-l*np.sin(theta),
                             l*np.cos(theta)*np.cos(phi),
                             l*np.cos(theta)*np.sin(phi)])

    def dP(self):
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return np.array([-l*np.sin(theta)*np.sin(phi),
                             l*np.sin(theta)*np.cos(phi),
                             0])
        else:
            return np.array([0,
                             -l*np.sin(theta)*np.sin(phi),
                             l*np.sin(theta)*np.cos(phi)])
    def dT2(self):
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return np.array([-l*np.sin(theta)*np.cos(phi),
                             -l*np.sin(theta)*np.sin(phi),
                             l*np.cos(theta)])
        else:
            return np.array([-l*np.cos(theta),
                             -l*np.sin(theta)*np.cos(phi),
                             -l*np.sin(theta)*np.sin(phi)])

    def dP2(self):
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return np.array([-l*np.sin(theta)*np.cos(phi),
                             -l*np.sin(theta)*np.sin(phi),
                             0])
        else:
            return np.array([0,
                             -l*np.sin(theta)*np.cos(phi),
                             -l*np.sin(theta)*np.sin(phi)])

    def dTP(self):
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return np.array([-l*np.cos(theta)*np.sin(phi),
                             l*np.cos(theta)*np.cos(phi),
                             0])
        else:
            return np.array([0,
                             -l*np.cos(theta)*np.sin(phi),
                             l*np.cos(theta)*np.cos(phi)])

    def drdT(self):
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return l * np.sin(theta)
        else:
            return l * np.cos(theta) * np.sin(phi)

    def drdP(self):
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return 0
        else:
            return l * np.sin(theta) * np.cos(phi)

    def fixCoord(self, force=False):
        if self.pos[0] >= (np.pi/4) and not force:
            return False

        theta = deepcopy(self.pos[0])
        phi = deepcopy(self.pos[1])

        cosTheta = np.cos(theta)
        sinTheta = np.sin(theta)
        cosPhi = np.cos(phi)
        sinPhi = np.sin(phi)

        thetadot = deepcopy(self.vel[0])
        phidot = deepcopy(self.vel[1])
        v = thetadot * self.dT() + phidot * self.dP()

        if self.zPolar:
            phi2 = np.arctan2(-cosTheta, sinTheta * sinPhi)
            theta2 = np.arctan2(np.sqrt(sinTheta**2 * sinPhi**2 + cosTheta**2),
                                sinTheta * cosPhi)
        else:
            phi2 = np.arctan2(sinTheta*cosPhi, cosTheta)
            theta2 = np.arctan2(np.sqrt(cosTheta**2 + sinTheta**2 * cosPhi**2),
                                -sinTheta * sinPhi)

        self.pos[0] = theta2
        self.pos[1] = phi2

        self.vel[0] = (1/self.length)**2 * v @ self.dT()
        self.vel[1] = (1/(self.length * np.sin(theta2)))**2 * v@self.dP()

        self.zPolar = not self.zPolar
        return True

    def toCartesian(self):
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]
        if self.zPolar:
            return np.array([l*np.sin(theta)*np.cos(phi),
                             l*np.sin(theta)*np.sin(phi),
                             -l*np.cos(theta)])
        else:
            return np.array([l*np.cos(theta),
                             l*np.sin(theta)*np.cos(phi),
                             l*np.sin(theta)*np.sin(phi)])

class Mat():
    def __init__(self, qList, g):
        self.qList = qList
        self.g = g

    def __mul__(self, other):
        try:
            other = np.float32(other)
        except TypeError:
            raise TypeError("Mat can only be multiplied with float-like types")
        return Mat(qList = [q*other for q in self.qList], g=self.g)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __add__(self, other):
        if type(other) is Mat:
            return Mat(qList = [q + p for q, p in
                                zip(self.qList,other.qList)], g = self.g)
        elif not other:
            return self
        else:
            raise TypeError("Mat can only be added to other Mat objects")
    def __radd__(self, other):
        return self.__add__(other)

    @staticmethod
    def _AElement(q1, q2):
        if type(q1) is not Pendulum or type(q2) is not Pendulum:
            raise TypeError("arguments must be of type Pendulum")
        return np.array([[q2.dT()@q1.dT(), q2.dP()@q1.dT()],
                         [q2.dT()@q1.dP(), q2.dP()@q1.dP()]])

    @staticmethod
    def _BElement(q1, q2):
        if type(q1) is not Pendulum or type(q2) is not Pendulum:
            raise TypeError("arguments must be of type Pendulum")
        return np.array([[q2.dT2()@q1.dT(),q2.dTP()@q1.dT(),q2.dP2()@q1.dT()],
                         [q2.dT2()@q1.dP(),q2.dTP()@q1.dP(),q2.dP2()@q1.dP()]])

    def Matrices(self):
        qList = self.qList
        MassList = [q.mass for q in qList]
        arrayA = []
        arrayB = []
        for i, vi in enumerate(qList):
            rowA = []
            rowB = []
            for j, vj in enumerate(qList):
                M = sum(MassList[max(i, j):])
                ElementA = self._AElement(vi, vj) * M
                ElementB = self._BElement(vi, vj) * M
                rowA.append(ElementA)
                rowB.append(ElementB)
            arrayA.append(rowA)
            arrayB.append(rowB)
        return np.block(arrayA), np.block(arrayB)

    def GVec(self):
        g = self.g
        qList = self.qList
        array = []
        massList = [q.mass for q in qList]
        for i, q in enumerate(qList):
            M = sum(massList[i:])
            array.append(g*M*q.drdT())
            array.append(g*M*q.drdP())
        return np.array(array)

    @staticmethod
    def _Func(Matr):
        if type(Matr) is not Mat:
            raise TypeError("argument must be of type Mat")
        qList = Matr.qList
        g = Matr.g

        A, B = Matr.Matrices()
        V3 = Matr.GVec()

        V2 = [f(q) for q in qList for f in (lambda q: q.vel[0]**2,
                                            lambda q: q.vel[0]*q.vel[1],
                                            lambda q: q.vel[1]**2)]

        V2 = np.array(V2)

        V1 = -np.linalg.inv(A) @ (B@V2 + V3)

        omegaThetaDot = [v for i, v in enumerate(V1) if i%2 == 0]
        omegaPhiDot = [v for i, v in enumerate(V1) if i%2 == 1]

        ThetaDot = [q.vel[0] for q in qList]
        PhiDot = [q.vel[1] for q in qList]

        newQList = [Pendulum(pos=[ThetaDot[i], PhiDot[i]],
                             vel=[omegaThetaDot[i], omegaPhiDot[i]],
                             acc=[0, 0], length=q.length,
                             mass=q.mass) for i, q in enumerate(qList)]

        return Mat(qList = newQList, g=g)

    def RK(self, timeDelta, weights, RKmatrix, function):
        try:
            weights = np.array(weights)
            RKmatrix = np.array(RKmatrix)
        except TypeError:
            raise TypeError("weights and RKmatrix must be array-like")

        try:
            timeDelta = float(timeDelta)
        except TypeError:
            raise TypeError("timeDelta must be float-like")

        if not callable(function):
            raise TypeError("function must be a function")

        if RKmatrix.ndim != 2:
            raise ValueError("RKmatrix must be a matrix")

        if len(RKmatrix[-1]) != (len(weights)-1):
            raise ValueError(
                "length of weights must be one more than last row of RKmatrix")

        if len(RKmatrix) != len(RKmatrix[-1]):
            raise ValueError(
                "number of rows and cols of RKmatrix must be equal")

        if round(sum(weights), 10) != 1.0:
            raise ValueError("sum of weights must be 1")

        k = [function(self)]

        for i, v in enumerate(RKmatrix):
            delta = sum([j*a for j, a in zip(k, v[:i+1])]) * timeDelta
            k.append(function(self + delta))

        self += timeDelta * sum([b * j for b, j in zip(weights, k)])

    def RK4(self, timeDelta):
        RKmatrix = np.array([[0.5,0,0], [0,0.5,0], [0,0,1]])
        weights = np.array([1/6, 1/3, 1/3, 1/6])
        self.RK(timeDelta, weights, RKmatrix, self._Func)

    def fixCoords(self):
        qList = self.qList
        for q in qList:
            q.fixCoord()
        self.qList = qList

if __name__ == '__main__':
    q1 = Pendulum(length=10, pos=[0, 0], vel=[0,0], acc=[0,0], mass=1)
    #q2 = Pendulum(length=10, pos=[0, 0], vel=[0,0], acc=[0,0], mass=1)
    q2 = Pendulum(length=10, pos=[np.pi, 0], vel=[0,0], acc=[0,0], mass=1)
    Matr = Mat(qList=[q1,q2], g=9.81)
    for i in range(600):
        Matr.fixCoords()
        Matr.RK4(timeDelta=0.1)
        print(", ".join([f"q{i}: {v.toCartesian()}" for i, v in
                         enumerate(Matr.qList)]))


