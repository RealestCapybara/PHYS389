import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.animation as ani

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
        #self.fixCoord()
        #other.fixCoord()
        #if self.zPolar != other.zPolar:
        #    other.fixCoord(force=True)

        Q = Pendulum(pos=self.pos+other.pos,vel=self.vel+other.vel,
                        acc=self.acc+other.acc,mass=self.mass,
                        length=self.length,zPolar=self.zPolar)

        #Q.fixCoord()
        return Q

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

    def fixCoord(self, release=True):

        theta = deepcopy(self.pos[0])
        phi = deepcopy(self.pos[1])

        phi = np.pi * ((phi/np.pi)%2)
        theta = np.pi * ((theta/np.pi)%2)

        if theta > np.pi:
            theta = 2*np.pi - theta
            phi += np.pi
            phi = np.pi * ((phi/np.pi)%2)

        if abs(np.sin(theta)) > 0.5:
            self.pos = np.array([theta, phi])
            return False

        cosTheta = np.cos(theta)
        sinTheta = np.sin(theta)
        cosPhi = np.cos(phi)
        sinPhi = np.sin(phi)

        thetadot = deepcopy(self.vel[0])
        phidot = deepcopy(self.vel[1])
        oldDT = deepcopy(self.dT())
        oldDP = deepcopy(self.dP())
        v = thetadot * oldDT + phidot * oldDP

        if self.zPolar:
            phi2 = np.arctan2(-cosTheta, sinTheta * sinPhi)
            theta2 = np.arctan2(np.sqrt(sinTheta**2*sinPhi**2+cosTheta**2),
                                sinTheta*cosPhi)

            #thetadot2 = (phidot*sinTheta*sinPhi -
            #             thetadot*cosTheta*cosPhi)/np.sin(theta2)

            #phidot2 = -thetadot2 + thetadot*sinTheta/(np.sin(phi2)*\
            #                                          np.cos(theta2))
        else:
            phi2 = np.arctan2(sinTheta * cosPhi, cosTheta)
            theta2 = np.arctan2(np.sqrt(cosTheta**2 + sinTheta**2 * cosPhi**2),
                                -sinTheta*sinPhi)

            #thetadot2 = (phidot+thetadot)*sinPhi*cosTheta/np.sin(theta2)

            #phidot2 = (thetadot2*np.cos(theta2)*np.cos(phi2) +\
            #           thetadot*sinTheta)/(np.sin(theta2)*np.sin(phi2))

        phi2 = np.pi * ((phi2/np.pi)%2)
        theta2 = np.pi * ((theta2/np.pi)%2)

        self.pos[0] = theta2
        self.pos[1] = phi2

        self.zPolar = not self.zPolar

        #self.vel[0] = thetadot2
        #self.vel[1] = phidot2
        self.vel[0] = (1/self.length)**2 * v @ self.dT()
        self.vel[1] = (1/(self.length * np.sin(theta2)))**2 * v@self.dP()

        if release:
            print(v)
            print(self.vel[0]*self.dT()+self.vel[1]*self.dP())
            print(v-self.vel[0]*self.dT()-self.vel[1]*self.dP())
            print(np.linalg.norm(v)-np.linalg.norm(self.vel[0]*self.dT()+
                                                   self.vel[1]*self.dP()))

        return True and release

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

    def __str__(self):
        return ",\n         ".join([f"\
q{i}: pos=({v.pos[0]/np.pi:.4f} pi, {v.pos[1]/np.pi:.4f} pi) \
vel=({v.vel[0]/np.pi:.4f} pi, {v.vel[1]/np.pi:.4f} pi) \
pole={(lambda zPole: 'z' if zPole else 'x')(v.zPolar)}"
                                    for i, v in enumerate(self.qList)])

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

        M = Mat(qList = newQList, g=g)
        #M.fixCoords(release=False)
        return M

    def RK(self, timeDelta, weights, RKmatrix, function):
        '''try:
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

        if RKmatrix.ndim != 2 and len(RKmatrix) !=0:
            raise ValueError("RKmatrix must be a matrix")

        if len(RKmatrix[-1]) != (len(weights)-1):
            raise ValueError(
                "length of weights must be one more than last row of RKmatrix")

        if len(RKmatrix) != len(RKmatrix[-1]):
            raise ValueError(
                "number of rows and cols of RKmatrix must be equal")

        if round(sum(weights), 10) != 1.0:
            raise ValueError("sum of weights must be 1")'''

        k = [function(self)]
        x = deepcopy(self)
        for i, v in enumerate(RKmatrix):
            delta = sum([j*a for j, a in zip(k, v[:i+1])]) * timeDelta
            k.append(function(x + delta))
        d = timeDelta * sum([b * j for b, j in zip(weights, k)])
        self.qList = (x + d).qList

    def RK4(self, timeDelta):
        RKmatrix = np.array([[0.5,0,0], [0,0.5,0], [0,0,1]])
        weights = np.array([1/6, 1/3, 1/3, 1/6])
        self.RK(timeDelta, weights, RKmatrix, self._Func)

    def fixCoords(self, release=True):
        qList = deepcopy(self.qList)
        val = True
        for q in qList:
            vi = q.fixCoord(release)
            val = val and vi
        self.qList = qList
        return val

if __name__ == '__main__':
    q1 = Pendulum(length=10, pos=[np.pi/2, 0], vel=[-np.pi,2.8887156], acc=[0,0], mass=1)
    q2 = Pendulum(length=10, pos=[0.1, 0.1], vel=[0,0], acc=[0,0], mass=1)
    q3 = Pendulum(length=10, pos=[np.pi, 0.1], vel=[0,0], acc=[0,0], mass=1)
    Matr = Mat(qList=[q1], g=9.81)
    x1 = []
    y1 = []
    z1 = []
    x2 = []
    y2 = []
    z2 = []
    theta = []
    phi = []
    thetadot = []
    phidot = []
    theta2 = []
    phi2 = []
    theta2dot = []
    phi2dot = []
    whichPole1 = []
    whichPole2 = []
    vx = []
    vy = []
    vz = []
    vels = []
    poss = []
    for i in range(1000):

        whichPole1.append((lambda x: 10 if x else -10)(Matr.qList[0].zPolar))
        #whichPole2.append((lambda x: 10 if x else -10)(Matr.qList[1].zPolar))
        q1 = Matr.qList[0]
        pos1 = q1.toCartesian()
        #q2 = Matr.qList[1]
        #pos2 = q2.toCartesian()

        x1.append(pos1[0])
        y1.append(pos1[1])
        z1.append(pos1[2])

        #x2.append(pos2[0])
        #y2.append(pos2[1])
        #z2.append(pos2[2])

        #vel = q.vel[0]*q.dT() + q.vel[1]*q.dP()
        #v = (1/q.length)**2 *vel @ q.dT()
        #vx.append(vel[0])
        #vy.append(vel[1])
        #vz.append(vel[2])
        #vels.append(q.vel)
        #poss.append(q.pos)
        theta.append(Matr.qList[0].pos[0])
        phi.append(Matr.qList[0].pos[1])
        thetadot.append(Matr.qList[0].vel[0])
        phidot.append(Matr.qList[0].vel[1])
        #theta2.append(Matr.qList[1].pos[0])
        #phi2.append(Matr.qList[1].pos[1])
        #theta2dot.append(Matr.qList[1].vel[0])
        #phi2dot.append(Matr.qList[1].vel[1])
        Matr.fixCoords()
        #Matr = Matr + 0.1**3 * Mat._Func(Matr)
        #Matr.RK(0.1**3, [1], np.array([]), Mat._Func)
        Matr.RK4(timeDelta=0.1**2)
        Matr.fixCoords()

    t = np.arange(0, 1000*0.1**2, 0.1**2)

    plt.plot(t, x1, label="x1")
    plt.plot(t, y1, label="y1")
    plt.plot(t, z1, label="z1")

    #plt.plot(t, x2, label="x2")
    #plt.plot(t, y2, label="y2")
    #plt.plot(t, z2, label="z2")

    #plt.plot(t, vx, label="vx")
    #plt.plot(t, vy, label="vy")
    #plt.plot(t, vz, label="vz")
    #plt.plot(t, whichPole1)
    #plt.plot(t, whichPole2)
    #plt.plot(t, theta, label="t1")
    #plt.plot(t, phi, label="pi")
    #plt.plot(t, thetadot, label="t1d")
    #plt.plot(t, phidot, label="p1d")
    #plt.plot(t, theta2, label="t2")
    #plt.plot(t, phi2, label="p2")
    #plt.plot(t, theta2dot, label="t2d")
    #plt.plot(t, phi2dot, label="p2d")


    #plt.plot(t[1:], [v[0]*0.1**3 + poss[i-1][0] for i, v in enumerate(vels[1:])], label="euler Position")
    #plt.plot(t, [i[0]*i[1]*(5/np.pi)**2 for i in vels], label="phidot * thetadot")
    #plt.plot(t, [(i[0]*5/np.pi)**2 for i in vels], label="thetadot ** 2")
    #plt.plot(t, [(i[1]*5/np.pi)**2 for i in vels], label="phidot ** 2")
    plt.legend()
    plt.show()

