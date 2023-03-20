import numpy as np
import matplotlib.pyplot as plt
from Pendulum import Pendulum
from Chain import Chain

if __name__ == '__main__':
    p1 = Pendulum(length=10, pos=[np.pi/2, 0], vel=[0,0], mass=1.5)
    p2 = Pendulum(length=10, pos=[2.8, 0.1], vel=[0,0], mass=10)
    p3 = Pendulum(length=10, pos=[np.pi, 0.1], vel=[0,0], mass=5)
    Chain1 = Chain(pList=[p1, p2, p3], g=9.81)
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
    KE = []
    PE = []
    Tot = []
    x3 = []
    y3 = []
    z3 = []
    for i in range(800):

        #whichPole1.append((lambda x: 10 if x else -10)(Chain1.pList[0].zPolar))
        #whichPole2.append((lambda x: 10 if x else -10)(Chain1.pList[1].zPolar))
        #whichPole3.append((lambda x: 10 if x else -10)(Chain1.pList[2].zPolar))
        p1 = Chain1.pList[0]
        pos1 = p1.toCartesian()
        p2 = Chain1.pList[1]
        pos2 = p2.toCartesian()
        p3 = Chain1.pList[2]
        pos3 = p3.toCartesian()

        x1.append(pos1[0])
        y1.append(pos1[1])
        z1.append(pos1[2])

        x2.append(pos2[0])
        y2.append(pos2[1])
        z2.append(pos2[2])

        x3.append(pos3[0])
        y3.append(pos3[1])
        z3.append(pos3[2])

        KE.append(Chain1.KE())
        PE.append(Chain1.PE())

        Tot.append(Chain1.KE()+Chain1.PE())

        #vel = p1.vel[0]*p1.dT() + p1.vel[1]*p1.dP()
        #v = (1/q.length)**2 *vel @ q.dT()
        #vx.append(vel[0])
        #vy.append(vel[1])
        #vz.append(vel[2])
        #vels.append(q.vel)
        #poss.append(q.pos)
        #theta.append(Chain1.pList[0].pos[0])
        #phi.append(Chain1.pList[0].pos[1])
        #thetadot.append(Chain1.pList[0].vel[0])
        #phidot.append(Chain1.pList[0].vel[1])
        #theta2.append(Matr.qList[1].pos[0])
        #phi2.append(Matr.qList[1].pos[1])
        #theta2dot.append(Matr.qList[1].vel[0])
        #phi2dot.append(Matr.qList[1].vel[1])
        Chain1.fixCoords()
        #Matr = Matr + 0.1**3 * Mat._Func(Matr)
        #Matr.RK(0.1**3, [1], np.array([]), Mat._Func)
        Chain1.RK4(timeDelta=0.1**2)
        Chain1.fixCoords()

    t = np.arange(0, 800*0.1**2, 0.1**2)

    plt.plot(t, x1, label="x1")
    plt.plot(t, y1, label="y1")
    plt.plot(t, z1, label="z1")

    plt.plot(t, x2, label="x2")
    plt.plot(t, y2, label="y2")
    plt.plot(t, z2, label="z2")

    plt.plot(t, x3, label="x3")
    plt.plot(t, y3, label="y3")
    plt.plot(t, z3, label="z3")

    plt.plot(t, KE, label="Kinetic")
    plt.plot(t, PE, label="Potential")
    plt.plot(t, Tot, label="Total")
    #plt.plot(t, vx, label="vx")
    #plt.plot(t, vy, label="vy")
    #plt.plot(t, vz, label="vz")
    #plt.plot(t, whichPole1)
    #plt.plot(t, whichPole2)
    #plt.plot(t, theta, label="t1")
    #plt.plot(t, phi, label="p1")
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

