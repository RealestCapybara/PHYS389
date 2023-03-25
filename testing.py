from Pendulum import Pendulum
from Chain import Chain
import numpy as np
from numpy.testing import assert_array_almost_equal as arrayAlmostEqual
from numpy.testing import assert_array_equal as arrayEqual
import unittest
import random as r
from random import random as rand
from copy import deepcopy

class TestPendulumMethods(unittest.TestCase):

    def test_eq(self):
        p1 = Pendulum(length=10, mass=1, pos=[0.3, 0.5], vel=[0.1, 0.1])
        p2 = Pendulum(length=10, mass=1, pos=[0.3, 0.5], vel=[0.1, 0.1])
        p3 = Pendulum(length=10, mass=1, pos=[0.1+0.2, 0.5], vel=[0.1, 0.1])
        self.assertEqual(p1, p2)
        self.assertAlmostEqual(p1, p3)

    def test_add(self):
        p1 = Pendulum(length=10, mass=1, pos=[0.5, 0.5], vel=[0.1, 0.1])
        p2 = Pendulum(length=10, mass=1, pos=[0.6, 0.7], vel=[0.2, 0.3])
        p3 = Pendulum(length=10.0, mass=1.0, pos=[1.1, 1.2], vel=[0.3, 0.4])
        self.assertAlmostEqual((p1 + p2), p3)
        self.assertEqual(p1 + p2, p2 + p1)

    def test_mul(self):
        p1 = Pendulum(length=5, mass=1, pos=[0.3, 0.4], vel=[0.1, 0.2])
        p2 = Pendulum(length=5, mass=1, pos=[0.6, 0.8], vel=[0.2, 0.4])
        p3 = Pendulum(length=5, mass=1, pos=[0.9, 1.2], vel=[0.3, 0.6])
        p4 = Pendulum(length=5, mass=1, pos=[1.2, 1.6], vel=[0.4, 0.8])
        self.assertEqual(2*p2, p4)
        self.assertEqual(p4*0.5, p1*2)
        self.assertAlmostEqual(p3*(1/3),p1)

    def test_dT(self):
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0])

        q1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0], zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0],
                      zPolar=False)


        self.assertTrue(
            arrayAlmostEqual(p1.dT(), np.array([10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p2.dT(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(p3.dT(), np.array([0, 10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p4.dT(), np.array([0, 0, 10])) is None)

        self.assertTrue(
            arrayAlmostEqual(q1.dT(), np.array([0, 10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q2.dT(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q3.dT(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(q4.dT(), np.array([-10, 0, 0])) is None)


    def test_dP(self):
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0])

        q1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0], zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0],
                      zPolar=False)

        self.assertTrue(
            arrayAlmostEqual(p1.dP(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p2.dP(), np.array([0, 10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p3.dP(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p4.dP(), np.array([-10, 0, 0])) is None)

        self.assertTrue(
            arrayAlmostEqual(q1.dP(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q2.dP(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(q3.dP(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q4.dP(), np.array([0, -10, 0])) is None)

    def test_dT2(self):
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0])

        q1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0], zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0],
                      zPolar=False)

        self.assertTrue(
            arrayAlmostEqual(p1.dT2(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(p2.dT2(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p3.dT2(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(p4.dT2(), np.array([0, -10, 0])) is None)

        self.assertTrue(
            arrayAlmostEqual(q1.dT2(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q2.dT2(), np.array([0, -10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q3.dT2(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q4.dT2(), np.array([0, 0, -10])) is None)

    def test_dP2(self):
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0])

        q1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0], zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0],
                      zPolar=False)

        self.assertTrue(
            arrayAlmostEqual(p1.dP2(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p2.dP2(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p3.dP2(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p4.dP2(), np.array([0, -10, 0])) is None)

        self.assertTrue(
            arrayAlmostEqual(q1.dP2(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q2.dP2(), np.array([0, -10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q3.dP2(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q4.dP2(), np.array([0, 0, -10])) is None)

    def test_dTP(self):
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0])

        q1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0], zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0],
                      zPolar=False)

        self.assertTrue(
            arrayAlmostEqual(p1.dTP(), np.array([0, 10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p2.dTP(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p3.dTP(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p4.dTP(), np.array([0, 0, 0])) is None)

        self.assertTrue(
            arrayAlmostEqual(q1.dTP(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(q2.dTP(), np.array([0, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q3.dTP(), np.array([0, -10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q4.dTP(), np.array([0, 0, 0])) is None)

    def test_dzdT(self):
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0])

        q1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0], zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0],
                      zPolar=False)


        self.assertAlmostEqual(p1.dzdT(), 0)
        self.assertAlmostEqual(p2.dzdT(), 10)
        self.assertAlmostEqual(p3.dzdT(), 0)
        self.assertAlmostEqual(p4.dzdT(), 10)

        self.assertAlmostEqual(q1.dzdT(), 0)
        self.assertAlmostEqual(q2.dzdT(), 0)
        self.assertAlmostEqual(q3.dzdT(), 10)
        self.assertAlmostEqual(q4.dzdT(), 0)

    def test_dzdP(self):
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0])

        q1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0], zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q4 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0],
                      zPolar=False)


        self.assertAlmostEqual(p1.dzdP(), 0)
        self.assertAlmostEqual(p2.dzdP(), 0)
        self.assertAlmostEqual(p3.dzdP(), 0)
        self.assertAlmostEqual(p4.dzdP(), 0)

        self.assertAlmostEqual(q1.dzdP(), 0)
        self.assertAlmostEqual(q2.dzdP(), 10)
        self.assertAlmostEqual(q3.dzdP(), 0)
        self.assertAlmostEqual(q4.dzdP(), 0)

    def test_Velocity(self):
        #pendulum is parallel to the x-axis, with perpendicular angular
        #velocity in the positive z direction
        p1 = Pendulum(length=10, mass=1, pos=[np.pi/2,0], vel=[1,0])
        #same position but with a perpendicular angular velocity in the
        #direction of positive y
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2,0], vel=[0,1])
        #same position but both azimuthal and polar velocity. should result in
        #combination of p1 and p2
        p3 = Pendulum(length=10, mass=1, pos=[np.pi/2,0], vel=[1,1])
        #pendulum is parallel to the y-axis, velocity in the negative z
        #direction
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2,np.pi/2], vel=[-1,0])
        #same position as before but velocity is in the positive x direction
        p5 = Pendulum(length=10, mass=1, pos=[np.pi/2,np.pi/2], vel=[0,-1])

        #pendulum is parallel to the negative z-axis, has velocity in the
        #direction of positive x
        q1 = Pendulum(length=10, mass=1, pos=[-np.pi/2,np.pi/2], vel=[1, 0],
                      zPolar=False)
        #same position as q2 with a velocity in the direction of positive y
        q2 = Pendulum(length=10, mass=1, pos=[-np.pi/2,np.pi/2], vel=[0, 1],
                      zPolar=False)

        self.assertTrue(
            arrayAlmostEqual(p1.Velocity(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(p2.Velocity(), np.array([0, 10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p3.Velocity(), np.array([0, 10, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(p4.Velocity(), np.array([0, 0, -10])) is None)
        self.assertTrue(
            arrayAlmostEqual(p5.Velocity(), np.array([10, 0, 0])) is None)

        self.assertTrue(
            arrayAlmostEqual(q1.Velocity(), np.array([10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q2.Velocity(), np.array([0, 10, 0])) is None)

    def test_toCartesian(self):
        #positive z
        p1 = Pendulum(length=10, mass=1, pos=[np.pi, 0], vel=[0, 0])
        #positive x
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0,0])
        #positive y
        p3 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0,0])
        #negative z
        p4 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        #negative x
        p5 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi], vel=[0,0])
        #negative y
        p6 = Pendulum(length=10, mass=1, pos=[np.pi/2, -np.pi/2], vel=[0,0])

        #positive z
        q1 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi/2], vel=[0, 0],
                      zPolar=False)
        #positive x
        q2 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0,0],
                      zPolar=False)
        #positive y
        q3 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0,0],
                      zPolar=False)
        #negative z
        q4 = Pendulum(length=10, mass=1, pos=[np.pi/2, -np.pi/2], vel=[0, 0],
                      zPolar=False)
        #negative x
        q5 = Pendulum(length=10, mass=1, pos=[np.pi, 0], vel=[0,0],
                      zPolar=False)
        #negative y
        q6 = Pendulum(length=10, mass=1, pos=[np.pi/2, np.pi], vel=[0,0],
                      zPolar=False)

        self.assertTrue(
            arrayAlmostEqual(p1.toCartesian(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(p2.toCartesian(), np.array([10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p3.toCartesian(), np.array([0, 10, 0])) is None)

        self.assertTrue(
            arrayAlmostEqual(p4.toCartesian(), np.array([0, 0, -10])) is None)
        self.assertTrue(
            arrayAlmostEqual(p5.toCartesian(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(p6.toCartesian(), np.array([0, -10, 0])) is None)

        self.assertTrue(
            arrayAlmostEqual(q1.toCartesian(), np.array([0, 0, 10])) is None)
        self.assertTrue(
            arrayAlmostEqual(q2.toCartesian(), np.array([10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q3.toCartesian(), np.array([0, 10, 0])) is None)

        self.assertTrue(
            arrayAlmostEqual(q4.toCartesian(), np.array([0, 0, -10])) is None)
        self.assertTrue(
            arrayAlmostEqual(q5.toCartesian(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q6.toCartesian(), np.array([0, -10, 0])) is None)

    def test_fixCoord_rangeCorrection(self):
        #pendulum in the x-position, ideal encoding
        p1 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        #pendulum is in the x-position with phase shifts of 2pi and 4pi
        p2 = Pendulum(length=10, mass=1, pos=[2*np.pi + np.pi/2, 4*np.pi],
                      vel=[0,0])
        #also in the x-position but theta is outside of the desired range
        #of [0, pi) (duplicate coordinates)
        p3 = Pendulum(length=10, mass=1, pos=[3*np.pi/2, -np.pi], vel=[0,0])

        #proving that both p2 and p3 do not have the same attributes as p1, and
        #that all 3 are different pendula.
        self.assertFalse(p1 == p2)
        self.assertFalse(p3 == p1)
        self.assertFalse(p2 == p3)

        #proving that despite being different all three pendula give the same
        #cartesian position when .toCartesian() is called.
        self.assertTrue(
            arrayAlmostEqual(p1.toCartesian(),p2.toCartesian()) is None)
        self.assertTrue(
            arrayAlmostEqual(p3.toCartesian(),p1.toCartesian()) is None)

        #ensures that fixCoord does nothing if the pendulum is sufficiently 
        #far from a pole, and theta and phi are in the ranges [0, 2pi) and 
        #[0, pi) respectively.
        testP = deepcopy(p1)
        p1.fixCoord()
        self.assertEqual(p1, testP)

        #checks if .fixCoord() corrects for phase shifts which are multiples of
        #2 pi. Checking if p2 is corrected to be equal to p1
        p2.fixCoord()
        self.assertEqual(p2, p1)

        #checks if .fixCoord() corrects for theta being outside of the range
        #[0, pi) such that the pendula gets corrected to the same 'effective'
        #value. Basically checks if fixCoord() corrects for duplicate values.
        p3.fixCoord()
        self.assertEqual(p3, p1)

    def test_fixCoord_Transformation(self):
        #as the specific position and velocity shouldn't matter, this loops
        #over 500 random Pendulum objects and tranforms them, then checks that
        #the cartesian position and velocity aren't changed by the transform-
        #ation.
        for i in range(50):
            #creates a random value of theta between 0 and 1/6pi or 5/6pi and
            #pi
            theta = np.pi*(lambda a: a if a > 0 else a + 1)(
                r.uniform(-1/6,1/6))
            #a random azimuthal angle between 0 and 2pi
            phi = r.uniform(0, 2*np.pi)
            #range matters a bit less for angular speeds, so its just from 0
            #to 10.
            thetadot = r.uniform(0, 10)
            phidot = r.uniform(0, 10)
            #randomly pick which coordinate basis is used.
            pole = r.choice([True, False])

            #random attribute pendulum. Mass is chosen to track the iteration
            #in the case of an error.
            P = Pendulum(length=10, mass=i+1, pos=[theta, phi],
                         vel=[thetadot, phidot], zPolar=pole)
            Q = deepcopy(P)
            P.fixCoord()

            #ensuring the cartesian positon and velocity are unchanged.
            self.assertTrue(
                arrayAlmostEqual(P.toCartesian(), Q.toCartesian()) is None)
            self.assertTrue(
                arrayAlmostEqual(P.Velocity(), Q.Velocity()) is None)
            #ensuring that the pendulum has transformed.
            self.assertEqual(Q.zPolar, not P.zPolar)


class TestChainMethods(unittest.TestCase):

    def test_eq(self):
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        c1 = Chain(pList = [p1, p2], g=9.81)

        p3 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        c2 = Chain(pList = [p3, p4], g=9.81)

        self.assertTrue(c1 == c2)

    def test_mul(self):
        p1 = Pendulum(length=10, mass=1, pos=[0.3, 0.5], vel=[0.7, 1.1])
        p2 = Pendulum(length=10, mass=1, pos=[1.3, 1.7], vel=[1.9, 2.3])
        c1 = Chain(pList = [p1, p2], g=9.81)

        p3 = Pendulum(length=10, mass=1, pos=[0.6, 1.0], vel=[1.4, 2.2])
        p4 = Pendulum(length=10, mass=1, pos=[2.6, 3.4], vel=[3.8, 4.6])
        c2 = Chain(pList = [p3, p4], g=9.81)

        self.assertEqual(c2, c1*2)
        self.assertEqual(0.5*c2, c1)

    def test_add(self):
        p1 = Pendulum(length=10, mass=1, pos=[0.03, 0.05], vel=[0.07, 0.11])
        p2 = Pendulum(length=10, mass=1, pos=[0.13, 0.17], vel=[0.19, 0.23])
        c1 = Chain(pList = [p1, p2], g=9.81)

        p3 = Pendulum(length=10, mass=1, pos=[0.06, 0.10], vel=[0.14, 0.22])
        p4 = Pendulum(length=10, mass=1, pos=[0.26, 0.34], vel=[0.38, 0.46])
        c2 = Chain(pList = [p3, p4], g=9.81)

        p5 = Pendulum(length=10, mass=1, pos=[0.09, 0.15], vel=[0.21, 0.33])
        p6 = Pendulum(length=10, mass=1, pos=[0.39, 0.51], vel=[0.57, 0.69])
        c3 = Chain(pList = [p5, p6], g=9.81)

        self.assertAlmostEqual((c1+c2), c3)

    def test_AElement(self):
        #three pendula with different positions, to be used as reference.
        p1 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[np.pi/4, np.pi], vel=[0, 0])

        #every combination of the pendula is tested against a pre-calculated
        #result.
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p1),
                             100*np.array([[1,0],[0,1]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p2, p2),
                             100*np.array([[1,0],[0,0]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p3, p3),
                             100*np.array([[1,0],[0,1/2]])) is None)

        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p2),
                             100*np.array([[0,0],[1,0]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p2, p1),
                             100*np.array([[0,1],[0,0]])) is None)

        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p3),
                             50*np.sqrt(2)*np.array([[1,0],[0,-1]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p3, p1),
                             50*np.sqrt(2)*np.array([[1,0],[0,-1]])) is None)

        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p2, p3),
                             -50*np.sqrt(2)*np.array([[0,1],[0,0]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p3, p2),
                             -50*np.sqrt(2)*np.array([[0,0],[1,0]])) is None)

        #three more pendula with the same values of theta and phi as their
        #respective 'p' pendula but in the other coordinate basis. If the
        #mathematics is consistent, the result of combining these pendula
        #should be the same as the other pendula, as the two coordinate basis
        #can be transformed between each other with an orthogonal matrix,
        #so inner products should be conserved.
        q1 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[np.pi/4, np.pi], vel=[0, 0],
                      zPolar=False)

        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p1),
                             Chain._AElement(q1, q1)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p2, p2),
                             Chain._AElement(q2, q2)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p3, p3),
                             Chain._AElement(q3, q3)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p2),
                             Chain._AElement(q1, q2)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p2, p1),
                             Chain._AElement(q2, q1)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p3),
                             Chain._AElement(q1, q3)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p3, p1),
                             Chain._AElement(q3, q1)) is None)

    def test_BElement(self):
        #three pendula with different theta and phi values, to be used as
        #reference.
        p1 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        p3 = Pendulum(length=10, mass=1, pos=[np.pi/4, np.pi], vel=[0, 0])

        #every combination of p1, p2, and p3 is tested against a pre-calculated
        #result
        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p1, p1),
                             100*np.array([[0,0,0],[0,0,0]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p2, p2),
                             100*np.array([[0,0,0],[0,0,0]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p3, p3),
                             50*np.array([[0,0,-1],[0,1,0]])) is None)

        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p1, p2),
                             100*np.array([[1,0,0],[0,0,0]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p2, p1),
                             100*np.array([[0,0,0],[0,0,0]])) is None)

        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p1, p3), 50*np.sqrt(2)*
                             np.array([[1,0,0],[0,-1,0]]))is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p3, p1),50*np.sqrt(2)*\
                             np.array([[1,0,1],[0,0,0]])) is None)

        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p2, p3),-50*np.sqrt(2)*
                             np.array([[0,1,0],[0,0,0]])) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._BElement(p3, p2),50*np.sqrt(2)*
                             np.array([[1,1,0],[0,0,0]])) is None)

        #three more pendula with the same values of theta and phi as their
        #respective 'p' pendula but in the other coordinate basis. If the
        #mathematics is consistent, the result of combining these pendula
        #should be the same as the other pendula, as the two coordinate basis
        #can be transformed between each other with an orthogonal matrix,
        #so inner products should be conserved.
        q1 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0],
                      zPolar=False)
        q2 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0],
                      zPolar=False)
        q3 = Pendulum(length=10, mass=1, pos=[np.pi/4, np.pi], vel=[0, 0],
                      zPolar=False)

        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p1),
                             Chain._AElement(q1, q1)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p2, p2),
                             Chain._AElement(q2, q2)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p3, p3),
                             Chain._AElement(q3, q3)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p2),
                             Chain._AElement(q1, q2)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p2, p1),
                             Chain._AElement(q2, q1)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p1, p3),
                             Chain._AElement(q1, q3)) is None)
        self.assertTrue(
            arrayAlmostEqual(Chain._AElement(p3, p1),
                             Chain._AElement(q3, q1)) is None)


    def test_Matrices(self):
        p1 = Pendulum(length=10, mass=100, pos=[3.45610, 0.0365], vel=[0, 0])
        p2 = Pendulum(length=10, mass=50, pos=[1, 4], vel=[0, 0])
        p3 = Pendulum(length=10, mass=30, pos=[np.pi/4, np.pi], vel=[0, 0])
        c1 = Chain(pList=[p1, p2, p3], g=9.81)

        #creates a two matrices out of _AElement and _BElement
        A, B = c1._Matrices()

        #checks every block of A and B to check that they're a matrix block
        #element scaled by the sum of the masses of pendula with equal or
        #greater index

        #this element only includes p1 and so the block elements are scaled by
        #the sum of p1.mass, p2.mass, and p3.mass
        self.assertTrue(
            arrayEqual(A[0:2,0:2], 180*Chain._AElement(p1,p1)) is None)
        self.assertTrue(
            arrayEqual(B[0:2,0:3], 180*Chain._BElement(p1,p1)) is None)

        #this element includes p2, so the mass of p1 cannot be included, so the
        #scaling factor is only 80
        self.assertTrue(
            arrayEqual(A[2:4,0:2], 80*Chain._AElement(p2,p1)) is None)
        self.assertTrue(
            arrayEqual(B[2:4,0:3], 80*Chain._BElement(p2,p1)) is None)

        #as this element includes p3, it can only be scaled by the mass of p3
        self.assertTrue(
            arrayEqual(A[4:6,0:2], 30*Chain._AElement(p3,p1)) is None)
        self.assertTrue(
            arrayEqual(B[4:6,0:3], 30*Chain._BElement(p3,p1)) is None)

        #second 'row' of block elements
        self.assertTrue(
            arrayEqual(A[0:2,2:4], 80*Chain._AElement(p1,p2)) is None)
        self.assertTrue(
            arrayEqual(B[0:2,3:6], 80*Chain._BElement(p1,p2)) is None)

        self.assertTrue(
            arrayEqual(A[2:4,2:4], 80*Chain._AElement(p2,p2)) is None)
        self.assertTrue(
            arrayEqual(B[2:4,3:6], 80*Chain._BElement(p2,p2)) is None)

        self.assertTrue(
            arrayEqual(A[4:6,2:4], 30*Chain._AElement(p3,p2)) is None)
        self.assertTrue(
            arrayEqual(B[4:6,3:6], 30*Chain._BElement(p3,p2)) is None)

        #third 'row' of block elements
        self.assertTrue(
            arrayEqual(A[0:2,4:6], 30*Chain._AElement(p1,p3)) is None)
        self.assertTrue(
            arrayEqual(B[0:2,6:9], 30*Chain._BElement(p1,p3)) is None)

        self.assertTrue(
            arrayEqual(A[2:4,4:6], 30*Chain._AElement(p2,p3)) is None)
        self.assertTrue(
            arrayEqual(B[2:4,6:9], 30*Chain._BElement(p2,p3)) is None)

        self.assertTrue(
            arrayEqual(A[4:6,4:6], 30*Chain._AElement(p3,p3)) is None)
        self.assertTrue(
            arrayEqual(B[4:6,6:9], 30*Chain._BElement(p3,p3)) is None)

    def test_GVec(self):
        p1 = Pendulum(length=10, mass=5, pos=[np.pi/2, 0], vel=[0, 0])
        p2 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0, 0])
        c1 = Chain(pList=[p1, p2], g=9.81)

        V1 = c1._GVec()
        V2 = np.array([9.81*6*10, 0, 0, 0])

        self.assertTrue(
            arrayAlmostEqual(V1, V2) is None)

        p3 = Pendulum(length=10, mass=5, pos=[np.pi/2, 0], vel=[0, 0])
        p4 = Pendulum(length=10, mass=1, pos=[np.pi/4, np.pi/2], vel=[0, 0])
        c2 = Chain(pList=[p3, p4], g=9.81)

        V3 = c2._GVec()
        V4 = np.array([9.81*6*10, 0, 9.81*1*10/np.sqrt(2), 0])

        self.assertTrue(
            arrayAlmostEqual(V3, V4) is None)

        p5 = Pendulum(length=10, mass=5, pos=[np.pi/4, np.pi/4], vel=[0, 0],
                      zPolar=False)
        p6 = Pendulum(length=10, mass=1, pos=[np.pi/2, 0], vel=[0, 0])
        c3 = Chain(pList=[p5, p6], g=9.81)

        V5 = c3._GVec()
        V6 = np.array([9.81*6*10/2, 9.81*6*10/2, 9.81*1*10, 0])

        self.assertTrue(
            arrayAlmostEqual(V5, V6) is None)

    def test_Func(self):
        #this function calculates the lagrangian differential equations from
        #matrices A, B, and GVec and then pipes the acceleration result back
        #into the new velocities, and simultaneously takes the old velocities
        #and pipes them into the new positions. This function is intended to
        #work in harmondy with runge-kutte methods.

        #I have independently worked out the differential equations for a 
        #single spherical pendulum, and so I have calculated the formulae for
        #the expected polar and azimuthal acceleration for given positions and
        #velocities.

        #formula is given by ddt = sin(t)cos(t)dt^2 + g/l * sin(t) where t is
        #the polar angle, and dt is the polar acceleration
        PolarAcc = (lambda t,dt,dp,g,l: np.sin(t)*(np.cos(t)*dp**2-g/l))
        #formula is given by ddp = -cot(t)*2*dt*dp where t is the polar angle,
        #and dt and dp are the polar and azimuthal velocities.
        AzimuthAcc = (lambda t,dt,dp,g,l: -(np.cos(t)/np.sin(t))*2*dt*dp)

        p1 = Pendulum(pos=[3.61923, 0.11111], vel=[2.2222, 5.3451], length=10,
                      mass=1)
        c1 = Chain(pList=[p1], g=9.81)

        c2 = Chain._Func(c1)
        p2 = c2.pList[0]

        self.assertEqual(p1.vel[0], p2.pos[0])
        self.assertEqual(p1.vel[1], p2.pos[1])

        p1PolarAcc = PolarAcc(p1.pos[0], p1.vel[0], p1.vel[1], 9.81, 10)
        p1AzimuthAcc = AzimuthAcc(p1.pos[0], p1.vel[0], p1.vel[1], 9.81, 10)

        self.assertAlmostEqual(p1PolarAcc, p2.vel[0])
        self.assertAlmostEqual(p1AzimuthAcc, p2.vel[1])
    
    def test_RK(self):
        #arbitrary chain
        p1 = Pendulum(length=10, mass=1, pos=[0, 0], vel=[0,0])
        p2 = Pendulum(length=10, mass=1, pos=[0, np.pi/2], vel=[0,0])
        c1 = Chain(pList=[p1, p2], g=9.81)
        c1.fixCoords()

        #explicit euler method
        c2 = c1 + 0.5 * Chain._Func(c1)

        #using the RK method to do the euler method
        c1.RK(timeDelta=0.5, RKmatrix=np.array([]), weights=[1],\
            function=Chain._Func)

        self.assertEqual(c1,c2)

        #arbitrary chain
        p3 = Pendulum(length=10, mass=1, pos=[0.11111, 0.15325], vel=[5,2])
        p4 = Pendulum(length=10, mass=1, pos=[3, 4], vel=[0.321,0.123])
        c3 = Chain(pList=[p3, p4], g=9.81)
        c3.fixCoords()

        #explicit euler method
        c4 = c3 + 0.3333 * Chain._Func(c3)

        #using the RK method to do the euler method
        c3.RK(timeDelta=0.3333, RKmatrix=np.array([]), weights=[1],\
            function=Chain._Func)

        self.assertEqual(c3,c4)

        #arbitrary chain
        p5 = Pendulum(length=10, mass=1, pos=[0.3333, 7.15325], vel=[5,2])
        p6 = Pendulum(length=10, mass=1, pos=[5.2, 4.0101], vel=[10,0])
        c5 = Chain(pList=[p5, p6], g=9.81)
        c5.fixCoords()

        #explicit midpoint method
        c6 = c5 + 0.1234 * Chain._Func(c5 + 0.5*0.1234*Chain._Func(c5))

        #using the RK method to do the midpoint method
        c5.RK(timeDelta=0.1234, RKmatrix=np.array([[0.5]]), weights=[0,1],\
            function=Chain._Func)

        self.assertEqual(c5,c6)


if __name__ == "__main__":
    unittest.main()
