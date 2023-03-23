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



if __name__ == "__main__":
    unittest.main()
