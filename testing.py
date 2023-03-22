from Pendulum import Pendulum
from Chain import Chain
import numpy as np
from numpy.testing import assert_array_almost_equal as arrayAlmostEqual
from numpy.testing import assert_array_equal as arrayEqual
import unittest
import random as rand

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
        #self.assertTrue(
        #    arrayAlmostEqual(p4.dT(), np.array([0, 0, 10])) is None)

        self.assertTrue(
            arrayAlmostEqual(q1.dT(), np.array([0, 10, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q2.dT(), np.array([-10, 0, 0])) is None)
        self.assertTrue(
            arrayAlmostEqual(q3.dT(), np.array([0, 0, 10])) is None)
        #self.assertTrue(
        #    arrayAlmostEqual(q4.dT(), np.array([-10, 0, 0])) is None)


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


if __name__ == "__main__":
    unittest.main()
