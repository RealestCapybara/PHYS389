from Pendulum import Pendulum
from Chain import Chain
import numpy as np
import unittest

class TestPendulumMethods(unittest.TestCase):

    def test_eq(self):
        p1 = Pendulum(length=10, mass=1, pos=[0.5, 0.5], vel=[0.1, 0.1])
        p2 = Pendulum(length=10, mass=1, pos=[0.5, 0.5], vel=[0.1, 0.1])
        self.assertEqual(p1, p2)

    def test_add(self):
        p1 = Pendulum(length=10, mass=1, pos=[0.5, 0.5], vel=[0.1, 0.1])
        p2 = Pendulum(length=10, mass=1, pos=[0.6, 0.7], vel=[0.2, 0.3])
        p3 = Pendulum(length=10.0, mass=1.0, pos=[1.1, 1.2], vel=[0.3, 0.4])
        print(p1+p2)
        self.assertEqual(p1 + p2, p3)
        self.assertEqual(p1 + p2, p2 + p1)
        self.assertEqual((p1+p2).length, 10)
        self.assertEqual((p1+p2).mass, 1)
        self.assertEqual(p1+0, p1)
        self.assertEqual(None+p2, p2)

    def test_mul(self):
        p1 = Pendulum(length=5, mass=1, pos=[0.6, 0.8], vel=[0.2, 0.4])
        p2 = Pendulum(length=5, mass=1, pos=[0.3, 0.4], vel=[0.1, 0.2])
        p3 = Pendulum(length=5, mass=1, pos=[1.2, 1.6], vel=[0.4, 0.8])
        self.assertEqual(2*p1, p3)
        self.assertEqual(p3*0.5, p2*2)

if __name__ == "__main__":
    unittest.main()
