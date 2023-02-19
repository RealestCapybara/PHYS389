import numpy as np
from numpy.linalg import norm

class Quat():
    def __init__(self, q0, q1, q2, q3):
        self.q0 = q0
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3

    def __mul__(self, other):
        if type(other) == type(self):
            q0 = self.q0
            q1 = self.q1
            q2 = self.q2
            q3 = self.q3

            p0 = other.q0
            p1 = other.q1
            p2 = other.q2
            p3 = other.q3
            return Quat(q0*p0-q1*p1-q2*p2-q3*p3,
                        q0*p1+q1*p0+q2*p3-q3*p2,
                        q0*p2-q1*p3+q2*p0+q3*p1,
                        q0*p3+q1*p2-q2*p1+q3*p0)
        else:
            try:
                return Quat(self.q0*other, self.q1*other,
                             self.q2*other, self.q3*other)
            except TypeError:
                raise TypeError(f"\
Quaternion cannot be multiplied with object of type {type(other)}\n\
Can only be multiplied by another quaternion or scaled")

    def __rmul__(self, other):
        if type(other) == type(self):
            return other.__mul__(self)
        else:
            return self.__mul__(other)

    def __add__(self, other):
        if type(other) == type(self):
            return Quat(self.q0+other.q0, self.q1+other.q1,
                         self.q2+other.q2, self.q3+other.q3)
        else:
            raise TypeError("\
Quaternion cannot be summed with non-Quaternion object")
    def __str__(self):
        return f"{self.q0} + {self.q1} i + {self.q2} j + {self.q3} k"

    def __repr__(self):
        return f"Quat({self.q0}, {self.q1}, {self.q2}, {self.q3})"

    def __format__(self, spec):
        return str(Quart(format(self.q0, spec), format(self.q1, spec),
                     format(self.q2, spec), format(self.q3, spec)))

    def __getitem__(self, key):
        if key == 0:
            return self.q0
        elif key == 1:
            return self.q1
        elif key == 2:
            return self.q2
        elif key == 3:
            return self.q3
        else:
            raise KeyError("Invalid key for indexing Quaternion object")

    def __abs__(self):
        return np.sqrt(self.q0**2 + self.q1**2 + self.q2**2 + self.q3**2)

    def __iter__(self):
        for val in [self.q0, self.q1, self.q2, self.q3]:
            yield val

    def conjugate(self):
        return Quat(self.q0, -self.q1, self.q2, -self.q3)

class Pend():
    def __init__(self, mass=1.0, pos = Quat(0, 0, 0, 1),
                 vel = Quat(0, 0, 0, 0), acc = Quat(0, 0, 0, 0)):
        """Pendulum Object"""
        for name, val in [('pos', type(pos)), ('vel', type(vel)),
                          ('acc', type(acc))]:
            if val != Quat:
                raise TypeError(f"\
{name} cannot be of type {val}, must be of type Quat")

        self.pos = pos
        self.vel = vel
        self.acc = acc
        self.mass = mass

    def rotate(self, delta):
        pos = self.pos
        q = self.pos * delta

        posA = abs(pos)
        qA = abs(q)

        q = q * (1/(2*posA**2))
        q = q + Quat(1-(qA**2/(8*posA**4)), 0, 0, 0)
        self.pos = q * pos * q.conjugate()

if __name__ == "__main__":
    X = Pend(pos=Quat(0, 0, 0, 10))
    Delta = Quat(0, 0.01, 0.01, 0.0)
    X.rotate(Delta)
    print(X.pos)
    print(abs(X.pos))
