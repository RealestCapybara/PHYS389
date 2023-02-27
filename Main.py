import numpy as np
from numpy.linalg import norm

class Quat():
    """Quaternion object. Requires 4 components, 1 real, 3 imaginary.
    Vector Form: q = 0 + x i + y j + z k
    Rotation Form: q = cos(A/2) + sin(A/2)[n_x i + n_y j + n_z k]
    """
    def __init__(self, q0, q1, q2, q3):
        try:
            self.q0 = np.float64(q0)
            self.q1 = np.float64(q1)
            self.q2 = np.float64(q2)
            self.q3 = np.float64(q3)
        except:
            raise TypeError("Quaternion can only take numeric input data")

    def __mul__(self, other):
        """Defines Quaternion self multiplication and scaling."""
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
        """Defined in case a scaler is on the other side of the product"""
        if type(other) == type(self):
            return other.__mul__(self)
        else:
            return self.__mul__(other)

    def __add__(self, other):
        """Defines the termwise addition of Quaternions"""
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
        return str(Quat(format(self.q0, spec), format(self.q1, spec),
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

    def __matmul__(self, other):
        """Encodes a 'dot product' functionality to quaternions."""
        return self.q0 * other.q0 + self.q1 * other.q1 + self.q2 * other.q2 +
                    self.q3 * other.q3

    def conjugate(self):
        """Returns the conjugate of the Quaterion object"""
        return Quat(self.q0, -self.q1, -self.q2, -self.q3)

class Pend():
    """Pendulum Object.
    mass (numeric) default=1.0 Encodes mass (both inertial and gravitational)
    pos (Quat) default=Quat(0, 0, 0, 1) Encodes length and orientation
    vel (Quat) default=Quat(0, 0, 0, 0) Encodes velocity
    acc (Quat) default=Quat(0, 0, 0, 0) Encodes acceleration
    """
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
        self.parent = None
        self.child = None
        self.sibling = None

    def rotate(self, delta):
        """Rotates the position by a given Delta (change in position)
        Needs very small delta due to use of trig approximations"""
        pos = self.pos
        #for two "vector form" quaternions, their product is analogous to a 
        #cross product. Therefore by taking the quaternion product of the 
        #position and the delta, a new quaternion parallel to the axis of 
        #rotation and equal to the product of the radius and the angular 
        #component of delta
        q = self.pos * delta

        posA = abs(pos)
        qA = abs(q)

        #to turn q into a rotation quaternion, the angle can be calculated
        #by using the tan(theta)=theta approximation, which results in 
        #theta = |q|/(|pos|^2) which can then be used to scale q by 
        #sin(theta/2)/|q| (1/|q| to account for the necessary normalisation)
        #and added to cos(theta/2). Quadratic trig approximations are used for
        #both. As the delta is expected to be small, their use is tolerable.
        #This means that errors will be expected to scale cubicly.
        q = q * (1/(2*posA**2))
        q = q + Quat(1-(qA**2/(8*posA**4)), 0, 0, 0)

        #Quaternion is then normalised to ensure that it is a unitary
        #transformation and doesn't modify the length of the pendulum
        q = q * (1/abs(q))

        #once the rotation quaternion is calculated, all that is needed to
        #apply it is simple quaternion multiplication
        self.pos = q * pos * q.conjugate()

    def TensionOnParent(self, inplace=True):
        x = self.pos
        a = self.acc
        T = (x @ a)/(abs(x)**2) * x
        if inplace == True:
            self.parent.acc += T
        else:
            return T

class GenList():
    def __init__(self, root):
        if type(root) != Pend:
            raise TypeError("root must be Pendulum")

        currentGen = 0
        List = [[root]]

        while True:
            NextGen = []
            for pend in List[currentGen]:
                child = pend.child
                while child is not None:
                    NextGen.append(child)
                    child = child.sibling
            if len(NextGen) == 0:
                break
            else:
                List.append(NextGen)
                currentGen += 1

        self._List = np.array(List)
        self._rList = np.array(List.reverse())

    def __iter__(self):
        return self._List
    def __getitem__(self, idx):
        return self._List[idx]
    def __str__(self):
        return f"GenList({str(self._List)})"
    def __repr__(self):
        return f"GenList({repr(self._List)})"

    def UpdateAcc(self):
        for Gen in self._rList:
            for P in Gen:
                P.acc = Quat(0, 0, 0, self.mass * 9.81)

        for Gen in self._rList:
            for P in Gen:
                P.TensionOnParent()

if __name__ == "__main__":
    pass
