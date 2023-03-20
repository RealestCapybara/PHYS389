import numpy as np
from Pendulum import Pendulum

class Chain():
    """
    Represents a chain of connected pendulums.

    Parameters
    ----------
    pList : list of Pendulum objects
        List of Pendulum objects that form the chain
    g : float or int
        The acceleration due to gravity

    Raises
    ------
    TypeError:
        If pList contains objects that are not instances of Pendulum or if g is
        not a scalar value.

    Attributes
    ----------
    pList : list of Pendulum objects
        List of Pendulum objects that form the chain
    g : numpy float64
        The acceleration due to gravity

    Methods
    -------
    RK(timeDelta, weights, RKmatrix, function)
        A generic Runge-Kutta method function to calculate the next position
    RK4(timeDelta)
        A 4th order Runge-Kutte method function to calculate the next position
    fixCoords()
        Iterates through the pendula and corrects their coordinate system
    _Matrices()
        Returns the A and B matrices that are used in the equations of motion
    _GVec()
        Returns the gravitational force vector for each pendulum in the chain

    Static Methods
    --------------
    _AElement(p1, p2)
        Returns an element of the block matrix A
    _BElement(p1, p2)
        Returns an element of the block matrix B
    _Func(Matr)
        Calculates equations of motion using A and B matrix and GVec

    Dunder Methods
    --------------
    __mul__(self, other)
        Scalar multiplication of the chain
    __rmul__(self, other)
        Accounts for a multiplication of a scalar from the other side
    __add__(self, other)
        Addition of the chain with another Chain object
    __str__(self)
        Returns a string representation of the chain

    """
    def __init__(self, pList, g):
        """
        Initializes a Chain object with the given list of pendula and a
        gravitational constant.

        Parameters
        ----------
        pList : list of Pendulum objects
            A list of Pendulum objects representing the pendulums in the chain.
        g : float or int
            The gravitational constant to use for the simulation.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If pList is not a list of Pendulum objects, or if g is not a
            scalar value.
        """
        if not all(isinstance(p, Pendulum) for p in pList):
            raise TypeError("pList must be a list of Pendulum objects")
        if not isinstance(g, (int, float, np.float64)):
            raise TypeError("g must be a scalar value")
        self.pList = pList
        self.g = np.float64(g)

    def __mul__(self, other):
        """Multiply a Chain object with a scalar value.

        This method creates a new Chain object with the same pendulums as the
        original object, but with each pendulum's angular position and velocity
        scaled by the given value.

        Parameters
        ----------
        other : float or int
            The scalar value to multiply each Pendulum by.

        Returns
        -------
        Chain
            A new Chain object with the same pendulums as the original object,
            but with each Pendulum object's angular position and velocity
            multiplied by the given value.

        Raises
        ------
        TypeError
            If the other parameter is not a float-like object
        """
        try:
            other = np.float64(other)
        except TypeError:
            raise TypeError("\
Chain can only be multiplied with float-like types")
        return Chain(pList = [p*other for p in self.pList], g=self.g)

    def __rmul__(self, other):
        """Multiplication by a scalar value if the Chain object is on the right

        Returns the same result as would be obtained with the Chain object on
        the left side of a multiplication (self.__mul__(other))

        Parameters
        ----------
        other : float or int
            The scalar value to multiply each Pendulum by.

        Returns
        -------
        Chain
            A new Chain object with the same pendulums as the original object,
            but with each Pendulum object's angular position and velocity
            multiplied by the given value.

        Raises
        ------
        TypeError
            If the other parameter is not a float-like object
        """
        return self.__mul__(other)

    def __add__(self, other):
        """
        Adds this Chain to another Chain, element-wise.

        adds self.pList to other.pList element wise and returns a Chain object
        containing the result, with the same value of g as the Chain object on
        the left (self). Upon each termwise Pendulum addition returns a new
        Pendulum with a combined (sum) angular position and velocity. If the
        other parameter is False-type (0, "", None, False, etc) then the Chain
        itself is returned without modification.

        Parameters
        ----------
        other: Chain or False-type object
            Chain to add to this Chain element-wise. If a False-type object is
            passed, this method returns this Chain itself.

        Returns
        -------
        Chain:
            A new Chain object, with the pList attribute contains a termwise
            sum of the two composing Chain objects' pList attributes. If other
            is a False-type object then the Chain itself is returned.

        Raises
        ------
        TypeError:
            If other is not a Chain object or a False-type object.
        """
        if type(other) is Chain:
            return Chain(pList = [p + q for p, q in
                                zip(self.pList,other.pList)], g = self.g)
        elif not other:
            return self
        else:
            raise TypeError("Chain can only be added to other Mat objects")

    def __radd__(self, other):
        return self.__add__(other)

    def __str__(self):
        """
        Returns a string representation of the Chain object.

        The string contains information about the position and velocity of each
        pendulum in the chain, as well as which . The position and velocity
        values are converted from radians to a product of pi and rounded to 4
        decimal places.

        Returns
        -------
        str
            A string representation of the Chain object, with each pendulum's
            information separated by commas and newlines for readability.
        """
        return ",\n         ".join([f"\
p{i}: pos=({v.pos[0]/np.pi:.4f} pi, {v.pos[1]/np.pi:.4f} pi) \
vel=({v.vel[0]/np.pi:.4f} pi, {v.vel[1]/np.pi:.4f} pi) \
pole={(lambda zPole: 'z' if zPole else 'x')(v.zPolar)}"
                                    for i, v in enumerate(self.pList)])

    @staticmethod
    def _AElement(p1, p2):
        """Compute an A matrix element from two Pendulum objects.

        Given two Pendulum objects p1 and p2, this method computes a block
        matrix element. This block element is used to create the A matrix
        used in the computation of the dynamics of the chain of Pendulum
        objects.

        Parameters
        ----------
        p1, p2 : Pendulum objects
            The two Pendulum objects used for the block element.

        Returns:
        --------
        A : np.ndarray, shape (2, 2)
            A block element of the 'A' matrix.

        Raises:
        -------
        TypeError
            If either p1 or p2 is not a Pendulum object.
        """
        if type(p1) is not Pendulum or type(p2) is not Pendulum:
            raise TypeError("arguments must be of type Pendulum")
        return np.array([[p2.dT()@p1.dT(), p2.dP()@p1.dT()],
                         [p2.dT()@p1.dP(), p2.dP()@p1.dP()]])

    @staticmethod
    def _BElement(p1, p2):
        """Compute a B matrix element from two Pendulum objects.

        Given two Pendulum objects p1 and p2, this method computes a block
        matrix element. This block element is used to create the B matrix
        used in the computation of the dynamics of the chain of Pendulum
        objects.

        Parameters
        ----------
        p1, p2 : Pendulum objects
            The two Pendulum objects used for the block element.

        Returns:
        --------
        A : np.ndarray, shape (2, 3)
            A block element of the 'B' matrix.

        Raises:
        -------
        TypeError
            If either p1 or p2 is not a Pendulum object.
        """
        if type(p1) is not Pendulum or type(p2) is not Pendulum:
            raise TypeError("arguments must be of type Pendulum")
        return np.array([[p2.dT2()@p1.dT(),p2.dTP()@p1.dT(),p2.dP2()@p1.dT()],
                         [p2.dT2()@p1.dP(),p2.dTP()@p1.dP(),p2.dP2()@p1.dP()]])

    def _Matrices(self):
        """
        Calculate the matrices A and B used to calculate the equations of
        motion for the chain of pendula.

        Uses both _AElement(p1, p2) and _BElement(p1, p2) over a range of
        Pendulum objects in the pList, and uses np.block() to combine the block
        elements into a larger matrix, with a block element for each pair of
        Pendulum objects.

        Returns
        -------
        Tuple
            A 2D numpy array representing the matrix A, and a 2D numpy array
            representing the matrix B.
        """
        pList = self.pList
        MassList = [p.mass for p in pList]
        arrayA = []
        arrayB = []
        for i, vi in enumerate(pList):
            rowA = []
            rowB = []
            for j, vj in enumerate(pList):
                M = sum(MassList[max(i, j):])
                ElementA = self._AElement(vi, vj) * M
                ElementB = self._BElement(vi, vj) * M
                rowA.append(ElementA)
                rowB.append(ElementB)
            arrayA.append(rowA)
            arrayB.append(rowB)
        return np.block(arrayA), np.block(arrayB)

    def _GVec(self):
        """
        Computes a vector that accounts for the force due to gravity on pendula
        in the chain in the differential equations.

        The vector has 2n elements, where n is the number of pendulums in the
        chain. The elements with an even index correspond to the partial
        derivative of the z-component of the pendulum's cartesian position with
        respect to the polar angle, and the odd elements correspond to partial
        derivatives of the pendulum's cartesian position with respect to the
        azimuthal angle.

        Returns:
        ---------
        numpy.ndarray
            A 1D array of shape (2n,), where n is the number of pendula in the
            chain.
        """
        g = self.g
        pList = self.pList
        array = []
        massList = [p.mass for p in pList]
        for i, q in enumerate(pList):
            M = sum(massList[i:])
            array.append(g*M*q.dzdT())
            array.append(g*M*q.dzdP())
        return np.array(array)

    @staticmethod
    def _Func(chain):
        """
        This static method takes in a Chain object and returns a new Chain
        object after updating the positions and velocities of the pendulums in
        the chain.

        This function is derived from calculating the lagrangian of the
        pendulum chain where the pendula are thought of as cartesian vectors
        with the x, y, and z components represented as functions of theta and
        phi. It can be determined that the differential equations created via
        this method can be represented in matrix form, with two matrices and
        three vectors. The matrices are composed of block elements created by
        taking the inner products of the partial derivatives of the cartesian
        pendulum vectors with respect to theta and phi.

        Of the three vectors,
        one is a vector that accounts for the gravitational acceleration (V3),
        one is a vector that accounts for the current angular velocity values
        (V2), and one is the output that accounts for the angular acceleration
        (V1). The vector V3 is calculated with the _GVec() method, the vector
        V2 has length 3n where n is the number of pendula, the 0 mod 3 index
        values are the square of the theta velocity, the 1 mod 3 index values
        are double the product of the phi velocity and the theta velocity, and
        the 2 mod 3 index values are the square of the phi velocity.

        These matrices and vectors forms a matrix equation that can be
        expressed by the following:

        V1 = A^-1 (B V2 + V3)

        The current values of the pendulum velocities are then taken as updated
        position values and the values from V1 are taken as updated velocity
        values. These are put into a Chain object which is then returned.

        Parameters:
        -----------
        chain : Chain object
            The Chain object whose positions and velocities need to be updated.

        Returns:
        --------
        M : Chain object
            The new Chain object with updated positions and velocities of the
            pendula.

        Raises:
        -------
        TypeError:
            If the argument is not of type Chain.

        """
        if type(chain) is not Chain:
            raise TypeError("argument must be of type Mat")
        pList = chain.pList
        g = chain.g

        A, B = chain._Matrices()
        V3 = chain._GVec()

        V2 = [f(p) for p in pList for f in (lambda p: p.vel[0]**2,
                                            lambda p: 2*p.vel[0]*p.vel[1],
                                            lambda p: p.vel[1]**2)]

        V2 = np.array(V2)

        V1 = -np.linalg.inv(A) @ (B@V2 + V3)

        omegaThetaDot = [v for i, v in enumerate(V1) if i%2 == 0]
        omegaPhiDot = [v for i, v in enumerate(V1) if i%2 == 1]

        ThetaDot = [p.vel[0] for p in pList]
        PhiDot = [p.vel[1] for p in pList]

        newPList = [Pendulum(pos=[ThetaDot[i], PhiDot[i]],
                             vel=[omegaThetaDot[i], omegaPhiDot[i]],
                             length=p.length, mass=p.mass)
                    for i, p in enumerate(pList)]

        M = Chain(pList = newPList, g=g)
        return M

    def RK(self, timeDelta, weights, RKmatrix, function):
        """
        Uses a generic Runge-Kutta method to numerically solve a system of
        differential equations for the Chain object. Using a suitable butcher
        tableau and function, this method can be used to update the system by
        one timestep.

        Parameters
        ----------
        timeDelta : float or int
            The time step size for the integration.
        weights : array-like
            An array of weights used to scale the test points before summing.
        RKmatrix : 2D array or empty array
            Representation of the main matrix in a butcher tableau formulation
            of a specific Runge-Kutte method.
        function : function
            A function that takes the Chain object and updates it. This needs
            to encode the system of differential equations desired.

        Raises
        ------
        TypeError
            If weights and RKmatrix are not array-like, or if timeDelta is not
            float-like, or if function is not a function.
        ValueError
            If the length of weights is not one more a row of
            RKmatrix, or if RKmatrix is not a square matrix, or if the sum of
            the weights is not 1.

        Returns
        -------
        None
        """
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

        if len(RKmatrix[-1]) != (len(weights)-1):
            raise ValueError(
                "length of weights must be one more than last row of RKmatrix")

        if len(RKmatrix) != 0:
            if RKmatrix.ndim != 2:
                raise ValueError("RKmatrix must be a matrix")
            if len(RKmatrix) != len(RKmatrix[0]):
                raise ValueError(
                    "number of rows and cols of RKmatrix must be equal")

        if round(sum(weights), 10) != 1.0:
            raise ValueError("sum of weights must be 1")

        k = [function(self)]
        x = self
        for i, v in enumerate(RKmatrix):
            delta = sum([j*a for j, a in zip(k, v[:i+1])]) * timeDelta
            k.append(function(x + delta))
        d = timeDelta * sum([b * j for b, j in zip(weights, k)])
        self.pList = (x + d).pList

    def RK4(self, timeDelta):
        """
        Uses the fourth-order Runge-Kutta method to update the internal state
        of the system by `timeDelta` units of time.

        Parameters
        ----------
        timeDelta : float or int
            The time step for the RK4 method.

        Returns
        -------
        None

        Raises
        ------
        TypeError
            If `timeDelta` is not a float.
        ValueError
            If the Runge-Kutta matrix or weights are invalid.

        Notes
        -----
        This method calls the `RK` method with the appropriate Runge-Kutta
        matrix and weights for the fourth-order method. The internal state of
        the system is updated in-place.
        """
        RKmatrix = np.array([[0.5,0,0], [0,0.5,0], [0,0,1]])
        weights = np.array([1/6, 1/3, 1/3, 1/6])
        self.RK(timeDelta, weights, RKmatrix, self._Func)

    def fixCoords(self):
        """
        Fixes the coordinates and coordinate system of each Pendulum object in
        the Chain. If any Pendulum object changes coordinate system, returns
        True

        Returns
        -------
            True if one Pendulum changed coordinate systems, False otherwise.

        Notes
        -----
        Each Pendulum object's phi and theta values are corrected such that
        they fit in the ranges [0, 2pi) and [0, pi) respectively, and if the
        pendulum is too close to a pole in its coordinate system, the
        coordinate system is swapped to the other coordinate system, with the
        position and velocity values properly updated.
        """
        pList = self.pList
        val = False
        for p in pList:
            vi = p.fixCoord()
            val = val or vi
        self.pList = pList
        return val

    def KE(self):
        pList = self.pList
        massList = [p.mass for p in pList]
        xdots = [np.array([0, 0, 0])]
        for p in pList:
            xdots.append(p.Velocity()+xdots[-1])
        T = 0
        for i, xdot in enumerate(xdots[1:]):
            T += 0.5*massList[i]*xdot@xdot
        return T

    def PE(self):
        pList = self.pList
        g = self.g
        massList = [p.mass for p in pList]
        xs = [0]
        for p in pList:
            xs.append(p.toCartesian()[2]+p.length+xs[-1])
        V = 0
        for i, x in enumerate(xs[1:]):
            V += massList[i]*g*x
        return V
