import numpy as np

class Pendulum():
    """Pendulum(length, pos, vel, mass, zPolar=True)
    A class that represents a spherical pendulum.

    Parameters
    ----------
    length : float or int
        The length of the pendulum.
    pos : array-like of shape (2,)
        The position of the pendulum in spherical coordinates.
    vel : array-like of shape (2,)
        The velocity of the pendulum in spherical coordinates.
    mass : float or int
        The mass of the pendulum.
    zPolar : bool, optional
        The coordinate system to use. If True, use z-polar coordinates.
        If False, use x-polar coordinates. Default is True.

    Attributes
    ----------
    length : numpy 64 bit float
        The length of the pendulum.
    pos : numpy array (2,)
        The angular position of the pendulum in spherical coordinates
        (theta, phi).
    vel : numpy array (2,)
        The angular velocity of the pendulum in spherical coordinates
        (dtheta, dphi).
    mass : numpy 64 bit float
        The mass of the pendulum.
    zPolar : bool
        The coordinate system being used.

    Raises
    ------
    TypeError
        If pos or vel is not array-like, or if mass or length is not
        float-like.
    ValueError
        If pos or vel does not have a length of 2, or if mass or length is not
        positive and non-zero.

    Methods
    -------
    dT(self)
        Calculate the partial derivative of the cartesian representation of
        the pendulum with respect to theta (the polar angle).
    dP(self)
        Calculate the partial derivative of the cartesian representation of
        the pendulum with respect to phi (the azimuthal angle).
    dT2(self)
        Calculate the second partial derivative of the cartesian representation
        of the pendulum with respect to theta.
    dP2(self)
        Calculate the second partial derivative of the cartesian representation
        of the pendulum with respect to phi.
    dTP(self)
        Calculate the partial derivative of the cartesian representation of
        the pendulum with respect to both theta and phi.
    drdT(self)
        Calculate the partial derivative of the z element of the cartesian
        representation of the pendulum with respect to theta.
    drdP(self)
        Calculate the partial derivative of the z element of the cartesian
        representation of the pendulum with respect to phi.
    Velocity(self)
        Calculate the cartesian velocity of the pendulum.
    toCartesian(self)
        Calculate the cartesian position of the pendulum.
    fixCoord(self, release=True)
        Swap the coordinate system if the pendulum gets too close to a pole.

    Dunder Methods
    --------------
    __mul__(self, other)
        Multiply the pendulum by a scalar.
    __rmul__(self, other)
        Multiply the pendulum by a scalar.
    __add__(self, other)
        Add the pendulum to another pendulum object.
    """

    def __init__(self, length, pos, vel, mass, zPolar=True):
        """
        Initialize a spherical pendulum object.

        Parameters
        ----------
        length : float or int
            Length of the pendulum.
        pos : array-like of shape (2,)
            Initial position of the pendulum, specified as a pair of polar
            angles (theta, phi).
        vel : array-like of shape (2,)
            Initial velocity of the pendulum, specified as a pair of angular
            velocities (dtheta, dphi).
        mass : float or int
            Mass of the pendulum.
        zPolar : bool, optional
            If True, use the z-polar coordinate system (default). If False, use
            the x-polar coordinate system.

        Raises
        ------
        TypeError
            If pos or vel are not array-like, or if mass or length are not
            float-like.
        ValueError
            If pos or vel do not have length 2, or if mass or length are not
            positive and non-zero.
        """
        try:
            self.pos = np.array(pos)
            self.vel = np.array(vel)
        except TypeError:
            raise TypeError("pos and vel must be array-like")

        if len(self.pos) != 2:
            raise ValueError("pos must have len == 2")
        if len(self.vel) != 2:
            raise ValueError("vel must have len == 2")

        try:
            self.mass = np.float64(mass)
            self.length = np.float64(length)
        except TypeError:
            raise TypeError("mass and length must be float-like")

        if self.mass <= 0 or self.length <=0:
            raise ValueError("mass and length must be positive and non-zero")

        self.zPolar = bool(zPolar)

    def __str__(self):
        return f"p: pos={self.pos}, vel={self.vel}, mass={self.mass}, \
length={self.length}, basis={(lambda v: 'z' if v else 'x')(self.zPolar)}"

    def __repr__(self):
        return f"Pendulum(length={self.length}, mass={self.mass}, \
pos={self.pos}, vel={self.vel}, zPolar={self.zPolar})"

    def __eq__(self, other):
        if not isinstance(other, Pendulum):
            return False
        else:
            length = self.length == other.length
            mass = self.mass == other.mass
            pos = all([i == j for i, j in zip(self.pos, other.pos)])
            vel = all([i == j for i, j in zip(self.vel, other.vel)])
            pole = self.zPolar == other.zPolar
            val = all([length, mass, pos, vel, pole])
            return val

    def __mul__(self, other):
        """
        Multiply a pendulum object by a scalar value.

        Parameters
        -----------
        other : float or int
            The scalar value to multiply the pendulum object by.

        Returns
        --------
        Pendulum
            A new pendulum object with position and velocity multiplied by the
            scalar value. The mass, length, and zPolar attributes are copied
            from the original pendulum object.

        Raises
        -------
        TypeError
            If the input value is not a float or integer, a TypeError is raised
            with the message "Pendulum can only be multiplied by a scalar".
        """
        try:
            other = np.float32(other)
        except TypeError:
            raise TypeError("Pendulum can only be multiplied by a scalar")
        return Pendulum(pos=self.pos*other,vel=self.vel*other, mass=self.mass,
                        length=self.length, zPolar=self.zPolar)

    def __rmul__(self, other):
        """
        Called when the Pendulum object is on the right side of the
        multiplication operator (*) and the left operand is not a Pendulum
        object.

        Parameters
        ----------
        other : float or int
            The scalar value to multiply the Pendulum object with.

        Returns
        -------
        Pendulum
            A new Pendulum object with its position and velocity scaled by the
            scalar value 'other' while the mass, length, and zPolar remain
            unchanged.

        Raises
        ------
        TypeError
            If the input 'other' is not a scalar value (float or int).
        """
        return self.__mul__(other)

    def __add__(self, other):
        """Add two Pendulum objects and return a new Pendulum object.

        Parameters
        ----------
        other : Pendulum
            The Pendulum object to be added to the current object.

        Returns
        -------
        Pendulum
            A new Pendulum object with the sum of positions and velocities of
            the two objects. The length, mass, and zPolar are taken from the
            first pendulum in the addition.

        Raises
        ------
        TypeError
            If the 'other' parameter is not a Pendulum object.
        """
        if type(other) is not Pendulum:
            raise TypeError(
"Pendulum object can only be added to other Pendulum objects")

        Q = Pendulum(pos=self.pos+other.pos,vel=self.vel+other.vel,
                     mass=self.mass, length=self.length,zPolar=self.zPolar)

        return Q

    def __sub__(self, other):
        if not isinstance(other, Pendulum):
            raise TypeError("Pendulum objects can only be subtracted from \
other Pendulum objects")
        return self + -1*other

    def __abs__(self):
        val = abs(self.pos[0]) + abs(self.vel[0]) + abs(self.pos[1]) +\
            abs(self.vel[1])
        return val

    def dT(self):
        """
        Calculates the partial derivative of the cartesian position of the
        pendulum with respect to theta, the polar angle.

        Returns
        -------
        numpy.ndarray, shape=(3,)
            A 3-element numpy array representing the partial derivative of
            the cartesian position of the pendulum with respect to theta.

        Notes
        -----
        The formula used and the definition of theta depends on the coordinate
        system. If zPolar is True then theta is the polar angle from the
        (negative) z-axis, but if zPolar is False then theta is the polar angle
        from the x-axis. The formula in either case is [dx/dtheta, dy/dtheta,
        dz/dtheta], but the definition of x, y, and z, depends on the spherical
        coordinate basis. The formulae for x, y, and z were chosen such that in
        either coordinate system, they would return the same cartesian values.
        """
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
        """
        Calculates the partial derivative of the cartesian position of the
        pendulum with respect to phi, the azimuthal angle.

        Returns
        -------
        numpy.ndarray, shape=(3,)
            A 3-element numpy array representing the partial derivative of
            the cartesian position of the pendulum with respect to phi.

        Notes
        -----
        The formula used and the definition of phi depends on the coordinate
        system. If zPolar is True then phi is the azimuthal angle from the
        x-axis going to the y-axis as phi increases, but if zPolar is False
        then phi is the azimuthal angle going from the y-axis to the z-axis as
        phi increases. The formula in either case is [dx/dphi, dy/dphi,
        dz/dphi], but the definition of x, y, and z, depends on the spherical
        coordinate basis. The formulae for x, y, and z were chosen such that in
        either coordinate system, they would return the same cartesian values.
        """
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
        """
        Calculates the second partial derivative of the cartesian position of
        the pendulum with respect to theta, the polar angle.

        Returns
        -------
        numpy.ndarray, shape=(3,)
            A 3-element numpy array representing the second partial derivative
            of the cartesian position of the pendulum with respect to theta.

        Notes
        -----
        The formula used and the definition of theta depends on the coordinate
        system. If zPolar is True then theta is the polar angle from the
        (negative) z-axis, but if zPolar is False then theta is the polar angle
        from the x-axis. The formula in either case is [d^2x/dtheta^2,
        d^2y/dtheta^2, d^2z/dtheta^2], but the definition of x, y, and z,
        depends on the spherical coordinate basis. The formulae for x, y, and z
        were chosen such that in either coordinate system, they would return
        the same cartesian values.
        """
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
        """
        Calculates the second partial derivative of the cartesian position of
        the pendulum with respect to phi, the azimuthal angle.

        Returns
        -------
        numpy.ndarray, shape=(3,)
            A 3-element numpy array representing the second partial derivative
            of the cartesian position of the pendulum with respect to phi.

        Notes
        -----
        The formula used and the definition of phi depends on the coordinate
        system. If zPolar is True then phi is the azimuthal angle from the
        x-axis going to the y-axis as phi increases, but if zPolar is False
        then phi is the azimuthal angle going from the y-axis to the z-axis as
        phi increases. The formula in either case is [d^2x/dphi^2, d^2y/dphi^2,
        d^2z/dphi^2], but the definition of x, y, and z, depends on the
        spherical coordinate basis. The formulae for x, y, and z were chosen
        such that in either coordinate system, they would return the same
        cartesian values.
        """
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
        """
        Calculates the second partial derivative of the cartesian position of
        the pendulum with respect to both theta and phi, the polar and
        azimuthal angles respectively.

        Returns
        -------
        numpy.ndarray, shape=(3,)
            A 3-element numpy array representing the second partial derivative
            of the cartesian position of the pendulum with respect to both
            theta and phi.

        Notes
        -----
        The formula used and the definition of theta and phi depends on the
        coordinate basis used. If zPolar is True then theta is the polar angle
        from the (negative) z-axis and phi is the azimuthal angle from the
        x-axis going to the y-axis as phi increases, but if zPolar is False
        then theta is the polar angle from the x-axis and phi is the azimuthal
        angle going from the y-axis to the z-axis as phi increases. The formula
        in either case is [d^2x/dphidtheta, d^2y/dphidtheta, d^2z/dphidtheta],
        but the definition of x, y, and z, depends on the spherical coordinate
        basis. The formulae for x, y, and z were chosen such that in either
        coordinate system, they would return the same cartesian values.
        """
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

    def dzdT(self):
        """Calculates the partial derivative of the z-component of the
        cartesian representation of the position of the pendulum with respect
        to theta, the polar angle.

        Returns
        -------
        numpy 64 bit float
            The value of the partial derivative of the z-component of the
            cartesian representation of the pendulum position with respect to
            theta.

        Notes
        -----
        The formula used and the definition of theta depends on the coordinate
        system. If zPolar is True then theta is the polar angle from the
        (negative) z-axis, but if zPolar is False then theta is the polar angle
        from the x-axis. The formula in either case is dz/dtheta, but the
        definition of z depends on the spherical coordinate basis. The formulae
        for z was chosen such that in either coordinate system, it would return
        the same cartesian value.
        """
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return l * np.sin(theta)
        else:
            return l * np.cos(theta) * np.sin(phi)

    def dzdP(self):
        """Calculates the partial derivative of the z-component of the
        cartesian representation of the position of the pendulum with respect
        to phi, the azimuthal angle.

        Returns
        -------
        float
            The value of the partial derivative of the z-component of the
            cartesian representation of the pendulum position with respect to
            phi.

        Notes
        -----
        The formula used and the definition of theta depends on the coordinate
        system. If zPolar is True then phi is the azimuthal angle from the
        x-axis going to the y-axis as phi increases, but if zPolar is False
        then phi is the azimuthal angle going from the y-axis to the z-axis as
        phi increases. The formula in either case is dz/dphi, but the
        definition of z depends on the spherical coordinate basis. The formulae
        for z was chosen such that in either coordinate system, it would return
        the same cartesian value.
        """
        l = self.length
        theta = self.pos[0]
        phi = self.pos[1]

        if self.zPolar:
            return 0
        else:
            return l * np.sin(theta) * np.cos(phi)

    def Velocity(self):
        """
        Returns the cartesian velocity of the pendulum.

        Returns
        -------
        numpy.ndarray
            A 3 element numpy array representing the cartesian velocity of the
            pendulum.

        Notes
        -----
        The method calculates the partial derivatives of the cartesian position
        of the pendulum with respect to both theta and phi, and scales them by
        their respective angular velocities. It then takes the sum of the
        scaled derivatives to obtain the velocity vector in cartesian
        coordinates.

        This is derived using the chain rule, and thus applies in
        either coordinate systems.
        """
        return self.vel[0]*self.dT() + self.vel[1]*self.dP()

    def toCartesian(self):
        """Calculates the cartesian 3D position of the pendulum.

        Returns
        -------
        numpy.ndarray
            A 3 element numpy array representing the cartesian position of the
            pendulum.

        Notes
        -----
        formulae used depends on which spherical coordinate system the pendulum
        is in.

        One coordinate system (zPolar==True) has the pole be parallel to
        the z-axis such that a value of theta=0 results in z = -length and a
        value of phi=0 will place the pendulum on the plane parallel to the
        x-axis and y-axis.

        The other coordinate system (zPolar==False) has the pole be parallel to
        the x-axis, such that a value of theta=0 results in x=length and a
        value of phi=0 will place the pendulum on the plane parallel to the
        y-axis and x-axis.
        """
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

    def fixCoord(self):
        """
        Adjusts the current coordinate system to ensure that phi (azimuthal
        angle) is within the range [0, 2pi) and that theta (polar angle) is
        within the range [0, pi).
        Then, if the value of sin(theta) is less than 0.5, the angles and
        angular velocities are transformed to the other coordinate system, such
        that the Cartesian representation remains constant. The value of zPolar
        is also changed to match this.

        Returns
        -------
        bool
            True if the coordinate system has been changed, False otherwise.

        Notes
        -----
        This method is used to ensure that the pendulum doesn't get too close
        to a pole. This would be an issue, as in spherical coordinate systems,
        the poles are singularities which would lead to numeric instability in
        the value of phi and the value of phi angular velocity.

        The angular positions are calculated very simply, as in the z-polar
        coordinate system,
        phi = atan2(y, x)
        theta = atan2(sqrt(x^2 + y^2), -z)
        and in the x-polar coordinate system,
        phi = atan2(z, y)
        theta = atan2(sqrt(y^2 + z^2), x)

        For the angular velocities, inner products of the cartesian velocity
        with partial derivatives are used. The cartesian velocity v is given by

        v = thetadot * dq/dtheta + phidot * dq/dtheta

        where q is the cartesian position vector of the pendulum. This applies
        in both spherical coordinate systems as it is derived from the chain
        rule. Using this formulae, if you take its inner product with
        partial derivatives in theta and phi, the following can be derived,

        v dot dq/dtheta = l^2 * thetadot
        v dot dq/dphi   = l^2 * sin^2(theta) * phidot

        Thus if the value of v is calculated before the positions are
        transformed, it can be constant over transformation and the new partial
        derivatives and values of phi and theta can be used to calculate the
        new values of thetadot and phidot such that the cartesian velocity is
        conserved over transformation.
        """

        theta = self.pos[0]
        phi = self.pos[1]

        phi = np.pi * ((phi/np.pi)%2)
        theta = np.pi * ((theta/np.pi)%2)

        if theta > np.pi:
            theta = 2*np.pi - theta
            phi += np.pi
            phi = np.pi * ((phi/np.pi)%2)

        self.pos = np.array([theta, phi])

        if abs(np.sin(theta)) > 0.5:
            return False

        p = self.toCartesian()
        v = self.Velocity()

        if self.zPolar:
            phi2 = np.arctan2(p[2], p[1])
            theta2 = np.arctan2(np.sqrt(p[1]**2+p[2]**2), p[0])

        else:
            phi2 = np.arctan2(p[1], p[0])
            theta2 = np.arctan2(np.sqrt(p[0]**2 + p[1]**2), -p[2])

        self.pos[0] = theta2
        self.pos[1] = phi2

        self.zPolar = not self.zPolar

        self.vel[0] = (1/self.length)**2 * v@self.dT()
        self.vel[1] = (1/(self.length * np.sin(theta2)))**2 * v@self.dP()

        return True

if __name__ == "__main__":
    pass
