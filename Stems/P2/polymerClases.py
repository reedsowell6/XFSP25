# polymerClases.py

import math
import random as rnd
from copy import deepcopy as dc

class Position:
    """
    by Jim Smay, last edit: 04/27/2022
    A point in 3D space with vector arithmetic support.
    """
    def __init__(self, pos=None, x=None, y=None, z=None):
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        if pos is not None:
            self.x, self.y, self.z = pos
        self.x = x if x is not None else self.x
        self.y = y if y is not None else self.y
        self.z = z if z is not None else self.z

    def __add__(self, other):
        return Position((self.x + other.x,
                         self.y + other.y,
                         self.z + other.z))

    def __iadd__(self, other):
        if isinstance(other, (int, float)):
            self.x += other
            self.y += other
            self.z += other
        elif isinstance(other, Position):
            self.x += other.x
            self.y += other.y
            self.z += other.z
        return self

    def __sub__(self, other):
        return Position((self.x - other.x,
                         self.y - other.y,
                         self.z - other.z))

    def __isub__(self, other):
        if isinstance(other, (int, float)):
            self.x -= other
            self.y -= other
            self.z -= other
        elif isinstance(other, Position):
            self.x -= other.x
            self.y -= other.y
            self.z -= other.z
        return self

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Position((self.x * other,
                             self.y * other,
                             self.z * other))
        if isinstance(other, Position):
            return Position((self.x * other.x,
                             self.y * other.y,
                             self.z * other.z))

    def __rmul__(self, other):
        return self * other

    def __imul__(self, other):
        if isinstance(other, (int, float)):
            self.x *= other
            self.y *= other
            self.z *= other
        return self

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return Position((self.x / other,
                             self.y / other,
                             self.z / other))

    def __idiv__(self, other):
        if isinstance(other, (int, float)):
            self.x /= other
            self.y /= other
            self.z /= other
        return self

    def __round__(self, n=None):
        if n is not None:
            return Position(x=round(self.x, n),
                            y=round(self.y, n),
                            z=round(self.z, n))
        return self

    def mag(self):
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)

    def normalize(self):
        length = self.mag()
        if length > 0:
            self.__idiv__(length)

    def getRndDir(self):
        # No per-call reseeding—uses global rnd state for reproducibility
        d = Position(x=rnd.random(), y=rnd.random(), z=rnd.random())
        d -= 0.5
        d.normalize()
        return d

    def getRndPosOnSphere(self, radius=1.0):
        return self + self.getRndDir() * radius

    def distTo(self, other=None):
        if other is None:
            return self.mag()
        return (self - other).mag()

class molecule:
    """
    A single mer (CH2 unit) with a mass and position.
    """
    def __init__(self, molecularWeight=12, position=None):
        self.MW = molecularWeight
        self.position = position if position is not None else Position()

class macroMolecule:
    """
    A polymer chain modeled as a freely‐jointed chain of mers.
    """
    def __init__(self,
                 targetDegreeOfPolymerization=1000,
                 segmentLength=0.154E-9,
                 merWt=14):
        # sample actual N from N(mean, 0.1·mean)
        actualN = int(rnd.gauss(targetDegreeOfPolymerization,
                                 0.1 * targetDegreeOfPolymerization))
        self.N = max(1, actualN)
        self.merWt = merWt
        self.MW = self.N * merWt
        self.segmentLength = segmentLength
        self.mers = []
        self.centerOfMass = Position()
        self.radiusOfGyration = 0.0
        self.endToEndDistance = 0.0

    def freelyJointedChainModel(self):
        # start at origin
        lastPos = Position(x=0, y=0, z=0)
        self.mers = []

        # build chain
        for i in range(self.N):
            m = molecule(molecularWeight=self.merWt)
            # adjust end hydrogens
            if i == 0 or i == self.N - 1:
                m.MW += 1
            m.position = lastPos.getRndPosOnSphere(self.segmentLength)
            self.mers.append(m)
            lastPos = m.position

        # center of mass
        total_mass = sum(mer.MW for mer in self.mers)
        com = Position(x=0, y=0, z=0)
        for mer in self.mers:
            com += mer.MW * mer.position
        self.centerOfMass = com / total_mass

        # end‐to‐end distance
        self.endToEndDistance = (
            self.mers[0].position - self.mers[-1].position
        ).mag()

        # radius of gyration
        rg2 = sum(
            mer.MW * (mer.position.distTo(self.centerOfMass)**2)
            for mer in self.mers
        ) / total_mass
        self.radiusOfGyration = math.sqrt(rg2)
