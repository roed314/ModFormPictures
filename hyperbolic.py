# This file contains different models for the hyperbolic plane

from sage.structure.unique_representation import UniqueRepresentation

class HyperbolicPlane(UniqueRepresentation):
    pass

class UpperHalfPlane(HyperbolicPlane):
    xrange = (-1,1)
    yrange = (0,2)
    disk = False
    def convert(self, z, model=None):
        if model is None or model is self:
            return z
        if model is PoincareDisk():
            raise NotImplementedError
        if model is KleinDisk():
            raise NotImplementedError
        raise TypeError("Unsupported model")

class PoincareDisk(HyperbolicPlane):
    xrange = (-1,1)
    yrange = (-1,1)
    disk = True
    def convert(self, z, model=None):
        """
        INPUT:

        - ``z`` -- a complex number representing a point in Poincare disk
        - ``model`` -- a HyperbolicPlane object

        OUTPUT:

        The coordinates of `z` in the given model
        """
        if model is self:
            return z
        if model is None or model is UpperHalfPlane():
            return (1 - 1j*z)/(z - 1j)
        if model is KleinDisk():
            r = z.abs()
            return 2 * z / (1 + r**2)
        raise TypeError("Unsupported model")

class KleinDisk(HyperbolicPlane):
    xrange = (-1,1)
    yrange = (-1,1)
    disk = True
    def convert(self, z, model=None):
        if model is self:
            return z
        if model is None or model is UpperHalfPlane():
            raise NotImplementedError
        if model is PoincareDisk():
            r = z.abs()
            return z / (1 + (1 - r^2).sqrt())
        raise TypeError("Unsupported model")


def hyperbolic_model(name):
    if name == "halfplane":
        return UpperHalfPlane()
    elif name == "poincare":
        return PoincareDisk()
    elif name == "klein":
        return KleinDisk()
    raise ValueError("name must be halfplane, poincare or klein")
