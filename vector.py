from math import pi, acos, sqrt

class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = "Cannot normalize the zero vector"
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = "Cannot compute the component parallel with the zero vector"
    NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG = "Cannot compute the component orthogonal with the zero vector"
    CANNOT_COMPUTE_ANGLE_MSG = "Cannot compute an angle with the zero vector"
    ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG = "Only can compute in 2D or 3D vectors"
    INFINITE_SMALL = 1e-10

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(coordinates)
            self.dimension = len(self.coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def plus(self, v):
        new_coordinates = [x+y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def minus(self, v):
        new_coordinates = [x-y for x,y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def times_scalar(self, c):
        new_coordinates = [c*x for x in self.coordinates]
        return Vector(new_coordinates)

    def magnitude(self):
        coordinates_squared = [x**2 for x in self.coordinates]
        return sqrt(sum(coordinates_squared))

    def nomalized(self):
        try:
            magnitude = self.magnitude()
            return self.times_scalar(1./magnitude)

        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot(self, v):
        return sum([x*y for x,y in zip(self.coordinates, v.coordinates)])

    def angle(self, v, in_degrees = False):
        try:
            u1 = self.nomalized()
            u2 = v.nomalized()
            angle_in_radians = acos(u1.dot(u2))

            if in_degrees:
                degrees_per_radian = 180. / pi
                return angle_in_radians * degrees_per_radian
            else:
                return angle_in_radians

        except Exception as err:
            if str(err) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.CANNOT_COMPUTE_ANGLE_MSG)

    def parallel(self, v):
        return (self.zero() or
                v.zero() or
                self.angle(v) == 0 or
                self.angle(v) == pi or
                self.angle(v) == None)

    def orthogonal(self, v):
        return abs(self.dot(v)) < INFINITE_SMALL

    def zero(self):
        return self.magnitude() < INFINITE_SMALL

    def component_parallel(self, basis):
        try:
            u = basis.nomalized()
            weight = self.dot(u)
            return u.times_scalar(weight)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)

    def component_orthogonal(self, basis):
        try:
            projection = self.component_parallel(basis)
            return self.minus(projection)

        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)

    def cross(self, v):
        try:
            x1, y1, z1 = self.coordinates
            x2, y2, z2 = v.coordinates
            new_coordinates = [y1*z2 - y2*z1,
                            -(x1*z2 - x2*z1),
                                x1*y2 - x2*y1]
            return Vector(new_coordinates)

        except Exception as e:
            if str(e) == 'need more than 2 values to unpack':
                self_embedded_in_R3 = Vector(self.coordinates + ('0',))
                v_embedded_in_R3 = Vector(v.coordinates + ('0',))
                return self_embedded_in_R3.cross(v_embedded_in_R3)
            elif(str(e) == 'too many values to unpack' or
                str(e) == 'need more than 1 value to unpack'):
                raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)

    def area_parallelogram(self, v):
        cross_product = self.cross(v)
        return cross_product.magnitude()

    def area_triangle(self, v):
        return (self.area_parallelogram(v) / 2.)

