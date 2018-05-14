import math

def load(filename, usenumpy=False, xindex=0, yindex=1):
    if usenumpy:
        from numpy import loadtxt
        lines = loadtxt(filename, comments="#", delimiter=" ", unpack=False)
        x = [float(n) for n in lines[:, xindex]]
        y = [float(n) for n in lines[:, yindex]]
        return x, y
    else:
        lines = open(filename, "r").readlines()
        x = [float(n.split(' ')[xindex]) for n in lines if not n.startswith("#")]
        y = [float(n.split(' ')[yindex]) for n in lines if not n.startswith("#")]
        return x, y

def array_delta(array, delta):
    array[:] = [x + delta for x in array]
    return array

def area(x, y, signed=False):
    # https://en.wikipedia.org/wiki/Centroid
    A = 0
    # the last point must be the same as the first
    x.append(x[0])
    y.append(y[0])
    for i in range(0, len(x) - 1):
        A += x[i] * y[i+1] - x[i+1] * y[i]
    if signed:
        return A / 2
    else: return  abs(A / 2)

def centroid(x, y):
    # https://en.wikipedia.org/wiki/Centroid
    xC, yC = 0, 0
    A = area(x, y, signed=True)
    # the last point must be the same as the first
    x.append(x[0])
    y.append(y[0])
    for i in range(0, len(x) - 1):
        # TODO extract common term of Area
        xC += (x[i] + x[i + 1]) * (x[i] * y[i + 1] - x[i+1] * y[i])
        yC += (y[i] + y[i + 1]) * (x[i] * y[i + 1] - x[i+1] * y[i])
    return xC / 6 / A, yC / 6 / A

def solvequadratic(a, b, c):
    if a != 0:
        x1 = (-b + math.sqrt(b**2 - 4 * a * c)) / 2 / a
        x2 = (-b - math.sqrt(b**2 - 4 * a * c)) / 2 / a
        return x1, x2
    else:
        return -0.5*b/a, -0.5*b/a

def eig(matrix):
    a = matrix[0][0]
    b = matrix[0][1]
    c = matrix[1][0]
    d = matrix[1][1]
    lambda1, lambda2 = solvequadratic(1,-a-d, a*d - c*b)
    # TODO fix when c = 0 and in any particular case
    if abs(lambda1-a)>1e-10:
        vector1 = [1, b / (lambda1-a)]
    else: vector1 = [1, (lambda1-d)/c]

    if abs(lambda2-a)>1e-10:
        vector2 = [1, b / (lambda2-a)]
    else: vector2 = [1, (lambda2-d)/c]

    return lambda1, lambda2, vector1, vector2


def inertia(x, y, principal=False, onlyX=False, onlyY=False, onlyXY=False, debug=False, tupleInsteadOfTensor=False):
    # https://en.wikipedia.org/wiki/Second_moment_of_area#Any_cross_section_defined_as_polygon
    xC, yC = centroid(x, y)
    ICx, ICy, ICxy = 0, 0, 0
    # the last point must be the same as the first
    x.append(x[0])
    y.append(y[0])
    x = array_delta(x, -xC)
    y = array_delta(y, -yC)
    for i in range(0, len(x) - 1):
        # TODO extract common term of Area
        ICx += (y[i]**2 + y[i]*y[i + 1] + y[i + 1]**2) * (x[i] * y[i + 1] - x[i+1] * y[i])
        ICy += (x[i]**2 + x[i]*x[i + 1] + x[i + 1]**2) * (x[i] * y[i + 1] - x[i+1] * y[i])
        ICxy += (x[i]*y[i+1] + 2*x[i]*y[i] + 2*x[i+1]*y[i+1] + x[i+1]*y[i]) * (x[i] * y[i + 1] - x[i+1] * y[i])
    ICx /= 12
    ICy /= 12
    ICxy /= 24
    ICz = ICx + ICy
    if onlyX: return ICx
    if onlyY: return ICy
    if onlyXY: return ICxy
    I = [[ICx, ICxy, 0], [ICxy, ICy, 0], [0, 0, ICz]]
    I1, I2, v1, v2 = eig([[ICx, ICxy], [ICxy, ICy]])
    direction = math.atan(v1[1] / v1[0])
    if principal:
        if not debug:
            return [I1, I2, ICz], direction
        else:
            return [I1, I2, ICz], direction, ICx, ICy, ICxy, ICz, xC, yC, area(x, y, signed=False)
    else:
        if not debug:
            if tupleInsteadOfTensor:
                return ICx, ICy, ICxy, ICz
            return I
        else:
            if tupleInsteadOfTensor:
                ICx, ICy, ICxy, ICz, xC, yC, area(x, y, signed=False)
            return I, xC, yC, area(x, y, signed=False)


x, y = load("blade")
#I, xC, yC, A = inertia(x, y, debug=True, principal=True)
Iprincipal, direction, ICx, ICy, ICxy, ICz, xC, yC, A = inertia(x, y, debug=True, principal=True)
#print("x = ")
#print(x)
#print("y = ")
#print(y)
print("xC = ")
print(xC)
print("yC = ")
print(yC)
print('A = ')
print(A)
#print('I (m^4) = ')
#print(I)
print('I principal = ')
print(Iprincipal)
print('Direction = (deg)')
print(direction / math.pi * 180)
#print("ICx = ")
#print(ICx)
#print("ICy = ")
#print(ICy)
#print("ICxy = ")
#print(ICxy)