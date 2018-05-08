from numpy import loadtxt

def load(filename):
    lines = loadtxt(filename, comments="#", delimiter=" ", unpack=False)
    x = [float(n) for n in lines[:, 0]]
    y = [float(n) for n in lines[:, 1]]
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

def inertia(x, y, onlyX=False, onlyY=False, onlyXY=False, debug=False):
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
    if onlyX: return ICx
    if onlyY: return ICy
    if onlyXY: return ICxy
    if not debug:
        return ICx, ICy, ICxy
    else: return ICx, ICy, ICxy, xC, yC, area(x, y, signed=False)


x, y = load("blade")
ICx, ICy, ICxy, xC, yC, A = inertia(x, y, debug=True)
print("x = ")
print(x)
print("y = ")
print(y)
print("xC = ")
print(xC)
print("yC = ")
print(yC)
print("ICx = ")
print(ICx)
print("ICy = ")
print(ICy)
print("ICxy = ")
print(ICxy)