A = 0
B = 3

class Point_aff:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __str__(self):
        return "(%s, %s)" % (self.x, self.y)

    def aff(self):
        print(self.x, self.y)

    # def to_jacob(self) -> Point_jacob:
    #     return Point_jacob(self.x, self.y, 1)

    def e_courbe_aff(self, a):
        return self.y ** 2 == self.x ** 3 + a * self.x


def add(p1, p2, p_corps: int = 7):
    lamb = (p2.y - p1.y) * inversion_mult(p2.x - p1.x)
    x3 = (pow(lamb, 2, p_corps) - p1.x - p2.x) % p_corps
    return Point_aff(x3, (lamb * (p1.x - x3) - p1.y) % p_corps)


def double(p, p_corps: int = 7):
    lamb = (3 * pow(p.x, 2) + A) * inversion_mult(2 * p.y)
    x3 = (pow(lamb, 2) - 2 * p.x) % p_corps
    return Point_aff(x3, (lamb * (p.x - x3) - p.y) % p_corps)


def inversion_mult(nb, p_corps: int = 7):
    return pow(nb, p_corps-2, p_corps)


def square_multiply(a, k, p):
    res = a
    size = len(bin(k)[2:])
    for i in range(size - 2, -1, -1):
        res = (((res ** 2) % p) * a ** ((k >> i) & 0b1)) % p
    return res
