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


def add(crd_1, crd_2, p_corps: int = 7):
    lamb = (crd_2[1] - crd_1[1]) * inversion_mult(crd_2[0] - crd_1[0])
    x3 = (pow(lamb, 2) - crd_1[0] - crd_2[0]) % p_corps
    return [x3, (lamb * (crd_1[0] - x3) - crd_1[1]) % p_corps]


def double(crd, p_corps: int = 7):
    lamb = (3 * pow(crd[0], 2) + A) * inversion_mult(2 * crd[1])
    x3 = (pow(lamb, 2) - 2 * crd[0]) % p_corps
    return [x3, (lamb * (crd[0] - x3) - crd[1]) % p_corps]


def inversion_mult(nb, p_corps: int = 7):
    return pow(nb, p_corps-2) % p_corps


def square_multiply(a, k, p):
    res = 1
    # size = sys.getsizeof(k) * 8
    size = len(bin(k)[2:])
    for i in range(size + 1):
        res = (res ** 2) * a ** ((k >> size - i) & 0b1)
    return res % p
