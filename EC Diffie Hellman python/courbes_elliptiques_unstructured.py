from random import randint

A = 0
B = 3
CORPS = 7


def add(crd_1, crd_2, p_corps: int = 7):
    lamb = (crd_2[1] - crd_1[1]) * inversion_mult(crd_2[0] - crd_1[0])
    x3 = (pow(lamb, 2) - crd_1[0] - crd_2[0]) % p_corps
    return [x3, (lamb * (crd_1[0] - x3) - crd_1[1]) % p_corps]


def double(crd, p_corps: int = 7):
    lamb = (3 * pow(crd[0], 2) + A) * inversion_mult(2 * crd[1])
    x3 = (pow(lamb, 2) - 2 * crd[0]) % p_corps
    return [x3, (lamb * (crd[0] - x3) - crd[1]) % p_corps]


def inversion_mult(nb, p_corps: int = 7):
    return pow(nb, p_corps - 2) % p_corps


class Point_jacob:
    def __init__(self, x=1, y=1, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "(%s, %s, %s)" % (self.x, self.y, self.z)

    def aff(self):
        print(self.x, self.y, self.z)


class Point_aff:
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __str__(self):
        return "(%s, %s)" % (self.x, self.y)

    def aff(self):
        print(self.x, self.y)


def add_jacob(p1: Point_jacob, p2: Point_jacob, p_corps: int = 7) -> Point_jacob:
    p3 = Point_jacob()
    if p1 == p2:
        x1 = p1.x
        y1 = p1.y
        z1 = p1.z

        A = x1 ** 2
        B = y1 ** 2
        C = B ** 2
        D = 2 * ((x1 + B) ** 2 - A - C)
        E = 3 * A
        F = E ** 2
        p3.x = (F - 2 * D) % p_corps
        p3.y = (E * (D - p3.x) - 8 * C) % p_corps
        p3.z = (2 * y1 * z1) % p_corps

    else:
        x1 = p1.x
        y1 = p1.y
        z1 = p1.z

        x2 = p2.x
        y2 = p2.y
        z2 = p2.z

        u1 = x1 * z2 ** 2
        u2 = x2 * z1 ** 2
        s1 = y1 * z2 * z2 ** 2
        s2 = y2 * z1 ** 3
        h = u2 - u1
        i = (2 * h) ** 2
        j = h * i
        r = 2 * (s2 - s1)
        v = u1 * i

        p3.x = (r ** 2 - j - 2 * v) % p_corps
        p3.y = (r * (v - p3.x) - 2 * s1 * j) % p_corps
        p3.z = (((z1 + z2) ** 2 - z1 * z1 - z2 * z2) * h) % p_corps

    return p3


def square_multiply(a, k, p):
    res = 1
    # size = sys.getsizeof(k) * 8
    size = len(bin(k)[2:])
    for i in range(size + 1):
        res = (res ** 2) * a ** ((k >> size - i) & 0b1)
    return res % p


# Problème -> on voit le if sur la consommation (une étape en plus)
# Solution naïve : faire une multiplication en plus (idem que si la valeur du bit vaut 1) et la placer dans une poubelle


def mult_scal(p: Point_jacob, k, p_corps: int = 7):
    res = p
    n = len(bin(k)[2:])
    for i in range(n - 1, 0, -1):
        res = add_jacob(res, res, p_corps)
        if (k >> i) & 0b1:
            res = add_jacob(res, p, p_corps)
    return res


def e_courbe_jacob(p: Point_jacob, a, b, p_corps: int = 7):
    return pow(p.y, 2, p_corps) == (
                pow(p.x, 3, p_corps) + a * p.x * pow(p.z, 4, p_corps) + b * pow(p.z, 6, p_corps)) % p_corps


def montgomery(p: Point_jacob, k: int, p_corps: int = 7) -> Point_jacob:
    # Version plus intelligente du double and add
    p_res_0 = p
    p_res_1 = add_jacob(p, p, p_corps)
    n = len(bin(k)[2:])
    for i in range(n - 1, 0, -1):
        if ((k >> i) & 0b1) == 1:
            p_res_0 = add_jacob(p_res_0, p_res_1, p_corps)
            p_res_1 = add_jacob(p_res_1, p_res_1, p_corps)
        else:
            p_res_1 = add_jacob(p_res_1, p_res_0, p_corps)
            p_res_0 = add_jacob(p_res_0, p_res_0, p_corps)
    return p_res_0


def e_courbe_aff(p: Point_aff, a):
    return p.y ** 2 == p.x ** 3 + a * p.x


def homogene(p: Point_jacob, p_corps: int = 7) -> Point_aff:
    return Point_aff((p.x * inversion_mult(p.z) ** 2) % p_corps, (p.y * inversion_mult(p.z) ** 3) % p_corps)


def alice_bob(p: Point_aff, ordre: int):
    a = randint(1, ordre)
    b = randint(1, ordre)

    p_j = Point_jacob(p.x, p.y, 1)

    A = montgomery(p_j, a)
    B = montgomery(p_j, b)

    print(homogene(montgomery(B, a)))
    print(homogene(montgomery(A, b)))

if __name__ == "__main__":
    p1 = Point_jacob(5, 3, 1)
    p2 = Point_jacob()
    p3 = Point_jacob(1, 5, 1)
    print(e_courbe_jacob(add_jacob(p1, p3, 7), 0, 3, 7))
    print(e_courbe_jacob(p1, 0, 3, 7))
    print(add_jacob(p1, p3, 7))
    print(add_jacob(p3, p1, 7))
    print(mult_scal(p1, 7, 7))
    print(montgomery(p1, 7, 7))

    print(homogene(mult_scal(p3, 3, 7)))
    print(homogene(montgomery(p3, 3, 7)))

    print(montgomery(p3, 2, 7))

    print(homogene(add_jacob(p3, add_jacob(p3, p3))))
    print(homogene(add_jacob(add_jacob(p3, p3), p3)))
    print(add_jacob(p1, p3, 7))
    print(homogene(add_jacob(p1, p3, 7)))

    p4_a = Point_aff(5, 3)
    alice_bob(p4_a, 13)

