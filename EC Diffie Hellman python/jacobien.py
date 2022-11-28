from affine import *


class Point_jacob:
    pass


class Point_jacob:
    def __init__(self, x=1, y=1, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "(%s, %s, %s)" % (self.x, self.y, self.z)

    def aff(self):
        print(self.x, self.y, self.z)

    def to_affine(self, p_corps: int = 7) -> Point_aff:
        return Point_aff((self.x * pow(inversion_mult(self.z), 2, p_corps)) % p_corps,
                         (self.y * pow(inversion_mult(self.z), 3, p_corps)) % p_corps)

    def e_courbe_jacob(self, a, b, p_corps: int = 7):
        p = self
        return pow(p.y, 2, p_corps) == (
                    pow(p.x, 3, p_corps) + a * p.x * pow(p.z, 4, p_corps) + b * pow(p.z, 6, p_corps)) % p_corps

    def add_jacob(self, p2, p_corps: int = 7) -> Point_jacob:
        p1 = self
        p3 = Point_jacob()

        if p1 == p2:
            x1 = p1.x
            y1 = p1.y
            z1 = p1.z

            A = pow(x1, 2, p_corps)
            B = pow(y1, 2, p_corps)
            C = pow(B, 2, p_corps)
            D = (2 * (pow((x1 + B), 2, p_corps) - A - C)) % p_corps
            E = (3 * A) % p_corps
            F = pow(E, 2, p_corps)
            p3.x = (F - (2 * D) % p_corps) % p_corps
            p3.y = ((E * (D - p3.x)) % p_corps - (8 * C) % p_corps) % p_corps
            p3.z = (((2 * y1) % p_corps) * z1) % p_corps

        else:
            x1 = p1.x
            y1 = p1.y
            z1 = p1.z

            x2 = p2.x
            y2 = p2.y
            z2 = p2.z

            u1 = (x1 * pow(z2, 2, p_corps)) % p_corps
            u2 = (x2 * pow(z1, 2, p_corps)) % p_corps
            s1 = (((y1 * z2) % p_corps) * pow(z2, 2, p_corps)) % p_corps
            s2 = (y2 * pow(z1, 3, p_corps)) % p_corps
            h = u2 - u1
            i = pow((2 * h), 2, p_corps)
            j = (h * i) % p_corps
            r = (2 * (s2 - s1)) % p_corps
            v = (u1 * i) % p_corps

            p3.x = (pow(r, 2, p_corps) - j - (2 * v) % p_corps) % p_corps
            p3.y = ((r * (v - p3.x)) % p_corps - (((2 * s1) % p_corps) * j) % p_corps) % p_corps
            p3.z = ((pow((z1 + z2), 2, p_corps) - (z1 * z1) % p_corps - (z2 * z2) % p_corps) * h) % p_corps

        return p3

    def mult_scal(self, k, p_corps: int = 7) -> Point_jacob:
        # double add
        res = self
        n = len(bin(k)[2:])
        for i in range(n - 2, -1, -1):
            res = res.add_jacob(res, p_corps)
            if (k >> i) & 0b1:
                res = res.add_jacob(self, p_corps)
        return res

    # Problème -> on voit le if sur la consommation (une étape en plus)
    # Solution naïve : faire une multiplication en plus
    # (idem que si la valeur du bit vaut 1) et la placer dans une poubelle

    def montgomery(self, k: int, p_corps: int = 7) -> Point_jacob:
        # Version plus intelligente du double and add
        p_res_0 = self
        p_res_1 = self.add_jacob(self, p_corps)
        n = len(bin(k)[2:])
        for i in range(n-2, -1, -1):
            if ((k >> i) & 0b1) == 1:
                p_res_0 = p_res_0.add_jacob(p_res_1, p_corps)
                p_res_1 = p_res_1.add_jacob(p_res_1, p_corps)
            else:
                p_res_1 = p_res_1.add_jacob(p_res_0, p_corps)
                p_res_0 = p_res_0.add_jacob(p_res_0, p_corps)
        return p_res_0



