from jacobien import *
from random import randint


def alice_bob(p: Point_aff, ordre: int):
    a = randint(1, ordre)
    b = randint(1, ordre)

    p_j = Point_jacob(p.x, p.y, 1)

    A = p_j.montgomery(a)
    B = p_j.montgomery(b)

    print(B.montgomery(a).to_affine())
    print(A.montgomery(b).to_affine())


if __name__ == "__main__":
    # p1 = Point_jacob(5, 3, 1)
    # p2 = Point_jacob()
    # p3 = Point_jacob(1, 5, 1)
    # print(e_courbe_jacob(add_jacob(p1, p3, 7), 0, 3, 7))
    # print(e_courbe_jacob(p1, 0, 3, 7))
    # print(add_jacob(p1, p3, 7))
    # print(add_jacob(p3, p1, 7))
    # print(mult_scal(p1, 7, 7))
    # print(montgomery(p1, 7, 7))
    #
    # print(p3.mult_scal(3, 7).to_affine())
    # print(p3.montgomery(3, 7).to_affine())
    #
    # print(montgomery(p3, 2, 7))

    # print(homogene(add_jacob(p3, add_jacob(p3, p3))))
    # print(homogene(add_jacob(add_jacob(p3, p3), p3)))
    # print(add_jacob(p1, p3, 7))
    # print(homogene(add_jacob(p1, p3, 7)))

    p4_a = Point_aff(5, 3)
    alice_bob(p4_a, 13)
