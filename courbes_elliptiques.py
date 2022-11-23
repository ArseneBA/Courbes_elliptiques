from jacobien import *


if __name__ == "__main__":
    p1 = Point_jacob(5, 3, 1)
    p2 = Point_jacob()
    p3 = Point_jacob(1, 5, 1)
    # print(e_courbe_jacob(add_jacob(p1, p3, 7), 0, 3, 7))
    # print(e_courbe_jacob(p1, 0, 3, 7))
    # print(add_jacob(p1, p3, 7))
    # print(add_jacob(p3, p1, 7))
    # print(mult_scal(p1, 7, 7))
    # print(montgomery(p1, 7, 7))
    #
    print(p3.mult_scal(3, 7).to_affine())
    print(p3.montgomery(3, 7).to_affine())
    #
    # print(montgomery(p3, 2, 7))

    # print(homogene(add_jacob(p3, add_jacob(p3, p3))))
    # print(homogene(add_jacob(add_jacob(p3, p3), p3)))
    # print(add_jacob(p1, p3, 7))
    # print(homogene(add_jacob(p1, p3, 7)))
