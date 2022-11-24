from jacobien import *
from random import randint


def test_inversion_aff():
    print(inversion_mult(3, 7))


def test_add_aff():
    p1 = Point_aff(5, 3)
    p2 = Point_aff(1, 5)
    print(add(p1, p2))


def test_double_aff():
    p = Point_aff(1, 5)
    print(double(p))


def test_square_multiply():
    print(square_multiply(5, 10, 7))

def test_add_jacob():
    p1 = Point_jacob(5, 3, 1)
    p3 = Point_jacob(1, 5, 1)
    print("Test addition")
    print(p1.add_jacob(p3))
    print(p3.add_jacob(p3))


def test_mult_scal():
    p3 = Point_jacob(1, 5, 1)
    print("Test mult_scal")
    print(p3.mult_scal(3).to_affine())
    print(p3.add_jacob(p3).add_jacob(p3).to_affine())


def test_montgomery():
    p3 = Point_jacob(1, 5, 1)
    print("Test montgomery")
    print(p3.mult_scal(3).to_affine())
    print(p3.montgomery(3).to_affine())


def alice_bob(p: Point_aff, ordre: int):
    print("Alice et Bob")
    a = randint(1, ordre)
    b = randint(1, ordre)

    p_j = Point_jacob(p.x, p.y, 1)

    A = p_j.montgomery(a)
    B = p_j.montgomery(b)

    print(B.montgomery(a).to_affine())
    print(A.montgomery(b).to_affine())


if __name__ == "__main__":
    # test_inversion_aff()
    # test_add_aff()
    # test_double_aff()
    test_square_multiply()

    # test_add_jacob()
    # test_mult_scal()
    # test_montgomery()
    #
    # p4_a = Point_aff(5, 3)
    # alice_bob(p4_a, 13)
