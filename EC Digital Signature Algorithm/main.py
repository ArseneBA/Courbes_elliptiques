from jacobien import *
from random import randint
import hashlib


def hash(mes: str) -> int:
    mes_res = ""
    for i in range(0, len(mes)):
        mes_res += str(ord(mes[i]))
    return int(mes_res)


n = 13

# Préparation des clés, q est la clé publique, s la privée
s = randint(1, n - 1)
g = Point_jacob(6, 4, 1)
q = g.montgomery(s)
print("Génération de clés :\n\tclé privée : ", s, "\n\tclé publique : ", q.to_affine())

# Signature
m = "EC Digital Signature Algorithm"
x = 0
y = 0
k = 1
while x == 0 or y == 0:
    k = randint(1, n - 1)
    p = g.montgomery(k)
    i = p.to_affine().x
    j = p.to_affine().y
    x = i % n
    #y = (inversion_mult(k, n) * (int(hashlib.sha256(m.encode()).hexdigest(), 16)) % n + s * x) % n
    y = (inversion_mult(k, n) * (hash(m)) % n + s * x) % n
print("Signature :\n\t(x,y) = (", x, ",", y, ")\n\tk = ", k)

# Vérification
print("Vérification :\n\tQ appartient à la courbe : ", q.e_courbe_jacob(0, 3))
print("\tnQ = ", q.montgomery(n))
