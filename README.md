# Courbes_elliptiques
On se place sur une courbe elliptique de type Short Weierstrass curve (y^2 = x^3+ a*x + b) On base nos calculs sur les formules d'addition, de doublement et d'inversion de https://hyperelliptic.org/EFD/
Pour les coordonnées affines :
https://hyperelliptic.org/EFD/g1p/auto-shortw.html

Pour les coordonnées jacobiennes avec a = 0: 

- addition : https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
- doublement : https://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
## EC Diffie-Hellman python
Construction d'un secret commun. Implémenté en python. Lancer courbes_elliptiques.py pour tester.

## EC Diffie-Hellman C
Idem, implémenté en C en utilisant la librairie gmp.
GMP permet de travailler avec des entiers de tailles arbitraires.
On est donc plus limité à 64bits.
https://gmplib.org/manual/Integer-Functions
