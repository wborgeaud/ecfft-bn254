def to_limbs(x):
    mask = 2^64 - 1
    return [(x>>(64*i)) & mask for i in range(4)]

p = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
F = GF(p)
a, b = 1, 5612291247948481584627780310922020304781354847659642188369727566000581075360
E = EllipticCurve(F, [a,b])
n = E.order().p_primary_part(2)
log_n = n.log(2)

g = E.gens()[0]
G = (g.order()//n) * g
assert G.order() == n

R = E.random_element()
H = [R + i*G for i in range(2^log_n)]
L = [h.xy()[0] for h in H]
S = [L[i] for i in range(0, n, 2)]
S_prime = [L[i] for i in range(1, n, 2)]

s = "\n".join([str(l) for x in L for l in to_limbs(int(x))])
open('bn254_coset', 'w').write(s)

def isogenies(log_n, S, S_prime, E):
    isos = []
    for i in range(log_n, 0, -1):
        n = 1 << i
        nn = n // 2

        for iso in E.isogenies_prime_degree(2):
            psi = iso.x_rational_map()
            if len(set([psi(x) for x in S]))==nn:
                break
        isos.append(psi)
        S = [psi(x) for x in S[:nn]]
        S_prime = [psi(x) for x in S_prime[:nn]]
        E = iso.codomain()

    return isos

isos = isogenies(log_n-1, S, S_prime, E)
print(len(isos))
s = "\n".join([str(l) for psi in isos for coeff in list(psi.numerator())+list(psi.denominator()) for l in to_limbs(int(coeff))])
open('bn254_isogenies', 'w').write(s)


