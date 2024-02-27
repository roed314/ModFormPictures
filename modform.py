# This file contains functions for defining and evaluating modular forms

from sage.misc.lazy_attribute import lazy_attribute
from sage.arith.misc import primes_first_n, next_prime, prime_powers
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.matrix.constructor import matrix
from sage.functions.other import floor
import cmath
from lmfdb import db

def discrete_log(elts, gens, mod):
    # algorithm 2.2, page 16 of https://arxiv.org/abs/0903.2785
    def table_gens(gens, mod):
        T = [1]
        n = len(gens)
        r = [None]*n
        s = [None]*n
        for i in range(n):
            beta = gens[i]
            r[i] = 1
            N = len(T)
            while beta not in T:
                for Tj in T[:N]:
                    T.append((beta*Tj) % mod)
                beta = (beta*gens[i]) % mod
                r[i] += 1
            s[i] = T.index(beta)
        return T, r, s
    T, r, s = table_gens(gens, mod)
    n = len(gens)
    N = [ prod(r[:j]) for j in range(n) ]
    Z = lambda s: [ (floor(s/N[j]) % r[j]) for j in range(n)]
    return [Z(T.index(elt % mod)) for elt in elts]

class ModularForm:
    """
    We store modular forms internally in terms of their q-expansions at one or more cusps.

    INPUT:

    - ``label`` -- the LMFDB label for a classical modular form
    """
    def __init__(self, label, prec=None):
        self.label = label
        self.use_trace_form = (label.count(".") == 3)
        if label.count(".") == 3:
            level, weight, char_orbit, _ = label.split(".")
            self.char_orbit = f"{level}.{char_orbit}"
            self.use_trace_form = True
        else:
            level, weight, char_orbit, _, char, _ = label.split(".")
            self.char_label = f"{level}.{char}"
            self.use_trace_form = False
        self.level = ZZ(level)
        self.levelbound = 1/RR(self.level).sqrt()
        self.weight = ZZ(weight)
        self.prec = prec

        # We store a matrix that translates into the fundamental domain, since this is locally constant
        # so we can avoid recomputing it for every point as we iterate over a plotting grid
        self.to_fund = self.matID = matrix(ZZ, 2, 2, [1,0,0,1])
        self.matS = matrix(ZZ, 2, 2, [0,-1,1,0])

    def in_fundamental_domain(self, z):
        assert self.level == 1
        x, y = z.real(), z.imag()
        return -0.501 < x < 0.501 and x**2 + y**2 > 0.99

    def act(self, A, z):
        return (A[0,0]*z + A[0,1]) / (A[1,0]*z + A[1,1])

    def to_fundamental_domain(self, z):
        """
        INPUT:

        - ``z`` -- a complex number

        OUTPUT:

        - ``w`` -- a complex number in the fundamental domain (or close to it) equivalent to z under the action of SL(2,Z)
        - ``A`` a matrix in SL(2,Z) so that w = z*A
        """
        z0 = z
        assert self.level == 1
        z = self.act(self.to_fund, z)
        while not self.in_fundamental_domain(z):
            x, y = z.real(), z.imag()
            xshift = -floor(x + 0.5)
            self.to_fund = matrix(ZZ,2,2,[1,xshift,0,1]) * self.to_fund
            z += xshift
            if x**2 + y**2 < 1:
                self.to_fund = self.matS * self.to_fund
                z = -1/z
        assert (self.act(self.to_fund, z0) - z).abs() < 0.000001
        return z, self.to_fund

    @lazy_attribute
    def qexp(self):
        if self.use_trace_form:
            return [0] + db.mf_newforms.lookup(self.label, "traces")
        else:
            an_normalized = [0] + db.mf_hecke_cc.lookup(self.label, "an_normalized")
            k = self.weight
            return [an * float(i)**((k-1)/2) for (i,an) in enumerate(an_normalized)]

    # Old version
    @lazy_attribute
    def nf_qexp(self):
        # TODO: Need to set up lmfdb as a dependency appropriately
        data = db.mf_hecke_nf.lucky({'label': modform_label}, ["ap", "hecke_ring_cyclotomic_generator", "hecke_ring_character_values", "field_poly"])
        aps = data["ap"]
        ZZx = PolynomialRing(ZZ, "x")
        K = NumberField(ZZx(data["field_poly"]), "a") # TODO: if hecke_ring_cyclotomic_generator > 0, need to do something else
        primes = primes_first_n(len(aps))
        if self.prec is None:
            an_list_bound = next_prime(primes[-1])
        else:
            an_list_bound = self.prec
            primes = [p for p in primes if p < self.prec]
        good_primes = [p for p in primes if not p.divides(self.level)]
        if not data.get("hecke_ring_character_values"):
            # trivial character
            char_values = dict(zip(good_primes, [1]*len(good_primes)))
        else:
            gens = [elt[0] for elt in hecke_ring_character_values]
            gens_values = [convert_elt_to_field(elt[1]) for elt in hecke_ring_character_values]
            char_values = dict(
                [(
                    p,prod(g**k for g, k in zip(gens_values, elt)))
                 for p, elt in zip(good_primes, discrete_log(good_primes, gens, level))
                 ])
        PS = PowerSeriesRing(QQ, "q")
        an = [0]*an_list_bound
        an[1] = 1
        # We start by computing the values a[n] for n a prime power using Euler factors
        for p, ap in zip(primes, aps):
            if p.divides(level):
                euler_factor = [1, -ap[0]]
            else:
                euler_factor = [1, -ap[0], p**(weight - 1) * char_values[p]]
            k = an_list_bound.exact_log(p) + 1
            foo = (1/R(euler_factor)).padded_list(k)
            for i in range(1, k):
                an[p^i] = foo[i]
        # Now we extend multiplicatively to define a[n] for n not a prime power.
        for pp in prime_powers(len(an)-1):
            for k in range(1, (len(an) - 1)//pp + 1):
                if gcd(k, pp) == 1:
                    an[pp*k] = an[pp]*an[k]
        return an

    def evaluate(self, z, use_modularity=True):
        value = 0
        scale = 1
        N, k = self.level, self.weight
        if use_modularity and N == 1:
            z0 = z
            z, gamma = self.to_fundamental_domain(z)
            a,b,c,d = gamma.list()
            #c, d = -gamma[1,0], gamma[0,0] # bottom row for inverse of gamma
            #c, d = gamma[1,0], gamma[1,1]
            #scale = (c*z0 + d)**k
            scale = (-c*z + a)**k
        elif N.is_prime() and z.imag() < self.levelbound and self.label != "163.3.b.a":
            z = -1/(N*z)
            if self.label == "37.2.a.a":
                scale = (N*z)**k
            else:
                scale = -(N*z)**k
        for idx, coeff in enumerate(self.qexp):
            if coeff != 0:
                value += coeff * cmath.exp(2 * cmath.pi * complex(1j) * idx * z)
        return value * scale
