#include "converter_old.hpp"

#include <gmp.h>
#include <stdio.h>

void GMPQuotientRingTORNS(uint64_t *rns_poly, const mpz_t *quotient_ring_poly,
                          const Parameters &parms) {
    std::vector<uint64_t> modulus = parms.GetModulus();
    mpz_t mp_modulus[parms.GetModulusNum()];
    for (int l = 0; l < modulus.size(); ++l) {
        mpz_init(mp_modulus[l]);
        mpz_import(mp_modulus[l], 1, -1, sizeof(uint64_t), 0, 0,
                   modulus.data() + l);
    }

    // convert to RNS form
    mpz_t mp_remainder;
    mpz_init(mp_remainder);
    for (int i = 0; i < parms.GetPolyDegree(); ++i) {
        for (int q = 0; q < parms.GetModulusNum(); ++q) {
            mpz_mod(mp_remainder, quotient_ring_poly[i], mp_modulus[q]);
            uint64_t *coeff_p = rns_poly + (i + q * parms.GetPolyDegree());
            *coeff_p = 0;
            mpz_export(coeff_p, NULL, -1, sizeof(uint64_t), 0, 0, mp_remainder);
        }
    }

    // clear mpz values
    for (int l = 0; l < modulus.size(); ++l) {
        mpz_clear(mp_modulus[l]);
    }

    mpz_clear(mp_remainder);
}

void QuotientRingToRNS(uint64_t *rns_poly, uint64_t *quotient_ring_poly,
                       const Parameters &parms) {
    // convert numbers to GMP style
    mpz_t mp_coeff[parms.GetPolyDegree()];
    for (int i = 0; i < parms.GetPolyDegree(); ++i) {
        mpz_init(mp_coeff[i]);
        mpz_import(mp_coeff[i], parms.GetCoeffWordLength(), -1,
                   sizeof(quotient_ring_poly[0]), 0, 0,
                   quotient_ring_poly + (i * parms.GetCoeffWordLength()));
    }

    GMPQuotientRingTORNS(rns_poly, mp_coeff, parms);

    // clear mpz values
    for (int i = 0; i < parms.GetPolyDegree(); ++i) {
        mpz_clear(mp_coeff[i]);
    }
}

/***********************************************************/
// boost/integer/extended_euclidean.hpp

// From "The Joy of Factoring", Algorithm 2.7, with a small optimization to
// remove tmps from Wikipedia. Solves mx + ny = gcd(m,n). Returns tuple with
// (gcd(m,n), x, y).

template <class Z>
struct euclidean_result_t {
    Z gcd;
    Z x;
    Z y;
};

template <class Z>
euclidean_result_t<Z> extended_euclidean(Z m, Z n) {
    assert(m >= 1 && n >= 1);

    bool swapped = false;
    if (m < n) {
        swapped = true;
        std::swap(m, n);
    }
    Z u0 = m;
    Z u1 = 1;
    Z u2 = 0;
    Z v0 = n;
    Z v1 = 0;
    Z v2 = 1;
    Z w0;
    Z w1;
    Z w2;
    while (v0 > 0) {
        Z q = u0 / v0;
        w0 = u0 - q * v0;
        w1 = u1 - q * v1;
        w2 = u2 - q * v2;
        u0 = v0;
        u1 = v1;
        u2 = v2;
        v0 = w0;
        v1 = w1;
        v2 = w2;
    }

    euclidean_result_t<Z> result;
    result.gcd = u0;
    if (!swapped) {
        result.x = u1;
        result.y = u2;
    } else {
        result.x = u2;
        result.y = u1;
    }

    return result;
}
/***********************************************************/

inline uint64_t UMulMod64To64(uint64_t op1, uint64_t op2, uint64_t divisor) {
    uint64_t quotient;
    uint64_t remainder;
    asm("mulq %[op2];"
        "divq %[d];"
        : "=a"(quotient), "=d"(remainder)
        : [op2] "g"(op2), [d] "g"(divisor), "a"(op1));
    return remainder;
}

inline uint64_t UMod128To64(uint64_t low, uint64_t high, uint64_t divisor) {
    uint64_t result;
    uint64_t remainder;
    asm("divq %[d];"
        : "=a"(result), "=d"(remainder)
        : [d] "r"(divisor), "a"(low), "d"(high));
    return remainder;
}

uint64_t ModInvFermat(uint64_t x, uint64_t modulus) {
    uint64_t inverse = 1;
    uint64_t power = modulus - 2;
    while (power > 0) {
        if (power & 1 == 1) {
            inverse *= x;
        }

        x *= x;
        power >>= 1;
    }

    return inverse;
}

uint64_t ModInv(uint64_t x, uint64_t modulus) {
    if (x >= (1ULL << 63) || modulus >= (1ULL << 63)) {
        return ModInvFermat(x, modulus);
    }

    auto res = extended_euclidean(static_cast<__int128_t>(x),
                                  static_cast<__int128_t>(modulus));
    long long inverse = res.x % modulus;
    if (inverse < 0) {
        inverse += modulus;
    }
    return static_cast<uint64_t>(inverse);
}

void RNSToGMPInt(mpz_t result, const std::vector<uint64_t> &rns,
                 const std::vector<uint64_t> &modulus) {
    // Garner's algorithm
    // for each step, we solve "modulus_products[k] * t[k] + constants[k] = b[k]
    // (mod. q[k])"
    //      modulus_products[k] = q[0]q[1]...q[k-1] (mod. q[k])
    //      constants[k] = t[0] + t[1]q[0] + ... + t[k-1]q[0]q[1]...q[k-2] (mod.
    //      q[k])
    std::vector<uint64_t> modulus_products(modulus.size(), 1);
    std::vector<uint64_t> constants(modulus.size(), 0);
    std::vector<uint64_t> t(modulus.size());
    for (int k = 0; k < modulus.size(); ++k) {
        uint64_t tmp;
        if (rns[k] < constants[k]) {
            tmp = modulus[k] + rns[k] - constants[k];
        } else {
            tmp = rns[k] - constants[k];
        }

        t[k] = UMulMod64To64(tmp, ModInv(modulus_products[k], modulus[k]),
                             modulus[k]);

        for (int i = k + 1; i < modulus.size(); ++i) {
            constants[i] =
                (constants[i] +
                 UMulMod64To64(t[k], modulus_products[i], modulus[i])) %
                modulus[i];
            modulus_products[i] =
                UMulMod64To64(modulus_products[i], modulus[k], modulus[i]);
        }
    }

    mpz_t mp_modulus[modulus.size()];
    for (int l = 0; l < modulus.size(); ++l) {
        mpz_init(mp_modulus[l]);
        mpz_import(mp_modulus[l], 1, -1, sizeof(uint64_t), 0, 0,
                   modulus.data() + l);
    }

    mpz_t term;
    mpz_init(term);
    for (int k = 0; k < modulus.size(); ++k) {
        mpz_import(term, 1, -1, sizeof(uint64_t), 0, 0, t.data() + k);
        for (int j = 0; j < k; ++j) {
            mpz_mul(term, term, mp_modulus[j]);
        }

        mpz_add(result, result, term);
    }

    // clear mpz values
    for (int l = 0; l < modulus.size(); ++l) {
        mpz_clear(mp_modulus[l]);
    }

    mpz_clear(term);
}

void RNSToGMPQuotientRing(mpz_t *quotient_ring_poly, uint64_t *rns_poly,
                          const Parameters &parms) {
    std::vector<uint64_t> modulus = parms.GetModulus();

    for (int i = 0; i < parms.GetPolyDegree(); ++i) {
        std::vector<uint64_t> rns_coeff(modulus.size());
        for (int l = 0; l < modulus.size(); ++l) {
            rns_coeff[l] = rns_poly[i + l * parms.GetPolyDegree()];
        }

        RNSToGMPInt(quotient_ring_poly[i], rns_coeff, modulus);
    }
}

void RNSToQuotientRing(uint64_t *quotient_ring_poly, uint64_t *rns_poly,
                       const Parameters &parms) {
    mpz_t mp_coeff[parms.GetPolyDegree()];
    for (int i = 0; i < parms.GetPolyDegree(); ++i) {
        mpz_init(mp_coeff[i]);
    }

    RNSToGMPQuotientRing(mp_coeff, rns_poly, parms);

    size_t count = parms.GetCoeffWordLength();
    for (int i = 0; i < parms.GetPolyDegree(); ++i) {
        uint64_t *target =
            quotient_ring_poly + (i * parms.GetCoeffWordLength());
        *target = 0;
        mpz_export(target, &count, -1, sizeof(uint64_t), 0, 0, mp_coeff[i]);
    }
}