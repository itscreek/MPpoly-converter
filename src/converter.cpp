#include "converter.hpp"

#include <stdlib.h>

#include <algorithm>
#include <utility>

namespace MPPoly {
inline void UAdd128To128(uint64_t *result, uint64_t *op1, uint64_t *op2) {
    if (op1[0] > UINT64_MAX - op2[0]) {
        ++result[1];
    }
    result[0] = op1[0] + op2[0];
    result[1] = op1[1] + op2[1];
}

inline void UMul64To128(uint64_t *result, uint64_t op1, uint64_t op2) {
    asm("mulq %[op2]"
        : "=a"(result[0]), "=d"(result[1])
        : [op2] "g"(op2), "a"(op1));
}

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
    if (high > divisor) {
        high %= divisor;
    }
    uint64_t result;
    uint64_t remainder;
    asm("divq %[d];"
        : "=a"(result), "=d"(remainder)
        : [d] "g"(divisor), "a"(low), "d"(high));
    return remainder;
}

uint64_t ModInv(uint64_t x, uint64_t modulus) {
    if (x > modulus) {
        x %= modulus;
    }
    uint64_t u0 = modulus;
    uint64_t u2 = 0;
    uint64_t v0 = x;
    uint64_t v2 = 1;
    uint64_t w0;
    uint64_t w2;
    while (v0 > 0) {
        uint64_t q = u0 / v0;
        w0 = u0 - q * v0;
        uint64_t s = q * v2;
        if (v2 != 0 && q > UINT64_MAX / v2) {
            s = UMulMod64To64(q, v2, modulus);
        }
        if (u2 >= s) {
            w2 = u2 - s;
        } else {
            uint64_t t = (s - u2) / modulus + 1;
            w2 = modulus * t + u2 - s;
        }
        u0 = v0;
        u2 = v2;
        v0 = w0;
        v2 = w2;
    }

    return u2;
}

PolynomialFormConverter::PolynomialFormConverter(
    uint64_t poly_degree, const std::vector<uint64_t> &modulus,
    int coeff_word_length)
    : poly_degree_(poly_degree),
      num_modulus_(modulus.size()),
      modulus_(modulus),
      coeff_word_length_(coeff_word_length) {
    mp_coeff_ = (mpz_t *)malloc(poly_degree_ * sizeof(mpz_t));
    for (int i = 0; i < poly_degree_; ++i) {
        mpz_init(mp_coeff_[i]);
    }

    mp_modulus_ = (mpz_t *)malloc(num_modulus_ * sizeof(mpz_t));
    for (int i = 0; i < num_modulus_; ++i) {
        mpz_init(mp_modulus_[i]);
        mpz_import(mp_modulus_[i], 1, -1, sizeof(uint64_t), 0, 0,
                   modulus_.data() + i);
    }

    ComputeModulusProducts();
    ComputeInverseOfModulusProducts();
}

PolynomialFormConverter::~PolynomialFormConverter() {
    for (int i = 0; i < poly_degree_; ++i) {
        mpz_clear(mp_coeff_[i]);
    }

    for (int i = 0; i < num_modulus_; ++i) {
        mpz_clear(mp_modulus_[i]);
    }

    free(mp_coeff_);
    free(mp_modulus_);
}

void PolynomialFormConverter::QuotientRingToRNS(uint64_t *rns_poly,
                                                uint64_t *quotient_ring_poly) {
    if (coeff_word_length_ <= 2) {
        UInt128QuotientRingToRNS(rns_poly, quotient_ring_poly);
        return;
    }

    // set mp_coeff_
    for (int i = 0; i < poly_degree_; ++i) {
        mpz_set_ui(mp_coeff_[i], 0);
        mpz_import(mp_coeff_[i], coeff_word_length_, -1, sizeof(uint64_t), 0, 0,
                   quotient_ring_poly + (i * coeff_word_length_));
    }

    GMPQuotientRingToRNS(rns_poly);
}

void PolynomialFormConverter::UInt128QuotientRingToRNS(
    uint64_t *rns_poly, uint64_t *quotient_ring_poly) {
    uint64_t *input_poly = quotient_ring_poly;
    uint64_t copy[coeff_word_length_ * poly_degree_];
    if (rns_poly == quotient_ring_poly) {
        for (int i = 0; i < coeff_word_length_ * poly_degree_; ++i) {
            copy[i] = quotient_ring_poly[i];
        }
        input_poly = copy;
    }

    for (int i = 0; i < poly_degree_; ++i) {
        for (int l = 0; l < num_modulus_; ++l) {
            rns_poly[i + l * poly_degree_] = UMod128To64(
                input_poly[coeff_word_length_ * i],
                input_poly[coeff_word_length_ * i + 1], modulus_[l]);
        }
    }
}

void PolynomialFormConverter::GMPQuotientRingToRNS(uint64_t *output_rns_poly) {
    mpz_t mp_remainder;
    mpz_init(mp_remainder);
    for (int i = 0; i < poly_degree_; ++i) {
        for (int l = 0; l < num_modulus_; ++l) {
            mpz_mod(mp_remainder, mp_coeff_[i], mp_modulus_[l]);
            uint64_t *coeff_p = output_rns_poly + (i + l * poly_degree_);
            *coeff_p = 0;
            mpz_export(coeff_p, NULL, -1, sizeof(uint64_t), 0, 0, mp_remainder);
        }
    }

    mpz_clear(mp_remainder);
}

void PolynomialFormConverter::RNSToQuotientRing(uint64_t *quotient_ring_poly,
                                                uint64_t *rns_poly) {
    if (coeff_word_length_ <= 2) {
        RNSToUInt128QuotientRing(quotient_ring_poly, rns_poly);
        return;
    }

    for (int i = 0; i < poly_degree_; ++i) {
        RNSToGMPIntCoeff(rns_poly, i);
    }

    // export mp_coeff_ to quotient_ring_poly
    size_t count = coeff_word_length_;
    for (int i = 0; i < poly_degree_; ++i) {
        uint64_t *target = quotient_ring_poly + (i * coeff_word_length_);
        if (mpz_sgn(mp_coeff_[i]) == 0) {
            for (int j = 0; j < coeff_word_length_; ++j) {
                target[j] = 0;
            }
        }
        mpz_export(target, &count, -1, sizeof(uint64_t), 0, 0, mp_coeff_[i]);
    }
}

void PolynomialFormConverter::GarnersAlgorithm(uint64_t *t,
                                               uint64_t *input_rns_poly,
                                               const int degree) {
    // Garner's algorithm
    // for each step, we solve
    //  "modulus_products_[k][k - 1] * t[k] + constants[k] = b[k] (mod. q[k])"
    //
    // modulus_products_[k][k - 1] = q[0]q[1]...q[k-1] (mod. q[k])
    // constants[k] = t[0] + t[1]modulus_[0] + ... +
    //  t[k-1]modulus_[0]modulus_[1]...modulus_[k-2] (mod. modulus_[k])
    std::vector<uint64_t> constants(num_modulus_, 0);
    for (int k = 0; k < num_modulus_; ++k) {
        int index_of_rns = degree + k * poly_degree_;
        uint64_t tmp = input_rns_poly[index_of_rns];
        if (tmp < constants[k]) {
            tmp += (modulus_[k] - constants[k]);
        } else {
            tmp -= constants[k];
        }

        t[k] = UMulMod64To64(tmp, inverse_of_modulus_products_[k], modulus_[k]);

        for (int i = k + 1; i < num_modulus_; ++i) {
            if (k == 0) {
                constants[i] = t[k] % modulus_[i];
                continue;
            }
            constants[i] =
                (constants[i] + UMulMod64To64(t[k], modulus_products_[i][k - 1],
                                              modulus_[i])) %
                modulus_[i];
        }
    }
}

void PolynomialFormConverter::RNSToUInt128QuotientRing(
    uint64_t *quotient_ring_poly, uint64_t *rns_poly) {
    // Suppose coeff_word_length_ = 2
    constexpr int word_length = 2;

    uint64_t *input_rns_poly = rns_poly;
    uint64_t copy[poly_degree_ * modulus_.size()];
    if (rns_poly == quotient_ring_poly) {
        for (int i = 0; i < poly_degree_ * modulus_.size(); ++i) {
            copy[i] = rns_poly[i];
        }
        input_rns_poly = copy;
    }

    for (int degree = 0; degree < poly_degree_; ++degree) {
        uint64_t t[num_modulus_] = {0};
        GarnersAlgorithm(t, input_rns_poly, degree);

        quotient_ring_poly[word_length * degree] = 0;
        quotient_ring_poly[word_length * degree + 1] = 0;
        for (int k = 0; k < num_modulus_; ++k) {
            uint64_t term[word_length];
            term[0] = t[k];
            term[1] = 0;
            for (int j = 0; j <= k - 1; ++j) {
                // compute term = term[1]::term[0] * modulus_[j]
                uint64_t term_low_times_mod[word_length] = {0};
                UMul64To128(term_low_times_mod, term[0], modulus_[j]);
                uint64_t term_high_times_mod = term[1] * modulus_[j];
                term[0] = term_low_times_mod[0];
                term[1] = term_low_times_mod[1] + term_high_times_mod;
            }

            UAdd128To128(quotient_ring_poly + (word_length * degree),
                         quotient_ring_poly + (word_length * degree), term);
        }
    }
}

void PolynomialFormConverter::RNSToGMPIntCoeff(uint64_t *input_rns_poly,
                                               const int degree) {
    uint64_t t[num_modulus_] = {0};
    GarnersAlgorithm(t, input_rns_poly, degree);

    mpz_set_ui(mp_coeff_[degree], 0);
    mpz_t term;
    mpz_init(term);
    for (int k = 0; k < num_modulus_; ++k) {
        mpz_import(term, 1, -1, sizeof(uint64_t), 0, 0, t + k);
        for (int j = 0; j <= k - 1; ++j) {
            mpz_mul(term, term, mp_modulus_[j]);
        }

        mpz_add(mp_coeff_[degree], mp_coeff_[degree], term);
    }

    mpz_clear(term);
}

void PolynomialFormConverter::ComputeModulusProducts() {
    modulus_products_ = std::vector<std::vector<uint64_t>>(num_modulus_);

    int modulus_idx = 0;
    for (auto &vec : modulus_products_) {
        if (modulus_idx == 0) {
            ++modulus_idx;
            continue;
        }

        vec = std::vector<uint64_t>(modulus_idx);
        vec[0] = modulus_[0] % modulus_[modulus_idx];
        for (int i = 1; i < modulus_idx; ++i) {
            vec[i] =
                UMulMod64To64(vec[i - 1], modulus_[i], modulus_[modulus_idx]);
        }
        ++modulus_idx;
    }
}

void PolynomialFormConverter::ComputeInverseOfModulusProducts() {
    inverse_of_modulus_products_ = std::vector<uint64_t>(num_modulus_);

    int count = 0;
    for (auto &inv : inverse_of_modulus_products_) {
        if (count == 0) {
            inv = 1;
            ++count;
            continue;
        }

        inv = ModInv(modulus_products_[count][count - 1], modulus_[count]);

        ++count;
    }
}

}  // namespace MPPoly