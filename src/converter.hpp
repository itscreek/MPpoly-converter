#pragma once

#include <gmp.h>

#include <cassert>
#include <cstdint>
#include <vector>

namespace MPPoly {

uint64_t ModInv(uint64_t x, uint64_t modulus);

class PolynomialFormConverter {
   public:
    PolynomialFormConverter(uint64_t poly_degree,
                            const std::vector<uint64_t> &modulus,
                            int coeff_word_length);

    ~PolynomialFormConverter();

    void QuotientRingToRNS(uint64_t *rns_poly, uint64_t *quotient_ring_poly);

    void RNSToQuotientRing(uint64_t *quotient_ring_poly, uint64_t *rns_poly);

   private:
    uint64_t poly_degree_;

    mpz_t *mp_coeff_;

    int num_modulus_;

    std::vector<uint64_t> modulus_;

    mpz_t *mp_modulus_;

    int coeff_word_length_;

    // modulus_products_[l0][l1] = modulus_[0] * ... * modulus_[l1] mod
    // modulus_[l0] (l1 <= l0-1)
    std::vector<std::vector<uint64_t>> modulus_products_;

    // inverse_of_modulus_products_[l] = (modulus_[0] * ... *
    // modulus_[l-1])^ -1 mod modulus_[l] = modulus_products_[l][l-1] ^ -1
    std::vector<uint64_t> inverse_of_modulus_products_;

    void UInt128QuotientRingToRNS(uint64_t *rns_poly,
                                  uint64_t *quotient_ring_poly);

    void GMPQuotientRingToRNS(uint64_t *output_rns_poly);

    void RNSToGMPIntCoeff(uint64_t *input_rns_poly, int degree);

    void ComputeModulusProducts();

    void ComputeInverseOfModulusProducts();
};

}  // namespace MPPoly