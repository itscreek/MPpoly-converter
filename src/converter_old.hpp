#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

class Parameters {
   public:
    Parameters(uint64_t poly_degree, const std::vector<uint64_t> &modulus,
               int coeff_word_length)
        : poly_degree_(poly_degree),
          modulus_(modulus),
          coeff_word_length_(coeff_word_length) {}

    inline uint64_t GetPolyDegree() const { return poly_degree_; }

    inline int GetModulusNum() const { return modulus_.size(); }

    inline uint64_t GetModulus(int index) const {
        assert(index >= 0 && index < modulus_.size());
        return modulus_[index];
    }

    inline std::vector<uint64_t> GetModulus() const {
        return std::vector<uint64_t>(modulus_);
    }

    inline int GetCoeffWordLength() const { return coeff_word_length_; }

   private:
    uint64_t poly_degree_;

    std::vector<uint64_t> modulus_;

    int coeff_word_length_;
};

void QuotientRingToRNS(uint64_t *rns_poly, uint64_t *quotient_ring_poly,
                       const Parameters &parms);

void RNSToQuotientRing(uint64_t *quotinet_ring_poly, uint64_t *rns_poly,
                       const Parameters &parms);