#include "converter.hpp"
#include "gtest/gtest.h"

namespace MPPoly {

TEST(ModInvTest, Test) {
    constexpr uint64_t x = 2;
    constexpr uint64_t modulus = 61;

    // compute inverse using Fermat's little theorem
    uint64_t inverse = 1;
    uint64_t power = modulus - 2;
    uint64_t p = x;
    while (power > 0) {
        if (power & 1 == 1) {
            inverse *= p;
        }

        p *= p;
        power >>= 1;
    }

    ASSERT_EQ(inverse % modulus, ModInv(x, modulus));
}

TEST(ConverterTest, QuotientRingToRNSTest) {
    std::vector<uint64_t> modulus = {1ULL << 62, (1ULL << 62) + 1};
    PolynomialFormConverter converter(2, modulus, 2);

    std::vector<uint64_t> poly = {0, 1, 0, 1};
    std::vector<uint64_t> expect_rns_poly = {0, 0, (1ULL << 62) - 3,
                                             (1ULL << 62) - 3};

    converter.QuotientRingToRNS(poly.data(), poly.data());

    ASSERT_EQ(expect_rns_poly, poly);
}

TEST(CoverterTest, RNSToQuotientRingTest) {
    std::vector<uint64_t> modulus = {1ULL << 62, (1ULL << 62) + 1};
    PolynomialFormConverter converter(2, modulus, 2);

    std::vector<uint64_t> poly = {0, 0, (1ULL << 62) - 3,
                                             (1ULL << 62) - 3};
    std::vector<uint64_t> expect_quotient_ring_poly = {0, 1, 0, 1};

    converter.RNSToQuotientRing(poly.data(), poly.data());

    ASSERT_EQ(expect_quotient_ring_poly, poly);
}

TEST(ConverterTest, Test){
    std::vector<uint64_t> modulus = {102379527169, 204759054337, 409518112769};
    PolynomialFormConverter converter(2, modulus, 2);

    std::vector<uint64_t> poly = {2, 1, 3, 4};
    std::vector<uint64_t> expect_poly = poly;

    std::vector<uint64_t> rns_poly(2 * 3);
    converter.QuotientRingToRNS(rns_poly.data(), poly.data());
    converter.RNSToQuotientRing(poly.data(), rns_poly.data());

    ASSERT_EQ(expect_poly, poly);
}
}  // namespace MPPoly