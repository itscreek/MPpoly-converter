#include "converter_old.hpp"
#include "gtest/gtest.h"

TEST(ConverterOLDTest, QuotinetRingToRNSTest) {
    std::vector<uint64_t> modulus = {1ULL << 62, (1ULL << 62) + 1};
    Parameters parms(2, modulus, 2);

    std::vector<uint64_t> poly = {0, 1, 0, 1};
    std::vector<uint64_t> expect_rns_poly = {0, 0, (1ULL << 62) - 3,
                                             (1ULL << 62) - 3};

    QuotientRingToRNS(poly.data(), poly.data(), parms);

    ASSERT_EQ(expect_rns_poly, poly);
}

TEST(ConverterOLDTest, RNSToQuotientRingTest) {
    std::vector<uint64_t> modulus = {1ULL << 62, (1ULL << 62) + 1};
    Parameters parms(2, modulus, 2);

    std::vector<uint64_t> poly = {0, 0, (1ULL << 62) - 3, (1ULL << 62) - 3};
    std::vector<uint64_t> expect_quotient_ring_poly = {0, 1, 0, 1};

    RNSToQuotientRing(poly.data(), poly.data(), parms);

    ASSERT_EQ(expect_quotient_ring_poly, poly);
}

TEST(ConverterOLDTest, Test) {
    std::vector<uint64_t> modulus = {102379527169, 204759054337, 409518112769};
    Parameters parms(2, modulus, 2);

    std::vector<uint64_t> poly = {2, 1, 3, 4};
    std::vector<uint64_t> expect_poly = poly;

    std::vector<uint64_t> rns_poly(2 * 3);
    QuotientRingToRNS(rns_poly.data(), poly.data(), parms);
    RNSToQuotientRing(poly.data(), rns_poly.data(), parms);

    ASSERT_EQ(expect_poly, poly);
}