#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include "converter.hpp"
#include "hexl/hexl.hpp"

void QuotientRingToRNSTime(const uint64_t poly_degree,
                           const std::vector<uint64_t> &modulus,
                           const int coeff_word_length) {
    constexpr int num_times = 1000;
    uint64_t mp_poly[coeff_word_length * poly_degree];
    uint64_t rns_poly[poly_degree * modulus.size()];
    std::random_device rnd;
    std::mt19937_64 mt(rnd());
    for (int i = 0; i < coeff_word_length * poly_degree; ++i) {
        mp_poly[i] = mt();
    }

    MPPoly::PolynomialFormConverter converter(poly_degree, modulus,
                                              coeff_word_length);

    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < num_times; ++i) {
        converter.QuotientRingToRNS(rns_poly, mp_poly);
    }
    auto end = std::chrono::system_clock::now();
    auto dur = end - start;
    auto msec =
        std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();

    std::cout << num_times << " times convertion from Z_Q to RNS form: " << msec
              << "ms" << std::endl;
    std::cout << "  each convertion time: "
              << msec / static_cast<double>(num_times) << "ms" << std::endl;
}

void RNSToQuotientRingTime(const uint64_t poly_degree,
                           const std::vector<uint64_t> &modulus,
                           const int coeff_word_length) {
    constexpr int num_times = 1000;
    uint64_t mp_poly[coeff_word_length * poly_degree];
    uint64_t rns_poly[poly_degree * modulus.size()];
    std::random_device rnd;
    std::mt19937_64 mt(rnd());
    for (int i = 0; i < poly_degree * modulus.size(); ++i) {
        rns_poly[i] = mt();
    }

    MPPoly::PolynomialFormConverter converter(poly_degree, modulus,
                                              coeff_word_length);

    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < num_times; ++i) {
        converter.RNSToQuotientRing(mp_poly, rns_poly);
    }
    auto end = std::chrono::system_clock::now();
    auto dur = end - start;
    auto msec =
        std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();

    std::cout << num_times << " times convertion from RNS form to Z_Q: " << msec
              << "ms" << std::endl;
    std::cout << "  each convertion time: "
              << msec / static_cast<double>(num_times) << "ms" << std::endl;
}

void PolynomialConvertTime(const uint64_t poly_degree,
                           const std::vector<uint64_t> &modulus,
                           const int coeff_word_length) {
    constexpr int num_times = 1000;
    uint64_t mp_poly[coeff_word_length * poly_degree];
    uint64_t rns_poly[poly_degree * modulus.size()];
    std::random_device rnd;
    std::mt19937_64 mt(rnd());
    for (int i = 0; i < coeff_word_length * poly_degree; ++i) {
        mp_poly[i] = mt();
    }

    MPPoly::PolynomialFormConverter converter(poly_degree, modulus,
                                              coeff_word_length);

    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < num_times; ++i) {
        converter.QuotientRingToRNS(rns_poly, mp_poly);
        converter.RNSToQuotientRing(mp_poly, rns_poly);
    }
    auto end = std::chrono::system_clock::now();
    auto dur = end - start;
    auto msec =
        std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();

    std::cout << num_times << " times convertion: " << msec << "ms"
              << std::endl;
    std::cout << "  each convertion time: "
              << msec / static_cast<double>(num_times) << "ms" << std::endl;
}

void NTTTime(const uint64_t poly_degree, const std::vector<uint64_t> &modulus) {
    constexpr int num_times = 1000;

    uint64_t rns_poly0[poly_degree * modulus.size()];
    uint64_t rns_poly1[poly_degree * modulus.size()];
    std::random_device rnd;
    std::mt19937_64 mt(rnd());
    for (int i = 0; i < poly_degree * modulus.size(); ++i) {
        rns_poly0[i] = mt();
        rns_poly1[i] = mt();
    }

    std::vector<intel::hexl::NTT> ntts(modulus.size());
    int i = 0;
    for (auto &ntt : ntts) {
        ntt = intel::hexl::NTT(poly_degree, modulus[i]);
        ++i;
    }

    auto start = std::chrono::system_clock::now();
    for (int i = 0; i < num_times; ++i) {
        // compute multiplication for each modulus
        for (int l = 0; l < modulus.size(); ++l) {
            int diff = poly_degree * l;
            ntts[l].ComputeForward(rns_poly0 + diff, rns_poly0 + diff, 1, 1);
            ntts[l].ComputeForward(rns_poly1 + diff, rns_poly1 + diff, 1, 1);
            intel::hexl::EltwiseMultMod(rns_poly0 + diff, rns_poly0 + diff,
                                        rns_poly1 + diff, poly_degree,
                                        modulus[l], 1);
            ntts[l].ComputeInverse(rns_poly0 + diff, rns_poly0 + diff, 1, 1);
        }
    }
    auto end = std::chrono::system_clock::now();
    auto dur = end - start;
    auto msec =
        std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();

    std::cout << num_times << " times multiplication: " << msec << "ms"
              << std::endl;
    std::cout << "  each mult time: " << msec / static_cast<double>(num_times)
              << "ms" << std::endl;
}

// measure the time of convert and multiplication
void PolynomialMultiplicationTime(const uint64_t poly_degree,
                                  const std::vector<uint64_t> &modulus,
                                  const int coeff_word_length) {
    uint64_t small_poly[poly_degree];
    uint64_t small_poly_rns[poly_degree * modulus.size()];
    uint64_t large_poly[coeff_word_length * poly_degree];
    uint64_t large_poly_rns[poly_degree * modulus.size()];
    uint64_t res_poly_rns[poly_degree * modulus.size()];

    MPPoly::PolynomialFormConverter converter(poly_degree, modulus,
                                              coeff_word_length);

    std::vector<intel::hexl::NTT> ntts(modulus.size());
    int i = 0;
    for (auto &ntt : ntts) {
        ntt = intel::hexl::NTT(poly_degree, modulus[i]);
        ++i;
    }

    int64_t total_dur = 0;
    std::random_device rnd;
    std::mt19937_64 mt(rnd());
    constexpr int time0 = 10;
    constexpr int time1 = 100;
    for (int k = 0; k < time0; ++k) {
        for (int j = 0; j < poly_degree; ++j) {
            small_poly[j] = mt();
        }

        for (int j = 0; j < coeff_word_length * poly_degree; ++j) {
            large_poly[j] = mt();
        }

        auto start = std::chrono::system_clock::now();
        for (int i = 0; i < time1; ++i) {
            converter.QuotientRingToRNS(large_poly_rns, large_poly);

            for (int j = 0; j < modulus.size(); ++j) {
                int diff = poly_degree * j;
                ntts[j].ComputeForward(small_poly_rns + diff, small_poly, 1, 1);
                ntts[j].ComputeForward(res_poly_rns + diff,
                                       large_poly_rns + diff, 1, 1);
                intel::hexl::EltwiseMultMod(
                    res_poly_rns + diff, small_poly_rns + diff,
                    res_poly_rns + diff, poly_degree, modulus[j], 1);
                ntts[j].ComputeInverse(res_poly_rns + diff, res_poly_rns + diff,
                                       1, 1);
            }

            converter.RNSToQuotientRing(large_poly, large_poly_rns);
        }
        auto end = std::chrono::system_clock::now();
        auto dur = end - start;
        auto msec =
            std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
        total_dur += msec;
    }

    std::cout << time0 << "*" << time1
              << " execution(convert, NTT, reconvert): " << total_dur << "ms"
              << std::endl;
    std::cout << "  each execution: "
              << total_dur / static_cast<double>(time0 * time1) << "ms"
              << std::endl;
}

int main() {
    constexpr uint64_t poly_degree = 8192;
    constexpr int coeff_word_length = 2;
    std::vector<uint64_t> modulus = {51189761537, 102379527169, 204759054337,
                                     409518112769};

    QuotientRingToRNSTime(poly_degree, modulus, coeff_word_length);
    RNSToQuotientRingTime(poly_degree, modulus, coeff_word_length);
    NTTTime(poly_degree, modulus);
}