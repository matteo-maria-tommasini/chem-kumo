#include "chem/kumo/local_quartet_sampling.hpp"
#include "chem/kumo/local_pair_bounds.hpp"
#include "chem/kumo/local_shell_pair_envelope.hpp"
#include "chem/kumo/fchk.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    using namespace chem::kumo;
    using Clock = std::chrono::steady_clock;
    try {
        std::size_t target = 1000;
        int arg = 1;
        if (argc > 2 && std::string(argv[1]) == "--target") {
            target = std::stoull(argv[2]);
            arg = 3;
        }
        if (arg == argc || target == 0 || target > 10000 || target % 4 != 0)
            throw std::invalid_argument("usage: profile [--target 1000|5000|10000] file.fchk[.gz] ...");
        std::cout << std::scientific << std::setprecision(10);
        for (; arg < argc; ++arg) {
            const std::string path = argv[arg];
            const BasisSet basis{FormattedCheckpoint(path)};
            const auto shells = local_shell_components(basis);
            const auto t0 = Clock::now();
            const auto sampled = sampling::sample_quartets(shells);
            const auto t1 = Clock::now();
            const auto norms = local_ao_product_norms(basis);
            const auto t2 = Clock::now();
            std::vector<Eigen::Index> offsets(shells.size() + 1, 0);
            for (std::size_t s = 0; s < shells.size(); ++s)
                offsets[s + 1] = offsets[s] + static_cast<Eigen::Index>(shells[s].size());
            std::map<std::size_t, LocalShellPairEnvelope> cache;
            auto envelope = [&](std::size_t p) -> const LocalShellPairEnvelope& {
                auto it = cache.find(p);
                if (it == cache.end()) {
                    const auto& sp = sampled.pairs.at(p);
                    it = cache.emplace(p, prepare_local_shell_pair_envelope(
                        shells[sp.a], shells[sp.b], 0.7)).first;
                }
                return it->second;
            };
            const auto t3 = Clock::now();
            for (unsigned bin = 0; bin < 4; ++bin) {
                const auto count = std::min(target / 4, sampled.samples[bin].size());
                std::size_t improved = 0, cs_skips = 0, hybrid_skips = 0;
                std::size_t invalid = 0, ao_elements = 0, zero_cs = 0;
                long double work = 0, cs_saved = 0, hybrid_saved = 0;
                long double sum_log = 0;
                std::size_t positive_ratios = 0;
                const auto b0 = Clock::now();
                for (std::size_t n = 0; n < count; ++n) {
                    const auto& sample = sampled.samples[bin][n];
                    const auto [pi, pj] = sampling::decode(sample.id);
                    const auto& ab = sampled.pairs.at(pi);
                    const auto& cd = sampled.pairs.at(pj);
                    const auto& eab = envelope(pi);
                    const auto& ecd = envelope(pj);
                    double max_cs = 0.0, max_hybrid = 0.0;
                    bool valid = true;
                    for (std::size_t a = 0; a < eab.na; ++a)
                        for (std::size_t b = 0; b < eab.nb; ++b)
                            for (std::size_t c = 0; c < ecd.na; ++c)
                                for (std::size_t d = 0; d < ecd.nb; ++d) {
                                    const double cs = norms(offsets[ab.a] + static_cast<Eigen::Index>(a),
                                                            offsets[ab.b] + static_cast<Eigen::Index>(b)) *
                                                      norms(offsets[cd.a] + static_cast<Eigen::Index>(c),
                                                            offsets[cd.b] + static_cast<Eigen::Index>(d));
                                    const double spatial = local_shell_pair_envelope_bound(
                                        eab, a, b, ecd, c, d);
                                    if (!std::isfinite(cs) || !std::isfinite(spatial) || cs < 0 || spatial < 0) {
                                        valid = false;
                                        break;
                                    }
                                    const double hybrid = std::min(cs, spatial);
                                    max_cs = std::max(max_cs, cs);
                                    max_hybrid = std::max(max_hybrid, hybrid);
                                    ++ao_elements;
                                }
                    if (!valid) { ++invalid; continue; }
                    improved += (max_hybrid < max_cs);
                    work += sample.cost;
                    constexpr double threshold = 1e-8;
                    if (max_cs < threshold) { ++cs_skips; cs_saved += sample.cost; }
                    if (max_hybrid < threshold) { ++hybrid_skips; hybrid_saved += sample.cost; }
                    if (max_cs > 0 && max_hybrid > 0) {
                        sum_log += std::log(static_cast<long double>(max_hybrid) / max_cs);
                        ++positive_ratios;
                    } else if (max_cs == 0) ++zero_cs;
                }
                const auto b1 = Clock::now();
                std::cout << "system=" << path << " target=" << target << " bin=" << bin
                          << " sampled=" << count << " block_tighter=" << improved
                          << " block_tighter_fraction=" << (count ? double(improved) / count : 0.0)
                          << " ao_elements=" << ao_elements << " invalid=" << invalid
                          << " zero_cs_blocks=" << zero_cs
                          << " geometric_max_bound_ratio_positive="
                          << (positive_ratios ? std::exp(sum_log / positive_ratios) : 1.0L)
                          << " positive_ratios=" << positive_ratios
                          << " cs_skips=" << cs_skips << " hybrid_skips=" << hybrid_skips
                          << " work_proxy=" << work << " cs_saved_work_proxy=" << cs_saved
                          << " hybrid_saved_work_proxy=" << hybrid_saved
                          << " bound_seconds=" << std::chrono::duration<double>(b1 - b0).count() << '\n';
                if (invalid) throw std::runtime_error("invalid full-AO bound");
            }
            const auto t4 = Clock::now();
            std::cout << "system=" << path
                      << " sample_seconds=" << std::chrono::duration<double>(t1 - t0).count()
                      << " norms_seconds=" << std::chrono::duration<double>(t2 - t1).count()
                      << " envelope_cache_entries=" << cache.size()
                      << " profile_seconds=" << std::chrono::duration<double>(t4 - t3).count()
                      << " mode=full_AO_block_bound status=OK\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
