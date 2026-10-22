#pragma once
#include "chem/kumo/local_shell_integrals.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
namespace chem::kumo::sampling {
struct Pair {
    std::size_t a, b;
    std::array<double, 3> center;
    std::uint64_t primitive_pairs;
    std::uint64_t ao_pairs;
};
struct Sample {
    std::uint64_t id;
    unsigned bin;
    long double cost;
    double distance;
};
inline std::uint64_t triangular(std::uint64_t n) { return n*(n+1)/2; }
inline std::pair<std::uint64_t,std::uint64_t> decode(std::uint64_t id) {
    // Canonical pair indices: row j has j+1 entries, i <= j.
    auto j=static_cast<std::uint64_t>((std::sqrt(8.0L*id+1)-1)/2);
    while (triangular(j+1)<=id) ++j;
    while (triangular(j)>id) --j;
    return {id-triangular(j),j};
}
inline unsigned distance_bin(double d) { return d<2?0:d<5?1:d<10?2:3; }
inline std::uint64_t safe_mul(std::uint64_t a, std::uint64_t b) {
    if (b && a>std::numeric_limits<std::uint64_t>::max()/b)
        throw std::overflow_error("primitive work estimate overflow");
    return a*b;
}

struct DrawResult {
    std::vector<Pair> pairs;
    std::array<std::vector<Sample>,4> samples;
    std::uint64_t formal_quartets=0;
    std::size_t unique_draws=0;
};
inline DrawResult sample_quartets(const std::vector<chem::kumo::LocalShellComponents>& shells,
                                 std::size_t per_bin=2500,
                                 std::uint64_t seed=0xC16D1ULL) {
    DrawResult result;
    const std::uint64_t ns=shells.size(), np=triangular(ns);
    if(np>2000000) throw std::runtime_error("shell-pair table too large");
    result.formal_quartets=triangular(np);
    if(!result.formal_quartets) return result;
    result.pairs.reserve(np);
    for(std::size_t b=0;b<shells.size();++b)
        for(std::size_t a=0;a<=b;++a) {
            const auto& sa=shells[a]; const auto& sb=shells[b];
            if(sa.empty()||sb.empty()||sa.front().empty()||sb.front().empty())
                throw std::runtime_error("empty shell");
            std::uint64_t pa=0,pb=0;
            for(const auto& c:sa) pa+=c.size();
            for(const auto& c:sb) pb+=c.size();
            const auto& x=sa.front().front().center;
            const auto& y=sb.front().front().center;
            std::array<double,3> mid{};
            for(int k=0;k<3;++k) mid[k]=0.5*(x[k]+y[k]);
            result.pairs.push_back({a,b,mid,safe_mul(pa,pb),safe_mul(sa.size(),sb.size())});
        }
    std::set<std::uint64_t> seen;
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<std::uint64_t> draw(0,result.formal_quartets-1);
    for(std::uint64_t attempt=0;attempt<2000000;++attempt) {
        bool complete=true;
        for(const auto& v:result.samples) complete &= v.size()>=per_bin;
        if(complete) break;
        const auto id=draw(rng);
        if(!seen.insert(id).second) {
            if(seen.size()==result.formal_quartets) break;
            continue;
        }
        const auto [i,j]=decode(id);
        const auto& a=result.pairs.at(i); const auto& b=result.pairs.at(j);
        double d2=0;
        for(int k=0;k<3;++k) d2+=(a.center[k]-b.center[k])*(a.center[k]-b.center[k]);
        const double d=std::sqrt(d2);
        const unsigned bin=distance_bin(d);
        if(result.samples[bin].size()<per_bin) {
            const auto cost=safe_mul(a.primitive_pairs,b.primitive_pairs);
            result.samples[bin].push_back({id,bin,static_cast<long double>(cost),d});
        }
        if(seen.size()==result.formal_quartets) break;
    }
    result.unique_draws=seen.size();
    return result;
}
}
