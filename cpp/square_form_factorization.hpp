#ifndef SQUARE_FORM_FACTORIZATION_HPP
#define SQUARE_FORM_FACTORIZATION_HPP
#include "number_theorem.hpp"
#include "miller-rabin.hpp"
#include <vector>
#include <stdio.h>

namespace SquareFormFactorization{
    /// \lfoor\sqrt{n}\rfloor
    inline uint64_t isqrt(const uint64_t n) noexcept{
        uint64_t sn = floor(sqrt(n));
        if(sn*sn>n){
            sn--;
        }else if((sn+1)*(sn+1)<=n){
            sn++;
        }
        return sn;
    }

    inline bool is_parfect_square(uint64_t n) noexcept{
        constexpr uint64_t mod64 = 0b0000001000000010000000100001001000000010000000110000001000010011;
        if((1<<(n%64) | mod64) == 0){
            return false;
        }
        const uint64_t sn = isqrt(n);
        return sn*sn==n;
    }

    template<typename T> inline bool exist(const std::vector<T>& v, const T& val) noexcept{
        const auto it = std::find(v.cbegin(), v.cend(), val);
        return it != v.cend();
    }

    inline std::pair<uint64_t, uint64_t> foward_step(const uint64_t nn, const uint64_t k, const uint64_t limit, const uint64_t sn) noexcept{
        constexpr size_t max_remainder_size = 64;
        std::vector<uint64_t> q_remain;
        q_remain.reserve(max_remainder_size);
        const uint64_t up = limit*2*k;
        const uint64_t k2 = 2*k;
        uint64_t p0 = sn;
        uint64_t q0 = 1;
        uint64_t q1 = nn-sn*sn;
        for(uint64_t i=0; i<4*limit; i++){
            const auto b = (sn+p0)/q1;
            const auto p1 = b*q1-p0;
            const auto q2 = q0+b*(p0-p1);
            p0 = p1;
            q0 = q1; q1 = q2;
            if(is_parfect_square(q1)){
                const uint64_t sq = isqrt(q1);
                if(!exist(q_remain, sq)){
                    return std::make_pair(sq, p0);
                }
            }else if(q1<up){
                const uint64_t t = q1/NumberTheorem::gcd(k2, q1%(k2));
                if(t<limit){
                    q_remain.push_back(t);
                    if(q_remain.size()>=max_remainder_size){
                        return std::make_pair(1, p0);
                    }
                }
            }
        }
        return std::make_pair(1, p0);
    }

    inline uint64_t backward_step(uint64_t p0, uint64_t q0, uint64_t q1, const uint64_t limit, const uint64_t sn) noexcept{
        for(uint64_t i=0; i<4*limit; i++){
            const auto b = (sn+p0)/q1;
            const auto p1 = b*q1-p0;
            const auto q2 = q0+b*(p0-p1);
            if(p0==p1){
                return q1;
            }
            p0 = p1;
            q0 = q1; q1 = q2;
        }
        return 1;
    }

    inline uint64_t square_form_factorization_aux(const uint64_t n, const uint64_t k, const uint64_t limit) noexcept{
        const uint64_t nn = n*k;
        const uint64_t sn = isqrt(nn);
        const auto [sq, p] = foward_step(nn, k, limit, sn);
        if(sq==1){
            return 1;
        }
        const uint64_t p0 = ((sn-p)/sq)*sq+p;
        const uint64_t q1 = backward_step(p0, sq, (nn-p0*p0)/sq, limit, sn);
        const uint64_t g = NumberTheorem::gcd(q1, 2*k);
        const uint64_t t = q1/g;
        return t;
    }

    inline uint64_t square_form_factorization(const uint64_t n) noexcept{
        if(n==0){ return 0; }
        if(n==1){ return 1; }
        if(MillerRabinPrimalityTest::is_prime(n)){ return 1; }
        if(n%2==0){ return 2; }
        const uint64_t sn = isqrt(n);
        if(sn*sn==n){
            return sn;
        }
        const uint64_t l = floor(2*sqrt(2*sqrt(n)));
        constexpr uint64_t k[] = {1, 3, 5, 7, 11, 3*5, 3*7, 3*11, 5*7, 5*11, 7*11, 3*5*7, 3*5*11, 3*7*11, 5*7*11, 3*5*7*11};
        constexpr size_t len = sizeof(k)/sizeof(k[0]);
        for(size_t i=0; i<len; i++){
            uint64_t f = square_form_factorization_aux(n, k[i], l);
            if(f>1){
                return f;
            }
        }
        return 1;
    }

    template <uint64_t d> inline void small_factor(uint64_t &n, std::vector<uint64_t>& f) noexcept{
        while(n%d==0){
            f.push_back(d);
            n/=d;
        }
    }

    inline std::vector<uint64_t> prime_factorization_aux(uint64_t n) noexcept{
        if(MillerRabinPrimalityTest::is_prime(n)){
            return {n};
        }
        const auto f1 = square_form_factorization(n);
        const auto f2 = n/f1;
        const auto factors1 = prime_factorization_aux(f1);
        const auto factors2 = prime_factorization_aux(f2);
        std::vector<uint64_t> factors;
        factors.reserve(factors1.size()+factors2.size());
        std::merge(factors1.begin(), factors1.end(), factors2.begin(), factors2.end(), std::back_inserter(factors));
        return factors;
    }

    inline std::vector<uint64_t> prime_factorization(uint64_t n) noexcept{
        std::vector<uint64_t> factors;
        if(n==1 || MillerRabinPrimalityTest::is_prime(n)){
            factors.push_back(n);
            return factors;
        }
        small_factor<2>(n, factors);
        small_factor<3>(n, factors);
        small_factor<5>(n, factors);
        small_factor<7>(n, factors);
        small_factor<11>(n, factors);
        small_factor<13>(n, factors);
        if(n==1){
            return factors;
        }
        const auto f = prime_factorization_aux(n);
        std::vector<uint64_t> result;
        std::merge(factors.begin(), factors.end(), f.begin(), f.end(), std::back_inserter(result));
        return result;
    }
}
#endif
