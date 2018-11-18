#ifndef MILLER_RABIN_HPP
#define MILLER_RABIN_HPP
#include <algorithm>
#include <cmath>
#include <cstdint>

inline uint64_t mul_mod(const uint64_t a, const uint64_t b, const uint64_t m) noexcept{
    const __uint128_t c = static_cast<__uint128_t>(a) * static_cast<__uint128_t>(b);
    return c%m;
}

inline uint64_t sqr_mod(const uint64_t a, const uint64_t m) noexcept{
    const __uint128_t c = static_cast<__uint128_t>(a) * static_cast<__uint128_t>(a);
    return c%m;
}

// a^p \pmod m
inline uint64_t pow_mod(const uint64_t a, uint64_t p, const uint64_t m) noexcept{
    uint64_t x=1, y=a;
    while(p>0){
        x = p&1 ? mul_mod(x, y, m) : x;
        y = sqr_mod(y, m);
        p>>=1;
    }
    return x;
}

inline bool improved_felmat_test(const uint64_t a, const uint64_t n, const uint64_t s) noexcept{
    const uint64_t d = (n-1)/s;
    uint64_t y = pow_mod(a, d, n);
    if(y==1 || y==n-1){
        return true;
    }
    for(uint64_t j=1; j<s; j<<=1){
        y = sqr_mod(y, n);
        if(y==n-1){
            return true;
        }
    }
    return false;
}

inline bool miller_rabin_primality_test(const uint64_t n) noexcept{
    if(n==2){
        return true;
    }
    if(n<2 || n%2==0){
        return false;
    }
    const uint64_t d=n-1;
    const uint64_t s = d&(-d); // max_{s \in 2^k} (d%s==0)
    const uint64_t ua = std::min(n-1, static_cast<uint64_t>(floor(2*log(n)*log(n))));
    for(uint64_t a=2; a<=ua; a++){
        if(!improved_felmat_test(a, n, s)){
            return false;
        }
    }
    return true;
}
#endif
