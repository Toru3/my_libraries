#include <cstddef>
#include <climits>
#include <cstdint>

class Montogomery{
    private:
        static constexpr size_t r_bits = 64;
        // nx \equiv -1 \pmod r
        inline static uint64_t calc_n_prime(const uint64_t mod) noexcept{
            uint64_t res = 0, t = 0, s = 1;
            for(size_t i=0; i<r_bits; i++){
                if((t&1) == 0){
                    t += mod;
                    res |= s;
                }
                t >>= 1;
                s <<= 1;
            }
            return res;
        }
        uint64_t n;
        uint64_t np;
        uint64_t r2;
        uint64_t phi;
        uint64_t reduction(const __uint128_t t) const noexcept{
            const auto u = static_cast<uint64_t>(t*np);
            const auto v = static_cast<__uint128_t>(u)*static_cast<__uint128_t>(n);
            const auto w = static_cast<uint64_t>((t+v) >> r_bits);
            return w>=n ? w-n : w;
        }
    public:
        Montogomery(const uint64_t modulo, const uint64_t phi_modulo = 0) noexcept{
            n = modulo;
            np = calc_n_prime(modulo);
            const auto rn = (static_cast<__uint128_t>(1) << r_bits) % static_cast<__uint128_t>(n);
            r2 = static_cast<uint64_t>(rn * rn % n);
            phi = phi_modulo==0 ? modulo - 1 : phi_modulo;
        }
        uint64_t convert(const uint64_t a) const noexcept{
            return reduction(static_cast<__uint128_t>(a)*static_cast<__uint128_t>(r2));
        }
        uint64_t invert(const uint64_t a) const noexcept{
            return reduction(a);
        }
        uint64_t one() const noexcept{
            return reduction(static_cast<__uint128_t>(r2));
        }
        uint64_t add(const uint64_t a, const uint64_t b) const noexcept{
            const uint64_t c = a+b;
            return c>n ? c-n : c;
        }
        uint64_t sub(const uint64_t a, const uint64_t b) const noexcept{
            const uint64_t c = a-b;
            return a<b ? c+n : c;
        }
        uint64_t mul(const uint64_t a, const uint64_t b) const noexcept{
            return reduction(static_cast<__uint128_t>(a)*static_cast<__uint128_t>(b));
        }
        uint64_t sqr(const uint64_t a) const noexcept{
            return reduction(static_cast<__uint128_t>(a)*static_cast<__uint128_t>(a));
        }
        uint64_t pow(const uint64_t a, uint64_t n) const noexcept{
            uint64_t x=one(), y=a;
            while(n>0){
                if(n%2==1){
                    x = mul(x, y);
                }
                y = sqr(y);
            }
            return x;
        }
        uint64_t inv(const uint64_t a) const noexcept{
            return pow(a, phi-1);
        }
        uint64_t div(const uint64_t a, const uint64_t b) const noexcept{
            return mul(a, inv(b));
        }
};
