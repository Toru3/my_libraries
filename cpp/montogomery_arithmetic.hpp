#include <cstddef>
#include <climits>
#include <cstdint>

namespace detail{
    inline constexpr uint16_t double_width(uint8_t a){ return a; }
    inline constexpr uint32_t double_width(uint16_t a){ return a; }
    inline constexpr uint64_t double_width(uint32_t a){ return a; }
    inline constexpr __uint128_t double_width(uint64_t a){ return a; }
}

template<typename T> using double_width = decltype(detail::double_width(static_cast<T>(0)));

template<typename T> class Montogomery{
    private:
        static constexpr size_t r_bits = sizeof(T)*CHAR_BIT;
        using U = double_width<T>;
        // nx \equiv -1 \pmod r
        inline static T calc_n_prime(const T mod) noexcept{
            T res = 0, t = 0, s = 1;
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
        T n;
        T np;
        T r2;
        T phi;
        inline T reduction(const U t) const noexcept{
            const auto u = static_cast<T>(t*np);
            const auto v = static_cast<U>(u)*static_cast<U>(n);
            const auto w = static_cast<T>((t+v) >> r_bits);
            return w>=n ? w-n : w;
        }
    public:
        inline static bool is_modulo_too_large(const T modulo) noexcept {
            constexpr double r = static_cast<U>(1) << r_bits;
            const double m = modulo;
            return (m-1)*(m-1)+(r-1)*m >= r*r;
        }
        inline Montogomery(const T modulo, const T phi_modulo = 0) noexcept{
            n = modulo;
            np = calc_n_prime(modulo);
            const auto rn = (static_cast<U>(1) << r_bits) % static_cast<U>(n);
            r2 = static_cast<T>(rn * rn % n);
            phi = phi_modulo==0 ? modulo - 1 : phi_modulo;
        }
        inline T convert(const T a) const noexcept{
            return reduction(static_cast<U>(a)*static_cast<U>(r2));
        }
        inline T invert(const T a) const noexcept{
            return reduction(a);
        }
        inline T one() const noexcept{
            return reduction(static_cast<U>(r2));
        }
        inline T add(const T a, const T b) const noexcept{
            const T c = a+b;
            return c>n ? c-n : c;
        }
        inline T sub(const T a, const T b) const noexcept{
            const T c = a-b;
            return a<b ? c+n : c;
        }
        inline T mul(const T a, const T b) const noexcept{
            return reduction(static_cast<U>(a)*static_cast<U>(b));
        }
        inline T sqr(const T a) const noexcept{
            return reduction(static_cast<U>(a)*static_cast<U>(a));
        }
        inline T pow(const T a, T n) const noexcept{
            T x=one(), y=a;
            while(n>0){
                if(n%2==1){
                    x = mul(x, y);
                }
                y = sqr(y);
            }
            return x;
        }
        inline T inv(const T a) const noexcept{
            return pow(a, phi-1);
        }
        inline T div(const T a, const T b) const noexcept{
            return mul(a, inv(b));
        }
};
