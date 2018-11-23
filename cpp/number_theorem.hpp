#ifndef NUMBER_THEOREM_HPP
#define NUMBER_THEOREM_HPP
#include <algorithm>
#include <tuple>
namespace NumberTheorem{
    /// binary GCD algorithm (Stein's algorithm)
    template<typename T> inline T gcd(T a, T b) noexcept{
        if(a==0) return b;
        if(b==0) return a;
        a = a<0 ? -a : a;
        b = b<0 ? -b : b;
        const auto sa = a&(-a);
        const auto sb = b&(-b);
        const auto shift = std::min(sa, sb);
        a /= sa;
        b /= sb;
        if(a<b){
            const auto t = a;
            a = b;
            b = t;
        }
        while(a!=b){
            auto t = a-b;
            do{ t >>= 1; }while(t%2==0);
            if(t>b){
                a=t;
            }else{
                a=b;
                b=t;
            }
        }
        return shift*a;
    }
    /// solve ax+by=gcd(a,b)
    template<typename T> inline std::tuple<T, T, T> extended_euclidean_algorithm(const T& a, const T& b) noexcept{
        if(a==0){
            return std::make_tuple(b, 0, 1);
        }
        if(b==0){
            return std::make_tuple(a, 1, 0);
        }
        T r0=a, x0=1, y0=0;
        T r1=b, x1=0, y1=1;
        T r2  , x2  , y2  ;
        while(r1!=0){
            const auto q = r0/r1;
            r2 = r0%r1;   r0 = r1; r1 = r2; 
            x2 = x0-q*x1; x0 = x1; x1 = x2;
            y2 = y0-q*y1; y0 = y1; y1 = y2;
        }
        if(r0<0){
            return std::make_tuple(-r0, -x0, -y0);
        }else{
            return std::make_tuple(r0, x0, y0);
        }
    }
}
#endif
