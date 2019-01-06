pub fn gcd<T: num::Zero + for<'x> std::ops::RemAssign<&'x T>>(mut a: T, mut b: T) -> T {
    while !b.is_zero() {
        a %= &b;
        std::mem::swap(&mut a, &mut b);
    }
    a
}

/// solve a*x+b*y=g where g = gcd(a, b)
/// out (g, x, y)
pub fn extended_euclidean_algorithm<T>(a: T, b: T) -> (T, T, T)
where
    for<'x> T: num::Zero + num::One + std::ops::RemAssign<&'x T> + std::ops::SubAssign<&'x T>,
    for<'x> &'x T: std::ops::Mul<Output = T> + std::ops::Div<Output = T>,
{
    let mut old = (a, T::one(), T::zero());
    let mut now = (b, T::zero(), T::one());
    while !now.0.is_zero() {
        let q = &old.0 / &now.0;
        old.0 %= &now.0;
        old.1 -= &(&q * &now.1);
        old.2 -= &(&q * &now.2);
        std::mem::swap(&mut old, &mut now);
    }
    old
}

#[cfg(test)]
mod tests {
    #[test]
    fn gcd_test() {
        assert_eq!(crate::gcd(3, 5), 1);
        assert_eq!(crate::gcd(4, 6), 2);
        assert_eq!(crate::gcd(34, 55), 1);
        assert_eq!(crate::gcd(0, 0), 0);
        assert_eq!(crate::gcd(7, 0), 7);
        assert_eq!(crate::gcd(0, 42), 42);
    }
    #[test]
    fn eea_test() {
        let f = |a, b| {
            let (g, x, y) = crate::extended_euclidean_algorithm(a, b);
            assert_eq!(a * x + b * y, g);
        };
        f(5, 8);
        f(45645, 43276);
        f(21465, 31497);
        f(214654i64, 312497i64);
    }
}
