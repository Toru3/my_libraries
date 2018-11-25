fn gcd_ref_mut<'a, T>(mut a: &'a mut T, mut b: &'a mut T) -> T
where
    T: Clone + num::Zero + std::ops::RemAssign,
{
    while !b.is_zero() {
        *a %= b.clone();
        std::mem::swap(&mut a, &mut b);
    }
    (&*a).clone()
}

pub fn gcd<T>(a: &T, b: &T) -> T
where
    T: Clone + num::Zero + std::ops::RemAssign,
{
    if a.is_zero() {
        b.clone()
    } else if b.is_zero() {
        a.clone()
    } else {
        let mut ia = a.clone();
        let mut ib = b.clone();
        gcd_ref_mut(&mut ia, &mut ib)
    }
}

fn extended_euclidean_algorithm_ref_mut<'a, 'b, 'c, T>(
    mut old: (&'a mut T, &'b mut T, &'c mut T),
    mut now: (&'a mut T, &'b mut T, &'c mut T),
) -> (T, T, T)
where
    T: Clone + num::Zero + num::One + std::ops::RemAssign + std::ops::SubAssign,
    for<'x> &'x T: std::ops::Mul<Output = T> + std::ops::Div<Output = T>,
{
    while !now.0.is_zero() {
        let q = &*old.0 / &*now.0;
        *old.0 %= (&*now.0).clone();
        *old.1 -= &q * &*now.1;
        *old.2 -= &q * &*now.2;
        std::mem::swap(&mut old, &mut now);
    }
    ((&*old.0).clone(), (&*old.1).clone(), (&*old.2).clone())
}

/// solve a*x+b*y=g where g = gcd(a, b)
/// out (g, x, y)
pub fn extended_euclidean_algorithm<T>(a: &T, b: &T) -> (T, T, T)
where
    T: Clone + num::Zero + num::One + std::ops::RemAssign + std::ops::SubAssign,
    for<'x> &'x T: std::ops::Mul<Output = T> + std::ops::Div<Output = T>,
{
    if a.is_zero() {
        (b.clone(), T::zero(), T::one())
    } else if b.is_zero() {
        (a.clone(), T::one(), T::zero())
    } else {
        extended_euclidean_algorithm_ref_mut(
            (&mut a.clone(), &mut T::one(), &mut T::zero()),
            (&mut b.clone(), &mut T::zero(), &mut T::one()),
        )
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn gcd_test() {
        assert_eq!(crate::gcd(&3, &5), 1);
        assert_eq!(crate::gcd(&4, &6), 2);
        assert_eq!(crate::gcd(&34, &55), 1);
        assert_eq!(crate::gcd(&0, &0), 0);
        assert_eq!(crate::gcd(&7, &0), 7);
        assert_eq!(crate::gcd(&0, &42), 42);
    }
    #[test]
    fn eea_test() {
        let f = |a, b| {
            let (g, x, y) = crate::extended_euclidean_algorithm(&a, &b);
            assert_eq!(a * x + b * y, g);
        };
        f(5, 8);
        f(45645, 43276);
        f(21465, 31497);
        f(214654i64, 312497i64);
    }
}
