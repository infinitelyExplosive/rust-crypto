use rug::{Assign, Integer, Rational};
use std::process::Output;

trait RugNumber<'a>: std::ops::Mul<&'a Self>
where
    Self: 'a,
{
}

impl RugNumber<'_> for Integer {}
impl RugNumber<'_> for Rational {}

pub fn find_inverse(e: &Integer, n: &Integer) -> Integer {
    let (_r, _s, t) = extended_euclidian(n, e);

    if t < 0 {
        let t = t + n;
        return t;
    } else {
        return t;
    }
}

pub fn bezout(a: &Integer, b: &Integer) -> (Integer, Integer) {
    if a > b {
        let (_, s, t) = extended_euclidian(a, b);
        return (s, t);
    } else {
        let (_, s, t) = extended_euclidian(b, a);
        return (t, s);
    }
}

pub fn extended_euclidian(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    let a = Integer::from(a);
    let b = Integer::from(b);
    let mut qs = Vec::new();
    let mut rs = Vec::new();
    let mut ss = Vec::new();
    let mut ts = Vec::new();

    rs.push(a);
    rs.push(b);

    ss.push(Integer::from(1));
    ss.push(Integer::from(0));

    ts.push(Integer::from(0));
    ts.push(Integer::from(1));

    // println!("r:{:x} \ts:{:x} \tt:{:x}", rs[0], ss[0], ts[0]);
    // println!("r:{:x} \ts:{:x} \tt:{:x}", rs[1], ss[1], ts[1]);

    let mut step = 1;
    while *rs.last().unwrap() != 0 {
        let new_q = Integer::from(&rs[step - 1] / &rs[step]);
        let new_r = Integer::from(&rs[step - 1] % &rs[step]);
        let new_s = &ss[step - 1] - Integer::from(&new_q * &ss[step]);
        let new_t = &ts[step - 1] - Integer::from(&new_q * &ts[step]);

        // println!("r:{:x} \ts:{:x} \tt:{:x}", new_r, new_s, new_t);
        qs.push(new_q);
        rs.push(new_r);
        ss.push(new_s);
        ts.push(new_t);

        step += 1;
    }

    // println!();
    rs.pop();
    ss.pop();
    ts.pop();
    return (rs.pop().unwrap(), ss.pop().unwrap(), ts.pop().unwrap());
}

pub fn fast_power(x: &Integer, e: &Integer, n: &Integer) -> Integer {
    let mut res = Integer::from(1);

    let mut base = Integer::from(x % n);
    let mut e = Integer::from(e);

    while e > 0 {
        if e.get_bit(0) {
            res *= &base;
            res %= n;
        }
        e >>= 1;
        base.square_mut();
        base %= n;
    }

    return res;
}

pub fn crt<'a>(
    vals: impl Iterator<Item = &'a Integer>,
    mods: impl Iterator<Item = &'a Integer>,
) -> Integer {
    let mut vals = vals.peekable();
    let mut mods = mods.peekable();

    let mut n = Integer::from(mods.next().unwrap());
    let mut a = Integer::from(vals.next().unwrap());

    while vals.peek().is_some() {
        let n2 = mods.next().unwrap();
        let a2 = vals.next().unwrap();

        let (m1, m2) = bezout(&n, n2);

        a = &a * m2 * n2 + a2 * m1 * &n;
        n *= n2;

        a %= &n;
        if a < 0 {
            a += &n;
        }
    }

    return a;
}

pub fn lll(basis_integer: &Vec<Vec<Integer>>) -> Vec<Vec<Integer>> {
    let n = basis_integer[0].len() - 1;
    let delta = Rational::from((3, 4));

    let mut basis: Vec<Vec<Rational>> = basis_integer
        .iter()
        .map(|vec| vec.iter().map(|x| Rational::from((x, 1))).collect())
        .collect();

    // println!("rational basis:");
    print_basis(&basis, 0);

    let mut b_star = gsp(&basis);

    // println!("b*:");
    print_basis(&b_star, 0);

    let mut mu_matrix: Vec<Vec<Rational>> = compute_mus(&basis, &b_star);

    // println!("mu:");
    print_basis(&mu_matrix, 0);
    println!();

    let mut k = 1;
    while k <= n {
        for j in (0..k).rev() {
            // println!("trying k={} j={} mu={}", k, j, mu_matrix[k][j]);
            if mu_matrix[k][j] > 0.5 {
                // let mut to_subtract = Rational::new();
                // print!(" mu*b_j = [");
                for i in 0..basis[0].len() {
                    // print!("{:?}, ", Rational::from(&mu_matrix[k][j]).round());
                    let to_subtract = Rational::from(&mu_matrix[k][j]).round() * &basis[j][i];
                    basis[k][i] -= &to_subtract;
                }
                // println!("]");
                b_star = gsp(&basis);
                mu_matrix = compute_mus(&basis, &b_star);

                // println!(" rational basis:");
                // print_basis(&basis, 1);
                // println!(" b*:");
                // print_basis(&b_star, 1);
                // println!("mu:");
                // print_basis(&mu_matrix, 1);
                // println!();
            }
        }
        if inner_product(&b_star[k], &b_star[k])
            >= (&delta - Rational::from(mu_matrix[k][k - 1].square_ref()))
                * inner_product(&b_star[k - 1], &b_star[k - 1])
        {
            k += 1;
            // println!("increment k to {}\n", k);
        } else {
            basis.swap(k - 1, k);
            b_star = gsp(&basis);
            mu_matrix = compute_mus(&basis, &b_star);
            k = std::cmp::max(k - 1, 1);

            // println!("swap {} {}", k, k-1);
            // println!("rational basis:");
            // print_basis(&basis, 1);
            // println!("b*:");
            // print_basis(&b_star, 1);
            // println!("mu:");
            // print_basis(&mu_matrix, 1);
            // println!();
            // println!("k to {}\n", k);
        }
    }

    let integer_basis: Vec<Vec<Integer>> = basis
        .iter()
        .map(|v| {
            v.iter()
                .map(|elem| Integer::from(Rational::from(elem).round().numer()))
                .collect()
        })
        .collect();
    return integer_basis;
}

fn compute_mus(basis: &Vec<Vec<Rational>>, b_star: &Vec<Vec<Rational>>) -> Vec<Vec<Rational>> {
    return basis
        .iter()
        .map(|v| {
            b_star
                .iter()
                .map(|v_star| {
                    Rational::from(inner_product(v, v_star) / inner_product(v_star, v_star))
                })
                .collect()
        })
        .collect();
}

pub fn gsp(basis: &Vec<Vec<Rational>>) -> Vec<Vec<Rational>> {
    let mut new_basis = Vec::new();

    for (i, vector) in basis.iter().enumerate() {
        // println!("reducing v{} ({:?})", i, vector);
        let mut u_n = basis[i].clone();
        for j in 0..i {
            let sub = proj(&new_basis[j], vector);
            // println!(
            // " subtracting proj_{}{:?} ({:?}) = {:?}",
            // j, new_basis[j], u_n, sub
            // );
            for (k, u_val) in u_n.iter_mut().enumerate() {
                *u_val -= &sub[k];
            }
        }
        // println!("reduced v{} to {:?}", i, u_n);
        new_basis.push(u_n);
    }
    return new_basis;
}

fn proj(u: &Vec<Rational>, v: &Vec<Rational>) -> Vec<Rational> {
    let mut ret = u.clone();
    let uv = inner_product(u, v);
    let uu = inner_product(u, u);

    let mu = uv / uu;
    for val in ret.iter_mut() {
        *val *= &mu;
    }

    return ret;
}

fn inner_product<'a, T, U>(u: &'a [T], v: &'a [U]) -> Rational
where
    &'a T: std::ops::Mul<&'a U>,
    <&'a T as std::ops::Mul<&'a U>>::Output: Into<Rational>,
{
    return u.iter().zip(v).map(|(un, vn)| (un * vn).into()).sum();
}

fn print_basis(basis: &Vec<Vec<Rational>>, indent: i32) {
    for vec in basis {
        for _ in 0..indent {
            print!(" ");
        }
        println!("{:?}", vec);
    }
}
