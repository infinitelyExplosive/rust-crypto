// use ndarray::{s, Array, Array1, Array3, ArrayView3};
use rug::{ops::Pow, Float, Integer, Rational};
use std::{
    fmt::Debug,
    ops::{AddAssign, Mul},
};

pub fn eval_poly(x: &Integer, f: &Vec<Integer>, n: &Integer) -> Integer {
    let mut sum = Integer::from(0);
    let mut x_power = Integer::from(1);
    for term in f {
        sum += Integer::from(&x_power * term);
        x_power *= x;
        if *n > 0 {
            sum %= n;
            x_power %= n;
        }
    }
    return sum;
}

// fn eval_lattice_poly(x: &Integer, f: &Vec<Integer>, const_x: &Integer, n: &Integer) -> Integer {
//     let mut sum = Integer::from(0);
//     let mut x_power = Integer::from(1);
//     let mut const_x_power = Integer::from(1);
//     for term in f {
//         sum += Integer::from(&x_power * term) / &const_x_power;
//         x_power *= x;
//         const_x_power *= const_x;
//     }
//     if *n > 0 {
//         return sum % n;
//     } else {
//         return sum;
//     }
// }

fn eval_rational_lattice_poly(x: &Rational, f: &Vec<Integer>, const_x: &Integer) -> Rational {
    let mut sum = Rational::from(0);
    let mut x_power = Rational::from(1);
    let mut const_x_power = Integer::from(1);
    for term in f {
        sum += Rational::from(&x_power * term) / &const_x_power;
        x_power *= x;
        const_x_power *= const_x;
    }
    return sum;
}

pub fn multiply_poly<'a, T>(f: &'a [T], g: &'a [T]) -> Vec<T>
where
    &'a T: Mul<&'a T>,
    T: From<<&'a T as Mul<&'a T>>::Output>,
    T: From<i32>,
    T: AddAssign<T>,
    T: PartialEq<i32>,
    T: Debug,
{
    let debug = false;
    if debug {
        print!(" mult {:?}*{:?}", f, g);
    }
    let mut prod = Vec::new();
    for _ in 0..(degree(f) + degree(g) + 1) {
        prod.push(T::from(0));
    }

    for i in 0..(degree(f) as usize + 1) {
        for j in 0..(degree(g) as usize + 1) {
            prod[i + j] += T::from(&f[i] * &g[j]);
        }
    }
    if debug {
        println!("\t|\t   result {:?}", prod);
    }
    return prod;
}

pub fn multiply_poly_zn(f: &[Integer], g: &[Integer], n: &Integer) -> Vec<Integer> {
    let debug = false;
    if debug {
        print!(" mult {:?}*{:?}", f, g);
    }
    let mut prod = Vec::new();
    for _ in 0..(degree(f) + degree(g) + 1) {
        prod.push(Integer::from(0));
    }

    for i in 0..(degree(f) as usize + 1) {
        for j in 0..(degree(g) as usize + 1) {
            prod[i + j] += Integer::from(&f[i] * &g[j]);
            if *n > 0 {
                prod[i + j] %= n;
            }
        }
    }
    if debug {
        println!("\t|\t   result {:?}", prod);
    }
    return prod;
}

fn degree<T>(f: &[T]) -> i32
where
    T: PartialEq<i32>,
{
    for (i, val) in f.iter().enumerate().rev() {
        if *val != 0 {
            return i as i32;
        }
    }
    return 0;
}

fn lead<T>(f: &[T]) -> (T, usize)
where
    T: PartialEq<i32>,
    T: From<i32>,
    T: Clone,
{
    for (i, val) in f.iter().enumerate().rev() {
        if *val != 0 {
            return (val.clone(), i);
        }
    }
    return (T::from(0), 0);
}

// pub fn divide_poly(f: &Vec<Rational>, g: &Vec<Rational>) -> (Vec<Rational>, Vec<Rational>) {
//     println!(" div {:?}/{:?}", f, g);
//     assert!(degree(g) >= 0, "divide by 0");
//     let g = &g[0..=(degree(g) as usize)];
//     let mut q: Vec<Rational> = (0..f.len()).map(|_x| Rational::from(0)).collect();
//     let mut r = f.clone();

//     while r.iter().any(|x| *x != 0) && degree(&r) >= degree(g) {
//         let (r_lead, r_power) = lead(&r);
//         let (g_lead, g_power) = lead(&g);
//         let t = r_lead / g_lead;
//         let t_power = r_power - g_power;

//         q[t_power] += &t;

//         print!("  t={}, t_pow={} | ", t, t_power);
//         for (i, val) in g.iter().enumerate() {
//             r[i + t_power] -= Rational::from(val * &t);
//             print!("r[{}]={},  ", i + t_power, r[i + t_power]);
//         }
//         println!("\n  q: {:?}\n  r: {:?}", q, r);
//     }

//     println!(" result q:{:?}, \tr:{:?}", q, r);
//     return (q, r);
// }

pub fn divide_poly_zn(
    f: &Vec<Integer>,
    g: &Vec<Integer>,
    n: &Integer,
) -> (Vec<Integer>, Vec<Integer>) {
    // println!(" div {:?}/{:?}", f, g);
    assert!(degree(g) > 0, "divide by 0");
    let g = &g[0..=(degree(g) as usize)];
    let mut q: Vec<Integer> = (0..f.len()).map(|_x| Integer::from(0)).collect();
    let mut r = f.clone();

    while r.iter().any(|x| *x != 0) && degree(&r) >= degree(g) {
        let (r_lead, r_power) = lead(&r);
        let (g_lead, g_power) = lead(&g);
        let g_inv = find_inverse(&g_lead, n);
        let t = r_lead * g_inv % n;

        let t_power = r_power - g_power;

        q[t_power] += &t;
        q[t_power] %= n;

        // print!("  t={}, t_pow={} | ", t, t_power);
        for (i, val) in g.iter().enumerate() {
            r[i + t_power] -= Integer::from(val * &t);
            r[i + t_power] %= n;
            // print!("r[{}]={},  ", i + t_power, r[i + t_power]);
        }
        // println!("\n  q: {:?}\n  r: {:?}", q, r);
    }

    // println!(" result q:{:?}, \tr:{:?}", q, r);
    return (q, r);
}

pub fn exp_poly(f: &Vec<Integer>, e: &Integer) -> Vec<Integer> {
    if *e == 0 {
        let mut result = Vec::new();
        result.push(Integer::from(1));
        for _ in 0..(f.len() - 1) {
            result.push(Integer::from(0));
        }
        return result;
    }
    let mut result = f.clone();

    for _ in 0..Integer::from(e - 1).to_i32().unwrap() {
        result = multiply_poly(&result, f);
    }
    return result;
}

pub fn determinant(
    matrix: &Vec<Vec<Vec<Integer>>>,
    cols: &Vec<usize>,
    outsize: usize,
    n: &Integer,
) -> Vec<Integer> {
    let mut out = Vec::new();
    for _ in 0..outsize {
        out.push(Integer::from(0));
    }

    if cols.len() == 2 {
        let r0 = matrix.len() - 2;
        let r1 = matrix.len() - 1;
        let positive = multiply_poly_zn(&matrix[r0][cols[0]], &matrix[r1][cols[1]], n);
        let negative = multiply_poly_zn(&matrix[r1][cols[0]], &matrix[r0][cols[1]], n);

        for i in 0..(std::cmp::min(outsize, positive.len())) {
            out[i] += &positive[i];
        }
        for i in 0..(std::cmp::min(outsize, negative.len())) {
            out[i] -= &negative[i];
            if *n > 0 {
                out[i] %= n;
            }
        }
        return out;
    }

    let row = matrix.len() - cols.len();
    for (i, col) in cols.iter().enumerate() {
        let sub_cols = cols.iter().filter(|x| **x != *col).cloned().collect();
        let sub_determinant = determinant(matrix, &sub_cols, outsize, n);
        let product = multiply_poly_zn(&matrix[row][*col], &sub_determinant, n);
        for j in 0..(std::cmp::min(outsize, product.len())) {
            if i % 2 == 0 {
                out[j] += &product[j];
            } else {
                out[j] -= &product[j];
            }
            if *n > 0 {
                out[j] %= n;
            }
        }
    }

    return out;
}

pub fn resultant(f: &Vec<Vec<Integer>>, g: &Vec<Vec<Integer>>, n: &Integer) -> Vec<Integer> {
    let debug = false;

    let f_degree = f.len() - 1;
    let g_degree = g.len() - 1;
    let mut y_vec = Vec::new();
    y_vec.push(Integer::from(0));

    let mut s_matrix = Vec::new();

    for i in 0..g_degree {
        let mut col = Vec::new();
        for _ in 0..i {
            let val = y_vec.clone();
            col.push(val);
        }
        col.append(&mut f.clone());
        for _ in 0..(g_degree - i - 1) {
            let val = y_vec.clone();
            col.push(val);
        }
        s_matrix.push(col);
    }
    for i in 0..f_degree {
        let mut col = Vec::new();
        for _ in 0..i {
            let val = y_vec.clone();
            col.push(val);
        }
        col.append(&mut g.clone());
        for _ in 0..(f_degree - i - 1) {
            let val = y_vec.clone();
            col.push(val);
        }
        s_matrix.push(col);
    }

    // for col in &s_matrix {
    //     for val in col {
    //         print!("{:?}, ", val);
    //     }
    //     println!();
    // }

    // let mut max_len = 20;
    // for row in &s_matrix {
    //     for col in row {
    //         let string = format!("{:?}", col);
    //         if string.len() > max_len {
    //             max_len = string.len();
    //         }
    //     }
    // }
    if debug {
        for i in 0..s_matrix.len() {
            for j in 0..s_matrix.len() {
                let string = format!("{:?}", s_matrix[j][i]);
                // print!("{:?},  ", s_matrix[j][i]);
                // let pad = if string.len() >= max_len {0} else {max_len - string.len()};
                print!("{}, ", string);
                // print!("{:<1$}", "", 4 * (4 - s_matrix[j][i].len()));
            }
            println!();
        }
    }
    let determinant = determinant(&s_matrix, &(0..s_matrix.len()).collect(), 10, n);
    if debug {
        println!("result: {:?}", determinant);
    }
    return determinant;
}

pub fn gcd(a: &Integer, b: &Integer) -> Integer {
    let (r, _s, _t) = extended_euclidean(a, b);
    // println!("gcd {}", r);
    return r;
}

pub fn find_inverse(e: &Integer, n: &Integer) -> Integer {
    let e = ((e.clone() % n) + n) % n;
    let (_r, _s, t) = extended_euclidean(n, &e);

    if t < 0 {
        let t = t + n;
        return t;
    } else {
        return t;
    }
}

pub fn bezout(a: &Integer, b: &Integer) -> (Integer, Integer) {
    if a > b {
        let (_, s, t) = extended_euclidean(a, b);
        return (s, t);
    } else {
        let (_, s, t) = extended_euclidean(b, a);
        return (t, s);
    }
}

pub fn extended_euclidean(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    let a_negative = *a < 0;
    let b_negative = *b < 0;
    let a = Integer::from(a).abs();
    let b = Integer::from(b).abs();
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
    let result_s = if a_negative {-ss.pop().unwrap()} else {ss.pop().unwrap()};
    let result_t = if b_negative {-ts.pop().unwrap()} else {ts.pop().unwrap()};
    return (rs.pop().unwrap(), result_s, result_t);
}

pub fn poly_extended_euclidean_zn(
    a: &Vec<Integer>,
    b: &Vec<Integer>,
    n: &Integer,
) -> (Vec<Integer>, Vec<Integer>, Vec<Integer>) {
    let debug = false;

    let a = a.clone();
    let b = b.clone();
    let mut qs = Vec::new();
    let mut rs = Vec::new();
    let mut ss = Vec::new();
    let mut ts = Vec::new();

    rs.push(a);
    rs.push(b);

    let mut identity_func = Vec::new();
    identity_func.push(Integer::from(1));
    ss.push(identity_func);
    let mut zero_func = Vec::new();
    zero_func.push(Integer::from(0));
    ss.push(zero_func);

    let mut zero_func = Vec::new();
    zero_func.push(Integer::from(0));
    ts.push(zero_func);
    let mut identity_func = Vec::new();
    identity_func.push(Integer::from(1));
    ts.push(identity_func);

    if debug {
        println!("r:{:?} \ts:{:?} \tt:{:?}", rs[0], ss[0], ts[0]);
        println!("r:{:?} \ts:{:?} \tt:{:?}", rs[1], ss[1], ts[1]);
    }

    let mut step = 1;
    while rs.last().unwrap().iter().any(|x| *x != 0) {
        let (new_q, new_r) = divide_poly_zn(&rs[step - 1], &rs[step], n);

        let q_times_s = multiply_poly_zn(&new_q, &ss[step], n);
        let mut new_s = ss[step - 1].clone();
        for _ in 0..(q_times_s.len() - new_s.len()) {
            new_s.push(Integer::from(0));
        }
        for (elem_s, elem_q) in new_s.iter_mut().zip(q_times_s) {
            *elem_s -= elem_q;
            *elem_s %= n;
        }

        let q_times_t = multiply_poly_zn(&new_q, &ts[step], &n);
        let mut new_t = ts[step - 1].clone();
        for _ in 0..(q_times_t.len() - new_t.len()) {
            new_t.push(Integer::from(0));
        }
        for (elem_t, elem_q) in new_t.iter_mut().zip(q_times_t) {
            *elem_t -= elem_q;
            *elem_t %= n;
        }

        if debug {
            println!("r:{:?} \ts:{:?} \tt:{:?}", new_r, new_s, new_t);
        }
        qs.push(new_q);
        rs.push(new_r);
        ss.push(new_s);
        ts.push(new_t);

        step += 1;
    }

    if debug {
        println!();
    }
    rs.pop();
    ss.pop();
    ts.pop();
    return (rs.pop().unwrap(), ss.pop().unwrap(), ts.pop().unwrap());
}

// pub fn poly_extended_euclidean(
//     a: &Vec<Integer>,
//     b: &Vec<Integer>,
// ) -> (Vec<Rational>, Vec<Rational>, Vec<Rational>) {
//     let debug = true;

//     let a: Vec<Rational> = a.iter().map(|x| Rational::from(x)).collect();
//     let b = b.iter().map(|x| Rational::from(x)).collect();
//     let mut qs = Vec::new();
//     let mut rs = Vec::new();
//     let mut ss = Vec::new();
//     let mut ts = Vec::new();

//     rs.push(a);
//     rs.push(b);

//     let mut identity_func = Vec::new();
//     identity_func.push(Rational::from(1));
//     ss.push(identity_func);
//     let mut zero_func = Vec::new();
//     zero_func.push(Rational::from(0));
//     ss.push(zero_func);

//     let mut zero_func = Vec::new();
//     zero_func.push(Rational::from(0));
//     ts.push(zero_func);
//     let mut identity_func = Vec::new();
//     identity_func.push(Rational::from(1));
//     ts.push(identity_func);

//     if debug {
//         println!("r:{:?} \ts:{:?} \tt:{:?}", rs[0], ss[0], ts[0]);
//         println!("r:{:?} \ts:{:?} \tt:{:?}", rs[1], ss[1], ts[1]);
//     }

//     let mut step = 1;
//     while rs.last().unwrap().iter().any(|x| *x != 0) {
//         let (new_q, new_r) = divide_poly(&rs[step - 1], &rs[step]);

//         let q_times_s = multiply_poly(&new_q, &ss[step]);
//         let mut new_s = ss[step - 1].clone();
//         for _ in 0..(q_times_s.len() - new_s.len()) {
//             new_s.push(Rational::from(0));
//         }
//         for (elem_s, elem_q) in new_s.iter_mut().zip(q_times_s) {
//             *elem_s -= elem_q;
//         }

//         let q_times_t = multiply_poly(&new_q, &ts[step]);
//         let mut new_t = ts[step - 1].clone();
//         for _ in 0..(q_times_t.len() - new_t.len()) {
//             new_t.push(Rational::from(0));
//         }
//         for (elem_t, elem_q) in new_t.iter_mut().zip(q_times_t) {
//             *elem_t -= elem_q;
//         }

//         if debug {
//             println!("r:{:?} \ts:{:?} \tt:{:?}", new_r, new_s, new_t);
//         }
//         qs.push(new_q);
//         rs.push(new_r);
//         ss.push(new_s);
//         ts.push(new_t);

//         step += 1;
//     }

//     if debug {
//         println!();
//     }
//     rs.pop();
//     ss.pop();
//     ts.pop();
//     return (rs.pop().unwrap(), ss.pop().unwrap(), ts.pop().unwrap());
// }

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

fn lift(x0: &Integer, a: &Integer, j: u32) -> Integer {
    let y0 = ((Integer::from(a) - Integer::from(x0).pow(2)) / Integer::from(2).pow(j)) % 2;
    let x1 = Integer::from(x0) + Integer::from(2).pow(j - 1) * y0;
    return x1;
}

/*
 * Implemented cases so far:
 * x^2 + c = 0 (mod 2^n)
 * 2áx^2 + 2b́x + 2ć = 0 (mod 2^n)
 * ax^2 + 2b́x + c = 0 (mod 2^n) where p2(r)!=0 and q=1 (m 8)
*/
pub fn solve_quadratic(a: &Integer, b: &Integer, c: &Integer, n: u32) -> Vec<Integer> {
    let debug = true;
    if *a == 0 && *b == 0 {
        // if debug {
        //     println!("a: {} b: {} c: {} n: {}", a, b, c, n);
        // }
        let c = Integer::from(-c);
        let mut x = Integer::from(1);
        let mut n0 = 3;
        while n0 < n {
            x = lift(&x, &c, n0);
            n0 += 1;
        }
        let two_n = Integer::from(2).pow(n);
        let mut results = Vec::new();
        results.push(Integer::from(&x));
        results.push(-Integer::from(&x) + &two_n);
        results.push((Integer::from(&x) + Integer::from(2).pow(n - 1)) % &two_n);
        results.push(-(Integer::from(&x) + Integer::from(2).pow(n - 1)) + &two_n);

        return results;
    }

    let p2a = p2(a, n);
    let p2b = p2(b, n);
    let p2c = p2(c, n);
    let t = [p2a, p2b, p2c].iter().map(|x| x.clone()).min().unwrap();

    if debug {
        println!("a = {}, b = {}, c = {}, n = {}", a, b, c, n);
        println!(
            "t = {}, p2(a) = {}, p2(b) = {} p2(c) = {}",
            t,
            p2a - t,
            p2b - t,
            p2c - t
        );
    }

    if t > 0 {
        let small_n = n - t;
        let two_t = Integer::from(2).pow(t);
        let new_a = Integer::from(a / &two_t);
        let new_b = Integer::from(b / &two_t);
        let new_c = Integer::from(c / &two_t);

        let results = solve_quadratic(&new_a, &new_b, &new_c, small_n);

        let two_n = Integer::from(2).pow(n);
        let two_nt = Integer::from(2).pow(n - t);
        let mut solutions = Vec::new();
        let mut rt = Integer::from(0);

        while rt < two_n {
            for result in &results {
                solutions.push(Integer::from(result + &rt));
            }
            rt += &two_nt;
        }
        return solutions;
    }

    let b_pr = Integer::from(b / 2);
    let b_pr_squared = b_pr.clone().pow(2);
    let a_squared = Integer::from(a).pow(2);
    let two_n = Integer::from(2).pow(n);
    let a_inv = find_inverse(&a, &two_n);
    let a_squared_inv = find_inverse(&a_squared, &two_n);

    let ab = Integer::from(&b_pr_squared * &a_squared_inv);
    let ac = Integer::from(&a_inv * c);
    let s;
    // if p2c < 2 * (p2b - 1) {
    //     // s = ((Integer::from(&ab - &ac) % 8) + 8) % 8; // figure out why this is wrong?
    //     s = ((Integer::from(&ab - &ac) % &two_n) + &two_n) % &two_n;
    // } else {
    s = ((Integer::from(&ab - &ac) % &two_n) + &two_n) % &two_n;
    // }
    let r = p2(&s, n);
    let p2r = p2(&Integer::from(r), n);
    let q = o2(&s, n);
    let q_mod_8 = Integer::from(&q % 8);
    if debug {
        println!(
            " s = {}, r = {}, p2(r) = {}, q = {} = {} (mod 8)",
            s, r, p2r, q, q_mod_8
        );
    }
    if p2r > 0 && q_mod_8 == 1 {
        if debug {
            println!(">>>>>>>>>>>>>>>>>>>>>>>>>>>>");
            println!("(x + {} {})^2 = {}", a_inv, b_pr, s);
        }
        // let mut results = lemma_4(&s, n);
        let new_a = Integer::from(0);
        let new_b = Integer::from(0);
        let new_c = -Integer::from(s);
        let mut results = solve_quadratic(&new_a, &new_b, &new_c, n);

        let offset = Integer::from(&a_inv * &b_pr);
        for result in results.iter_mut() {
            *result -= &offset;
            *result %= &two_n;
            *result += &two_n;
            *result %= &two_n;
        }

        return results;
    }
    return Vec::new();
}

fn p2(a: &Integer, n: u32) -> u32 {
    if *a == 0 {
        return n;
    }
    let mut i = 0;
    let mut a = a.clone();
    while a.is_even() && a != 0 {
        a >>= 1;
        i += 1;
    }
    return i;
}

fn o2(a: &Integer, n: u32) -> Integer {
    let p2a = p2(a, n);
    return a.clone() / Integer::from(2).pow(p2a);
}

pub fn coppersmith(f: &Vec<Integer>, n: &Integer, m: u32, epsilon_denom: u32) -> Integer {
    let debug = false;

    let d = f.len() as u32 - 1;
    // let m = n.significant_bits() / d;
    let w = (d * (m + 1)) as usize;

    let x_exponent = epsilon_denom - d;
    let x_root = d * epsilon_denom;
    let mut x_powers = Vec::new();
    x_powers.push(Integer::from(1));
    x_powers.push(n.clone().pow(x_exponent).root(x_root));
    for _ in 2..(d * (m + 1)) {
        x_powers.push(Integer::from(x_powers.last().unwrap() * &x_powers[1]));
    }

    // if debug {
    println!(
        "d = {}  N = {}  1/e = {}  X = {} ",
        d, n, epsilon_denom, x_powers[1]
    );
    println!("m = {}  w = {}", m, w);
    println!("{}", x_powers.len());
    // }

    let mut basis = Vec::new();
    let mut det = Integer::from(1);

    if debug {
        println!("starting with basis");
    }
    for v in 0..(m + 1) {
        for u in 0..d {
            let mut g_uv = exp_poly(f, &Integer::from(v));

            for _ in 0..u {
                g_uv.insert(0, Integer::from(0));
            }
            if debug {
                print!("  g_{},{} (len {}) [", u, v, g_uv.len());
            }
            for (i, coef) in g_uv.iter_mut().enumerate() {
                if debug {
                    if *coef == 0 {
                        print!("0, ")
                    } else if *coef == 1 {
                        print!("m{}x{}, ", (m - v), i);
                    } else {
                        print!("{}m{}x{}, ", coef, (m - v), i);
                    }
                }
                *coef *= Integer::from(n.pow(m - v));
                *coef *= &x_powers[i];
            }
            det *= &g_uv[(u + v * d) as usize];
            for _ in g_uv.len()..w {
                g_uv.push(Integer::from(0));
                if debug {
                    print!("0, ");
                }
            }
            if debug {
                println!("]");
            }
            basis.push(g_uv);
        }
    }
    if debug {
        println!();
    }
    let exp_two_w4 = Float::with_val(1024, 2).pow(w as f64 / 4.0);
    let left_bound = det.root(w as u32) * exp_two_w4;

    let right_bound = Integer::from(n.pow(m)) / Float::with_val(1024, w).sqrt();
    println!(
        "condition:{} frac {:.3}",
        left_bound < right_bound,
        left_bound / right_bound
    );

    let (reduced_basis, _min_idx) = lll(&basis);

    if debug {
        for v in 0..(m + 1) {
            for u in 0..d {
                print!(" v_{},{}  [", u, v);
                for coef in &reduced_basis[(u + d * v) as usize] {
                    print!("{:}, ", coef);
                }
                let norm = inner_product(
                    &reduced_basis[(u + d * v) as usize],
                    &reduced_basis[(u + d * v) as usize],
                )
                .to_f64();
                println!("] norm {:.2e}", norm);
            }
        }
    }

    for reduced_poly in &reduced_basis[0..2] {
        let search_range = 10;
        let guesses = approximate_zero(reduced_poly, &x_powers[1]);
        for guess_x in guesses {
            // if debug {
            //     println!("guess x={}", guess_x);
            // }

            for i in -search_range..=search_range {
                let x = Integer::from(&guess_x + i);
                let f_of_x = eval_poly(&x, f, n);
                if f_of_x == 0 {
                    return x;
                }
            }
        }
    }
    return Integer::from(-1);
}

pub fn lll(basis_integer: &Vec<Vec<Integer>>) -> (Vec<Vec<Integer>>, usize) {
    let debug = false;

    let n = basis_integer[0].len() - 1;
    let delta = Rational::from((5, 6));

    let mut basis: Vec<Vec<Rational>> = basis_integer
        .iter()
        .map(|vec| vec.iter().map(|x| Rational::from(x)).collect())
        .collect();

    let min_norm = basis.iter().map(|v| l2_norm_squared(v)).min().unwrap();

    if debug {
        println!("rational basis:");
        print_basis(&basis, 0);
    }

    let (mut b_star, mut mu_matrix) = gsp(&basis);

    if debug {
        println!("b*:");
        print_basis(&b_star, 0);
        println!("mu:");
        print_basis(&mu_matrix, 0);
        println!();
    }

    let mut k = 1;
    // println!("{:.^1$}", "goal", n); // performance meter for large computations
    while k <= n {
        // print!("\x1b[2K\r");
        // print!("{:.<1$}", "", k);
        // io::stdout().flush().unwrap();
        for j in (0..k).rev() {
            if debug {
                println!("trying k={} j={} mu={}", k, j, mu_matrix[k][j]);
            }
            if mu_matrix[k][j] >= 0.5 || mu_matrix[k][j] <= -0.5 {
                if debug {
                    print!(" mu*b_j = [");
                }
                for i in 0..basis[0].len() {
                    if debug {
                        print!("{:?}, ", Rational::from(&mu_matrix[k][j]).round());
                    }
                    let to_subtract = Rational::from(&mu_matrix[k][j]).round() * &basis[j][i];
                    basis[k][i] -= &to_subtract;
                }
                if debug {
                    println!("]");
                }
                gsp_efficient(&basis, &mut b_star, &mut mu_matrix, k);

                if debug {
                    println!(" rational basis:");
                    print_basis(&basis, 1);
                    println!(" b*:");
                    print_basis(&b_star, 1);
                    println!("mu:");
                    print_basis(&mu_matrix, 1);
                    println!();
                }
            }
        }

        if inner_product(&b_star[k], &b_star[k])
            >= (&delta - Rational::from(mu_matrix[k][k - 1].square_ref()))
                * inner_product(&b_star[k - 1], &b_star[k - 1])
        {
            k += 1;
            if debug {
                println!("increment k to {}\n", k);
            }
        } else {
            basis.swap(k - 1, k);
            gsp_efficient(&basis, &mut b_star, &mut mu_matrix, k - 1);

            k = std::cmp::max(k - 1, 1);

            if debug {
                println!("swap {} {}", k, k - 1);
                println!("rational basis:");
                print_basis(&basis, 1);
                println!("b*:");
                print_basis(&b_star, 1);
                println!("mu:");
                print_basis(&mu_matrix, 1);
                println!();
                println!("k to {}\n", k);
            }
        }
    }
    if debug {
        println!();
    }

    // check ∀1≤i≤n, j<i. |μ_i,j|≤1/2
    for i in 0..basis.len() {
        for j in 0..i {
            if mu_matrix[i][j] > (1, 2) {
                println!("mu_{},{} was {:.3}", i, j, mu_matrix[i][j].to_f32());
            }
            assert!(mu_matrix[i][j] <= (1, 2));
        }
    }

    // check ∀1≤i<n. δ‖ ̃b_i‖2 ≤ ‖μ_i+1,i  ̃b_i +  ̃b_i+1‖^2
    for i in 0..(basis.len() - 1) {
        let lhs = &delta * l2_norm_squared(&b_star[i]);
        let mu_bi: Vec<Rational> = b_star[i]
            .iter()
            .map(|elem| Rational::from(&mu_matrix[i + 1][i] * elem))
            .collect();
        let mu_bi_bip1: Vec<Rational> = mu_bi
            .iter()
            .zip(&b_star[i + 1])
            .map(|(lhs, rhs)| Rational::from(lhs + rhs))
            .collect();
        let rhs = l2_norm_squared(&mu_bi_bip1);

        if lhs > rhs {
            println!("failed on b_{0} > ub_{0} + b_{0}+1", i)
        }
        assert!(lhs <= rhs);
    }

    let basis_output: Vec<Vec<Integer>> = basis
        .iter()
        .map(|v| {
            v.iter()
                .map(|elem| Rational::from(elem).round().into_numer_denom().0)
                .collect()
        })
        .collect();

    let mut v_norms: Vec<Rational> = Vec::new();
    if debug {
        println!("min norm^2 is      {}", min_norm);
    }
    for (i, v) in basis.iter().enumerate() {
        let v_norm = l2_norm_squared(v);
        if debug {
            println!("v{:2} norm^2 is      {}", i, v_norm);
        }
        v_norms.push(v_norm);
    }
    let min_norm = v_norms.iter().min().unwrap();
    let min_idx = v_norms.iter().position(|x| x == min_norm).unwrap();
    if debug {
        println!("min reduced norm^2 {} (idx {})", min_norm, min_idx);
        println!();
    }
    return (basis_output, min_idx);
}

pub fn gsp(basis: &Vec<Vec<Rational>>) -> (Vec<Vec<Rational>>, Vec<Vec<Rational>>) {
    let debug = false;

    let mut max_denom = Integer::from(0);

    let mut new_basis = Vec::new();
    let mut mus = Vec::new();
    for (i, vector) in basis.iter().enumerate() {
        let mut mu_row = Vec::new();
        if debug {
            println!("reducing v{} ({:?})", i, vector);
        }
        let mut u_n = vector.clone();
        for j in 0..i {
            let (sub, mu) = proj(&new_basis[j], vector);
            if debug {
                println!(
                    " subtracting proj_{}{:?} ({:?}) = {:?}",
                    j, new_basis[j], vector, sub
                );
            }
            for (k, u_val) in u_n.iter_mut().enumerate() {
                if *sub[k].denom() > max_denom {
                    max_denom = Integer::from(sub[k].denom());
                }
                *u_val -= &sub[k];
                // limit_precision_mut(u_val, 64);
            }
            mu_row.push(mu);
        }
        if debug {
            println!("reduced v{} to {:?}", i, u_n);
        }
        new_basis.push(u_n);
        mus.push(mu_row);
    }
    if debug {
        println!("max denom {}", max_denom);
    }
    return (new_basis, mus);
}

pub fn gsp_efficient(
    basis: &Vec<Vec<Rational>>,
    b_star: &mut Vec<Vec<Rational>>,
    mus: &mut Vec<Vec<Rational>>,
    mut updated_row: usize,
) {
    let debug = false;

    while updated_row < basis.len() {
        let vector = &basis[updated_row];
        let mu_row = &mut mus[updated_row];
        mu_row.clear();

        if debug {
            println!("reducing v{} ({:?})", updated_row, vector);
        }

        let mut u_n = vector.clone();
        for j in 0..updated_row {
            let (sub, mu) = proj(&b_star[j], vector);
            if debug {
                println!(
                    " subtracting proj_{}{:?} ({:?}) = {:?}",
                    j, basis[j], vector, sub
                );
            }
            for (k, u_val) in u_n.iter_mut().enumerate() {
                *u_val -= &sub[k];
                // limit_precision_mut(u_val, 1024);
            }
            mu_row.push(mu);
        }
        if debug {
            println!("reduced v{} to {:?}", updated_row, u_n);
        }
        b_star[updated_row] = u_n;
        updated_row += 1;
    }
}

fn proj(u: &Vec<Rational>, v: &Vec<Rational>) -> (Vec<Rational>, Rational) {
    let mut ret = u.clone();
    let uv = inner_product(u, v);
    let uu = inner_product(u, u);

    let mu = uv / uu;
    for val in ret.iter_mut() {
        *val *= &mu;
    }

    return (ret, mu);
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

fn l2_norm_squared(v: &Vec<Rational>) -> Rational {
    return inner_product(v, v);
}

fn derivative(f: &Vec<Integer>, const_x: &Integer) -> Vec<Integer> {
    let mut to_ret = Vec::new();
    for i in 1..f.len() {
        let val = Integer::from(i as u32 * &f[i]) / const_x;
        // println!("{}", val);
        to_ret.push(val);
    }
    return to_ret;
}

pub fn approximate_zero(f: &Vec<Integer>, const_x: &Integer) -> Vec<Integer> {
    let debug = false;

    let mut results = Vec::new();
    let f_prime = derivative(f, const_x);

    if debug {
        println!("f' = {:?}", f_prime);
    }

    let n_guesses = 50;

    let xs: Vec<Rational> = (0..n_guesses)
        .map(|i| {
            Rational::from((
                2 * i * const_x.clone() - n_guesses * const_x.clone(),
                n_guesses,
            ))
        })
        .collect();

    if debug {
        println!("\n{:?}", xs);
    }

    for mut x in xs {
        let mut f_of_x = eval_rational_lattice_poly(&x, f, const_x);

        if debug {
            println!(
                " starting x to {:.3} \tf(x) = {:.1e}",
                x.to_f64(),
                f_of_x.to_f64()
            );
        }
        let mut count = 0;
        loop {
            let denom = eval_rational_lattice_poly(&x, &f_prime, const_x);
            if denom == 0 {
                println!("breaking");
                return results;
            }
            let to_sub = f_of_x / denom;
            x -= &to_sub;
            x = limit_precision(x, 4096);

            if to_sub <= (1, 256) {
                count += 1;
                if debug {
                    println!("newton finished in {} steps", count);
                    println!();
                }
                let int_x = x.round().into_numer_denom().0;
                if !results.contains(&int_x) {
                    results.push(int_x);
                }
                break;
            }
            f_of_x = eval_rational_lattice_poly(&x, f, const_x);
            count += 1;
            if debug {
                println!(
                    " trying x to {:.3} \tf(x) = {:.1e}",
                    x.to_f64(),
                    f_of_x.to_f64()
                );
            }
        }
    }

    return results;
}

fn limit_precision(mut x: Rational, shift: i32) -> Rational {
    x <<= shift;
    x.round_mut();
    x >>= shift;
    return x;
}

// fn limit_precision_mut(x: &mut Rational, shift: i32) {
//     *x <<= shift;
//     x.round_mut();
//     *x >>= shift;
// }
