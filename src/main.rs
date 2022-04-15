#![allow(dead_code)]
// use ndarray::{Array2, Array3, ArrayView3};
use rug::integer::{IsPrime, Order};
use rug::ops::Pow;
use rug::rand::RandState;
use rug::{Assign, Integer, Rational};
use std::fs::File;
use std::io::BufRead;
// use std::num::IntErrorKind;
use std::str::FromStr;
use std::time::Instant;
use std::{io, str, vec};

mod cryptlib;
mod cryptlib_bv;

fn main() {
    // component tests
    // test_rsa();
    // test_crt();
    // test_gsp();
    // test_gsp_equivalence();
    // test_lll();
    // test_coppersmith();
    // test_div_poly_zn();
    // test_poly_euclid();
    // test_determinant();
    // test_inv_quad();
    // test_coppersmith_bv();
    // test_real_bv_polys();

    // attacks
    // test_hastad_broadcast();
    // test_franklin_reiter();
    test_short_pad();
    // test_partial_key();
}

fn test_real_bv_polys() {
    let data = File::open("polys.txt").unwrap();
    for line in io::BufReader::new(data).lines() {
        if let Ok(poly) = line {
            println!("{}", poly);
            let parts: Vec<i32> = poly.split(',').map(|x| x.parse::<i32>().unwrap()).collect();

            let mut f = Vec::new();
            let mut row = Vec::new();
            row.push(Integer::from(parts[0]));
            row.push(Integer::from(parts[2]));
            f.push(row);
            let mut row = Vec::new();
            row.push(Integer::from(parts[1]));
            row.push(Integer::from(parts[3]));
            f.push(row);

            let cap_x = Integer::from(parts[4] + 1);
            let cap_y = Integer::from(parts[5] + 1);

            if let Some((x_result, y_result)) = cryptlib_bv::coppersmith_bv(&f, &cap_x, &cap_y, 1) {
                println!("result {} {}", x_result, y_result);

                assert!(
                    cryptlib_bv::eval_poly_bv(&f, &x_result, &y_result, &Integer::from(-1)) == 0,
                    "f(x0, y0) != 0"
                );
                println!("{:-<1$}", "", 20);
            } else {
                println!("Failed to find root");
            }
        }
    }
}

fn test_coppersmith_bv() {
    let f: Vec<Vec<Integer>> = vec![vec![1, 433], vec![-28, 150]]
        .iter()
        .map(|row| row.iter().map(|val| Integer::from(*val)).collect())
        .collect();

    let cap_x = Integer::from(19);
    let cap_y = Integer::from(17);
    if let Some((x, y)) = cryptlib_bv::coppersmith_bv(&f, &cap_x, &cap_y, 1) {
        println!("{} {}", x, y);
    } else {
        println!("Failed to find root");
    }
}

fn test_inv_quad() {
    let goal = (Integer::from(9135123412323i64) * 8 + 1) * -1;
    let n = 512;
    let zero = Integer::from(0);
    let results = cryptlib::solve_quadratic(&zero, &zero, &goal, n);
    for result in results {
        println!("{}", result);
    }
}

fn test_partial_key() {
    let extra_bits = 14;// 1024:54, 512:28, 256:14, 128:8
    let n_bits = 256;

    // let p = Integer::from(4013);
    // let q = Integer::from(3851);
    // let p = Integer::from(2703229499i64);
    // let q = Integer::from(3240128051i64);

    let e = Integer::from(3);

    let mut p = Integer::new();
    let mut q = Integer::new();

    let mut rand = RandState::new();
    rand.seed(&Integer::from(116));
    while p.is_probably_prime(40) == IsPrime::No || Integer::from((p.clone() - 1) % &e) == 0 {
        p.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }

    while q.is_probably_prime(40) == IsPrime::No || Integer::from((q.clone() - 1) % &e) == 0 {
        q.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }

    let n = Integer::from(&p * &q);
    let phi_n = (p.clone() - 1) * (q.clone() - 1);
    let d = cryptlib::find_inverse(&e, &phi_n);

    let mut n_len = 0;
    let mut tmp = n.clone();
    while tmp > 0 {
        tmp >>= 1;
        n_len += 1;
    }
    let mask_len = if n_len % 4 == 0 {
        n_len / 4 + extra_bits
    } else {
        n_len / 4 + 1 + extra_bits
    };
    let mask = (Integer::from(1) << mask_len) - 1;
    let p0: Integer = p.clone() & &mask;
    let q0: Integer = q.clone() & &mask;
    let n0: Integer = n.clone() & &mask;
    let d0: Integer = d.clone() & &mask;
    let ed0: Integer = (e.clone() * &d0) & &mask;

    let now = Instant::now();

    println!("p: {:x} ({:x})\nq: {} ({})", p, p0, q, q0);
    println!("n: {:x} (len {})", n, n_len);
    println!("masklen: {}", mask_len);
    println!("n0 {:x}", n0);
    println!("d: {:x} ({:x})", d, d0);
    println!("ed mod 2^(n/4): {:x}", ed0);

    let mut p_candidates = Vec::new();

    for k in 0..(e.to_i32().unwrap()) {
        let mut f = Vec::new();
        f.push(n.clone() * k);
        f.push(ed0.clone() - Integer::from(&n * k) - k - 1);
        f.push(Integer::from(k));

        let modulus = Integer::from(2).pow(mask_len);
        let a = ((Integer::from(&f[2]) % &modulus) + &modulus) % &modulus;
        let b = ((Integer::from(&f[1]) % &modulus) + &modulus) % &modulus;
        let c = ((Integer::from(&f[0]) % &modulus) + &modulus) % &modulus;
        println!(
            " k: {}  {}x^2 + {}x + {} = 0 (mod 2^{})",
            k, a, b, c, mask_len
        );
        let mut candidates = cryptlib::solve_quadratic(&a, &b, &c, mask_len);
        candidates.sort();

        for candidate in candidates {
            let compliment = (n.clone()
                * cryptlib::find_inverse(&candidate, &Integer::from(2).pow(mask_len)))
                & &mask;
            let p_candidate = std::cmp::min(candidate, compliment);
            if !p_candidates.contains(&p_candidate) {
                p_candidates.push(p_candidate);
            }
        }
        println!(" p candidates {:?}", p_candidates);
    }
    // let l_correct = Float::with_val(128, &p).log2().ceil().to_integer().unwrap();
    // let k_correct: Float = Float::with_val(128, &n).log2() / 4;
    // let k_correct: Integer = k_correct.floor().to_integer().unwrap() + 1 + extra_bits;
    // println!("k correct {}  l correct {}", k_correct, l_correct);
    // let correct_x = p.clone() / Integer::from(2).pow(k_correct.to_u32().unwrap());
    // let correct_y = q.clone() / Integer::from(2).pow(k_correct.to_u32().unwrap());
    // println!("correct x {} correct y {}", correct_x, correct_y);

    let mut p_guess: Option<Integer> = None;
    let mut q_guess: Option<Integer> = None;
    'l_search: for l_offset in [0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 7, -7, 8, -8] {
        let l = (n_len as i32 / 2 + l_offset) as u32;
        println!("trying l {}", l);
        for (i, p0_guess) in p_candidates.iter().enumerate() {
            // let p0_guess = p0.clone();
            let q0_guess: Integer = (n.clone()
                * cryptlib::find_inverse(&p0_guess, &Integer::from(2).pow(mask_len)))
                & &mask;
            println!(" trying {} p0: {:x} q0: {:x}", i, p0_guess, q0_guess);

            let k = (n_len / 4) + 1 + extra_bits;
            let cap_x = Integer::from(2).pow(l - k);
            let cap_y = n.clone() / Integer::from(2).pow(l + k - 1);

            let mut f = Vec::new();
            let mut row = Vec::new();
            row.push((p0_guess.clone() * q0_guess.clone() - n.clone()) / Integer::from(2).pow(k));
            row.push(p0_guess.clone());
            f.push(row);
            let mut row = Vec::new();
            row.push(q0_guess.clone());
            row.push(Integer::from(2).pow(k));
            f.push(row);

            if let Some((x0, y0)) = cryptlib_bv::coppersmith_bv(&f, &cap_x, &cap_y, 1) {
                println!("  x0 {} y0 {}", x0, y0);
                p_guess = Some(Integer::from(2).pow(k) * x0 + p0_guess);
                q_guess = Some(Integer::from(2).pow(k) * y0 + &q0_guess);
                println!(
                    "\n  found p: {:x}\n        q: {:x}",
                    p_guess.as_ref().unwrap(),
                    q_guess.as_ref().unwrap()
                );
                if Integer::from(p_guess.as_ref().unwrap() * q_guess.as_ref().unwrap()) == n {
                    break 'l_search;
                }
            }
        }
    }

    if let (Some(recovered_p), Some(recovered_q)) = (p_guess, q_guess) {
        let recovered_phi = (recovered_p.clone() - 1) * (recovered_q.clone() - 1);
        let recovered_d = cryptlib::find_inverse(&e, &recovered_phi);
        let duration = now.elapsed();
        println!("recovered d {:x}", recovered_d);
        assert!(recovered_d == d);
        let mins = duration.as_secs() / 60;
        let secs = duration.as_secs() % 60;
        println!("in {} minutes {} seconds", mins, secs);
    } else {
        println!("failed");
    }
}

fn test_short_pad() {
    let n_bits = 512;
    // let e = Integer::from(3);
    let e = Integer::from(3);

    let mut p = Integer::new();
    let mut q = Integer::new();

    let mut rand = RandState::new();
    rand.seed(&Integer::from(1));
    while p.is_probably_prime(40) == IsPrime::No || Integer::from(&p % &e) == 0 {
        p.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }

    while q.is_probably_prime(40) == IsPrime::No || Integer::from(&q % &e) == 0 {
        q.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }
    let n = Integer::from(&p * &q);
    println!("n:{}", n);

    let m = Integer::from_digits("YELLOW SUBMARINE".as_bytes(), Order::Lsf);
    // let m = Integer::from(211601);

    let m1 = (m.clone() << 32) + Integer::from(12461);
    let m2 = (m.clone() << 32) + Integer::from(28911);
    let diff = Integer::from(&m2 - &m1);

    let c1 = cryptlib::fast_power(&m1, &e, &n);
    let c2 = cryptlib::fast_power(&m2, &e, &n);

    let mut g1: Vec<Vec<Integer>> = (0..4)
        .map(|_x| (0..4).map(|_y| Integer::from(0)).collect())
        .collect();
    g1[0][0] -= &c1; // x^3 - c1
    g1[3][0] += 1;

    let mut g2: Vec<Vec<Integer>> = (0..4)
        .map(|_x| (0..4).map(|_y| Integer::from(0)).collect())
        .collect();
    g2[0][0] -= &c2; // (x+y)^3 - c2
    g2[0][3] += 1;
    g2[1][2] += 3;
    g2[2][1] += 3;
    g2[3][0] += 1;

    println!("{:?}", g1);
    println!("{:?}", g2);
    println!();
    let resultant = cryptlib::resultant(&g1, &g2, &n);
    println!("{:?}", resultant);

    let inverted = resultant.iter().map(|x| Integer::from(-x)).collect();
    let result_inverted = cryptlib::eval_poly(&diff, &inverted, &n);
    println!("sanity: resultant({}) = {}", diff, result_inverted);
    let delta = cryptlib::coppersmith(&inverted, &n, 1, 18).unwrap();
    println!("delta {}", delta);

    let mut f = Vec::new();
    f.push(-delta);
    f.push(Integer::from(1));

    let mut fr_g1 = cryptlib::exp_poly(&f, &e);
    fr_g1[0] -= &c1;

    for val in fr_g1.iter_mut() {
        *val %= &n;
        *val += &n;
        *val %= &n;
    }
    let mut identity_func = Vec::new();
    identity_func.push(Integer::from(0));
    identity_func.push(Integer::from(1));
    let mut fr_g2 = cryptlib::exp_poly(&identity_func, &e);
    fr_g2[0] -= &c2;
    for val in fr_g2.iter_mut() {
        *val %= &n;
        *val += &n;
        *val %= &n;
    }

    println!("g1: {:?}\ng2: {:?}", fr_g1, fr_g2);

    let (mut r, _s, _t) = cryptlib::poly_extended_euclidean_zn(&fr_g1, &fr_g2, &n);
    // println!("r: {:?}\ns: {:?}\nt: {:?}", r, s, t);

    r[0] += &n;
    println!("r {:?}", r);
    let inv_x_term = cryptlib::find_inverse(&r[1], &n);
    let recovered_m2: Integer = (((Integer::from(&inv_x_term * -1) * &r[0] % &n) + &n) % &n) >> 32;

    println!("recovered m2: {}", recovered_m2);

    let mut msg_bytes = Vec::new();
    for offset in 0..=(recovered_m2.significant_bits() / 8) {
        let low_bytes: Integer = recovered_m2.clone() >> (offset * 8) & 0xff;
        let low_u8 = low_bytes.to_u8().unwrap();
        msg_bytes.push(low_u8);
    }
    println!("{}", String::from_utf8(msg_bytes).unwrap());
}

fn test_determinant() {
    let mut matrix = Vec::new();
    let values = [
        [[5, 0], [3, 1], [2, 1]],
        [[1, 1], [4, 0], [1, 0]],
        [[3, 0], [2, 0], [3, 1]],
    ];
    // let values = [[[1,1],[3,1]],[[2,1],[5,1]]];

    for row in values {
        let mut row_vec = Vec::new();
        for elem in row {
            let mut elem_vec = Vec::new();
            for y_term in elem {
                elem_vec.push(Integer::from(y_term));
            }
            row_vec.push(elem_vec);
        }
        matrix.push(row_vec);
    }

    for row in &matrix {
        println!("{:?}", row);
    }
    let cols = (0..matrix.len()).collect();
    let result = cryptlib::determinant(&matrix, &cols, 5, &Integer::from(97));
    println!("{:?}", result);
}

fn test_div_poly_zn() {
    let n_bits = 256;
    let e = Integer::from(3);

    let mut p = Integer::new();
    let mut q = Integer::new();

    let mut rand = RandState::new();
    rand.seed(&Integer::from(3));
    while p.is_probably_prime(40) == IsPrime::No || Integer::from(&p % &e) == 0 {
        p.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }

    while q.is_probably_prime(40) == IsPrime::No || Integer::from(&q % &e) == 0 {
        q.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }
    // let n = Integer::from(&p * &q);
    let n = Integer::from(7);
    println!("{}", n);

    let mut f = Vec::new();
    f.push(Integer::from(5));
    f.push(Integer::from(4));
    f.push(Integer::from(6));
    f.push(Integer::from(2));
    f.push(Integer::from(4));

    let mut g = Vec::new();
    g.push(Integer::from(2));
    g.push(Integer::from(0));
    g.push(Integer::from(3));

    let (mut r, q) = cryptlib::divide_poly_zn(&f, &g, &n);

    for elem in r.iter_mut() {
        *elem += &n;
        *elem %= &n;
    }

    println!("{:?}", r);
    println!("{:?}", q);

    let mut prod = cryptlib::multiply_poly_zn(&r, &g, &n);
    prod[0] += &q[0];
    prod[1] += &q[1];

    for elem in prod.iter_mut() {
        *elem %= &n;
        *elem += &n;
        *elem %= &n;
    }

    println!("prod {:?}", prod);
}

fn test_poly_euclid() {
    let p = Integer::from(47);
    let q = Integer::from(67);

    let n = Integer::from(&p * &q);

    let e = Integer::from(3);

    println!("n: {} e: {}", n, e);

    let mut f = Vec::new();
    f.push(Integer::from(5));
    f.push(Integer::from(3));

    let m2 = Integer::from(14);
    let m1 = cryptlib::eval_poly(&m2, &f, &n);

    println!("m1: {} m2: {}", m1, m2);

    let c1 = cryptlib::fast_power(&m1, &e, &n);
    let c2 = cryptlib::fast_power(&m2, &e, &n);

    println!("c1: {} c2: {}", c1, c2);

    let mut g1 = cryptlib::exp_poly(&f, &e);
    g1[0] -= &c1;
    g1[0] %= &n;
    g1[0] += &n;
    g1[0] %= &n;

    let mut g2 = Vec::new();
    g2.push(n.clone() - &c2);
    g2.push(Integer::from(0));
    g2.push(Integer::from(0));
    g2.push(Integer::from(1));

    println!("g1: {:?}", g1);
    println!("g2: {:?}", g2);

    let (r, s, t) = cryptlib::poly_extended_euclidean_zn(&g1, &g2, &n);
    println!("{:?}", r);
    println!("{:?}", s);
    println!("{:?}", t);

    let inv_x_term = cryptlib::find_inverse(&r[1], &n);
    let recovered_m2 = Integer::from(&inv_x_term * -1) * &r[0] % &n;

    println!("recovered m2: {}", recovered_m2);
}

fn test_franklin_reiter() {
    let n_bits = 2048;
    let e = Integer::from(3);

    let mut p = Integer::new();
    let mut q = Integer::new();

    let mut rand = RandState::new();
    rand.seed(&Integer::from(3));
    while p.is_probably_prime(40) == IsPrime::No || Integer::from(&p % &e) == 0 {
        p.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }

    while q.is_probably_prime(40) == IsPrime::No || Integer::from(&q % &e) == 0 {
        q.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }
    let n = Integer::from(&p * &q);
    println!("p:{} q:{}\nn:{}\n", p, q, n);

    let phi_n = Integer::from(&p - 1) * Integer::from(&q - 1);
    let _d = cryptlib::find_inverse(&phi_n, &e);

    let mut f = Vec::new();
    f.push(Integer::from(20));
    f.push(Integer::from(3));

    let msg2 = Integer::from_digits("really long message takes no time because this is a direct math method with no sampling (boring)".as_bytes(), Order::Lsf);
    let msg1 = cryptlib::eval_poly(&msg2, &f, &n);

    println!("msg1: {}\nmsg2: {}", msg1, msg2);

    let c1 = cryptlib::fast_power(&msg1, &e, &n);
    let c2 = cryptlib::fast_power(&msg2, &e, &n);

    println!("c1: {}\nc2: {}", c1, c2);

    let mut g1 = cryptlib::exp_poly(&f, &e);
    g1[0] -= &c1;
    for val in g1.iter_mut() {
        *val %= &n;
        *val += &n;
        *val %= &n;
    }
    let mut identity_func = Vec::new();
    identity_func.push(Integer::from(0));
    identity_func.push(Integer::from(1));
    let mut g2 = cryptlib::exp_poly(&identity_func, &e);
    g2[0] -= &c2;
    for val in g2.iter_mut() {
        *val %= &n;
        *val += &n;
        *val %= &n;
    }

    println!("g1: {:?}\ng2: {:?}", g1, g2);

    let (r, _s, _t) = cryptlib::poly_extended_euclidean_zn(&g1, &g2, &n);
    // println!("r: {:?}\ns: {:?}\nt: {:?}", r, s, t);

    let inv_x_term = cryptlib::find_inverse(&r[1], &n);
    let recovered_m2 = Integer::from(&inv_x_term * -1) * &r[0] % &n;

    println!("recovered m2: {}", recovered_m2);

    let mut msg_bytes = Vec::new();
    for offset in 0..=(recovered_m2.significant_bits() / 8) {
        let low_bytes: Integer = recovered_m2.clone() >> (offset * 8) & 0xff;
        let low_u8 = low_bytes.to_u8().unwrap();
        msg_bytes.push(low_u8);
    }
    println!("{}", String::from_utf8(msg_bytes).unwrap());
}

fn test_gsp_equivalence() {
    let mut basis = Vec::new();
    let mut other_basis = Vec::new();
    let data = vec![
        vec![3, 5, -2, 3, 3],
        vec![1, 2, 3, 4, 5],
        vec![-5, 5, -5, 5, 0],
        vec![3, 0, 2, 0, 1],
        vec![1, 1, 1, -1, 1],
    ];

    let other_data = vec![
        vec![3, 5, -2, 3, 3],
        vec![1, 2, 3, 4, 5],
        vec![-4, 4, -6, 6, 0],
        vec![3, 0, 2, 0, 1],
        vec![1, 1, 1, -1, 1],
    ];

    for line in data {
        let mut row = Vec::new();
        for val in line {
            row.push(Rational::from(val));
        }
        basis.push(row);
    }

    for line in other_data {
        let mut row = Vec::new();
        for val in line {
            row.push(Rational::from(val));
        }
        other_basis.push(row);
    }

    println!("{:?}", basis);
    let (b_star_1, _mu_matrix) = cryptlib::gsp(&basis);
    println!("gsp:");
    for line in &b_star_1 {
        print!(" ");
        for val in line {
            print!("{:+.3},  ", val.to_f32());
        }
        println!("")
    }

    let mut b_star_2 = (0..basis.len()).map(|_x| Vec::new()).collect();
    let mut mu_matrix = (0..basis.len()).map(|_x| Vec::new()).collect();
    cryptlib::gsp_efficient(&basis, &mut b_star_2, &mut mu_matrix, 0);

    println!("gsp2:");
    for line in &b_star_2 {
        print!(" ");
        for val in line {
            print!("{:+.3},  ", val.to_f32());
        }
        println!("")
    }

    for (line1, line2) in b_star_1.iter().zip(b_star_2.iter()) {
        for (val1, val2) in line1.iter().zip(line2) {
            assert!(*val1 == *val2);
        }
    }

    println!();

    let (b_star_1, _mu_matrix) = cryptlib::gsp(&other_basis);
    println!("gsp:");
    for line in &b_star_1 {
        print!(" ");
        for val in line {
            print!("{:+.3},  ", val.to_f32());
        }
        println!("")
    }

    cryptlib::gsp_efficient(&other_basis, &mut b_star_2, &mut mu_matrix, 2);

    println!("gsp2:");
    for line in &b_star_2 {
        print!(" ");
        for val in line {
            print!("{:+.3},  ", val.to_f32());
        }
        println!("")
    }

    for (line1, line2) in b_star_1.iter().zip(b_star_2.iter()) {
        for (val1, val2) in line1.iter().zip(line2) {
            assert!(*val1 == *val2);
        }
    }
}

fn test_coppersmith() {
    let mut f = Vec::new();

    // let COEFFS = vec![5609315825568i64, -18680690149i64, 18544119, -7299, 1]; // 511, 1807, 2133, 2848
    // let n = Integer::from_str("117129523791978766508").unwrap();
    // let epsilon_denom = 10;
    // let m = 10;

    // let COEFFS = vec![-1194078866, 6379973, -5620, 1]; // 233, 1234, 4153
    // let n = Integer::from(541_327_006_526_i64);
    // let epsilon_denom = 8;
    // let m = 8;

    // let m = 6; //O?(10/2)
    // let epsilon_denom = 6; // 38 = x ^ (1/2 - 1/6) = 54872
    // let n = Integer::from(100000);
    // let COEFFS = vec![26233,-746,1]; // 37, 709

    // let m = 5;
    // let epsilon_denom = 5;
    // let n = Integer::from(35000); // 23 = x ^ (1/2 - 1/5) = 35000
    // let COEFFS = vec![3013, -154, 1]; // 23, 131

    let p: i64 = 1073741827;
    let q: i64 = 4294967311;
    let n = Integer::from(p) * q;
    let coeffs: Vec<&str> = vec![
        "1942528644709637042",
        "1234567890123456789",
        "987654321987654321",
        "1",
    ];
    let m = 2;
    let epsilon_denom = 10;

    for value in coeffs {
        // f.push(Integer::from(value));
        f.push(Integer::from_str(value).unwrap());
    }

    println!("f is {:?}", f);
    let now = Instant::now();
    let value = cryptlib::coppersmith(&f, &n, m, epsilon_denom).unwrap();

    let duration = now.elapsed();
    if value > 0 {
        println!("solution is {}", value);
    }
    let mins = duration.as_secs() / 60;
    let secs = duration.as_secs() % 60;
    println!("in {} minutes {} seconds", mins, secs);
}

fn test_hastad_broadcast() {
    struct HastadRSAConfig {
        n: Integer,
        f: Vec<Integer>,
    }

    let n_bits = 256;
    let num_configs = 3; // degree (x+c)^3 = 3
    let m = 2;
    let epsilon_denom = 11;
    let e = Integer::from(3);
    let msg = Integer::from_digits("YELLOW SUBMARINE".as_bytes(), Order::Lsf);

    let mut configs: Vec<HastadRSAConfig> = Vec::new();

    let mut rand = RandState::new();
    // rand.seed(&Integer::from(1));
    for i in 0..num_configs {
        let mut p = Integer::new();
        let mut q = Integer::new();

        while p.is_probably_prime(40) == IsPrime::No
            || Integer::from(&p % &e) == 0
            || configs.iter().any(|config| &config.n % p.clone() == 0)
        {
            p.assign(Integer::random_bits(n_bits / 2, &mut rand));
        }

        while q.is_probably_prime(40) == IsPrime::No
            || Integer::from(&q % &e) == 0
            || configs.iter().any(|config| &config.n % q.clone() == 0)
        {
            q.assign(Integer::random_bits(n_bits / 2, &mut rand));
        }
        let n = Integer::from(&p * &q);
        println!("p:{} q:{}\nn:{}\n", p, q, n);

        let mut f = Vec::new();

        f.push(Integer::from(32 * i));
        f.push(Integer::from(1));
        // println!("f: {:?}\n", f);

        let config = HastadRSAConfig { n: n, f: f };

        configs.push(config);
    }

    println!("msg:{}", msg);

    let cs: Vec<Integer> = configs
        .iter()
        .map(|config| {
            let f_x = cryptlib::eval_poly(&msg, &config.f, &config.n);
            cryptlib::fast_power(&f_x, &e, &config.n)
        })
        .collect();

    let gs: Vec<Vec<Integer>> = configs
        .iter()
        .zip(&cs)
        .map(|(config, c)| {
            let mut g = cryptlib::exp_poly(&config.f, &e);
            g[0] -= c;
            g
        })
        .collect();

    let n1s: Vec<Integer> = (0..num_configs)
        .map(|i| {
            (&configs[0..i])
                .iter()
                .chain(&configs[i + 1..num_configs])
                .fold(Integer::from(1), |acc, config| acc * &config.n)
        })
        .collect();

    let ts: Vec<Integer> = configs
        .iter()
        .zip(&n1s)
        .map(|(config, n1)| {
            let (m1, _m2) = cryptlib::bezout(n1, &config.n);
            Integer::from(n1 * m1)
        })
        .collect();

    let n = Integer::from(&n1s[0] * &configs[0].n);
    let mut g = Vec::new();
    for _ in 0..(e.to_i32().unwrap() + 1) {
        g.push(Integer::from(0));
    }

    for (g_partial, t) in gs.iter().zip(ts.iter()) {
        for i in 0..g_partial.len() {
            g[i] += Integer::from(&g_partial[i] * t);
            g[i] %= &n;
        }
    }

    println!("overall:");
    println!("n:{}", n);
    println!("g: {:?}", g);
    println!("sanity: {}", cryptlib::eval_poly(&msg, &g, &n));

    let now = Instant::now();
    let x_0 = cryptlib::coppersmith(&g, &n, m, epsilon_denom).unwrap();
    let duration = now.elapsed();
    println!("{}", x_0);
    let mut msg_bytes = Vec::new();
    for offset in 0..=(x_0.significant_bits() / 8) {
        let low_bytes: Integer = x_0.clone() >> (offset * 8) & 0xff;
        let low_u8 = low_bytes.to_u8().unwrap();
        msg_bytes.push(low_u8);
    }
    println!("{}", String::from_utf8(msg_bytes).unwrap());
    let mins = duration.as_secs() / 60;
    let secs = duration.as_secs() % 60;
    println!("in {} minutes {} seconds", mins, secs);
}

fn test_lll() {
    let mut basis = Vec::new();

    // let data = vec![vec![1, 1, 1], vec![-1, 0, 2], vec![3, 5, 6]];
    let m = 10001;
    let x = 10;
    let data = vec![
        vec![m, 0, 0, 0],
        vec![0, m * x, 0, 0],
        vec![0, 0, m * x, 0],
        vec![-222, 5000 * x, 10 * x, x],
    ];
    // let data = vec![vec![m, 0, 0, 0], vec![0, m, 0, 0], vec![0, 0, m, 0], vec![-222, 5000, 10, 0]];
    for line in data {
        let mut row = Vec::new();
        for val in line {
            row.push(Integer::from(val));
        }
        basis.push(row);
    }
    println!("{:?}", basis);
    let reduced = cryptlib::lll(&basis);
    println!("lll:\n{:?}", reduced);
}

fn test_gsp() {
    let mut basis = Vec::new();
    // let data = vec![vec![1,-1,1], vec![1,0,1], vec![1,1,2]];
    let data = vec![
        vec![3, 5, -2, 3, 3],
        vec![1, 2, 3, 4, 5],
        vec![-5, 5, -5, 5, 0],
        vec![3, 0, 2, 0, 1],
        vec![1, 1, 1, -1, 1],
    ];
    for line in data {
        let mut row = Vec::new();
        for val in line {
            row.push(Rational::from(val));
        }
        basis.push(row);
    }

    println!("{:?}", basis);
    let reduced = cryptlib::gsp(&basis);
    println!("gsp:\n{:?}", reduced);
    // https://www.emathhelp.net/en/calculators/linear-algebra/gram-schmidt-calculator/?i=%5B%5B3%2C1%2C-5%2C3%2C1%5D%2C%5B5%2C2%2C5%2C0%2C1%5D%2C%5B-2%2C3%2C-5%2C2%2C1%5D%2C%5B3%2C4%2C5%2C0%2C-1%5D%2C%5B3%2C5%2C0%2C1%2C1%5D%5D
}

fn test_crt() {
    let mut values = Vec::new();
    // values.push(Integer::from(0));
    // values.push(Integer::from(3));
    // values.push(Integer::from(4));
    values.push(Integer::from(6));
    values.push(Integer::from(13));
    values.push(Integer::from(9));
    values.push(Integer::from(19));

    let mut mods = Vec::new();
    // mods.push(Integer::from(3));
    // mods.push(Integer::from(4));
    // mods.push(Integer::from(5));
    mods.push(Integer::from(11));
    mods.push(Integer::from(16));
    mods.push(Integer::from(21));
    mods.push(Integer::from(25));

    let x = cryptlib::crt(values.iter(), mods.iter());
    println!("{}", x);
    println!("{} mod 11", Integer::from(&x % 11));
    println!("{} mod 16", Integer::from(&x % 16));
    println!("{} mod 21", Integer::from(&x % 21));
    println!("{} mod 25", Integer::from(&x % 25));
}

fn test_rsa() {
    let n_bits: u32 = 512;
    let e: Integer = Integer::from(65537);

    cryptlib::extended_euclidean(&Integer::from(240), &Integer::from(47));
    let mut p = Integer::new();
    let mut q = Integer::new();

    let mut rand = RandState::new();
    while p.is_probably_prime(40) == IsPrime::No || Integer::from(&p % &e) == 0 {
        p.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }

    while q.is_probably_prime(40) == IsPrime::No || Integer::from(&q % &e) == 0 {
        q.assign(Integer::random_bits(n_bits / 2, &mut rand));
    }
    println!("p:{}\nq:{}", p, q);

    let n = Integer::from(&p * &q);
    let phi_n = Integer::from(&p - 1) * &(Integer::from(&q - 1));
    println!("     N:{}\nphi(N):{}", n, phi_n);

    let d = cryptlib::find_inverse(&e, &phi_n);

    println!("e:{}\nd:{}", e, d);

    let product = Integer::from(&e * &d);
    let remainder = product % &phi_n;
    println!("remainder (should be 1):{}", remainder);

    let msg = "test message";
    let digits = msg.as_bytes();
    let msg_int = Integer::from_digits(digits, Order::Lsf);

    println!("msg as int: {}", msg_int);

    let c = cryptlib::fast_power(&msg_int, &e, &n);

    println!("ciphertext: {}", c);

    let d = cryptlib::fast_power(&c, &d, &n);

    println!("recovered: {}", d);

    let recovered_digits = d.to_digits::<u8>(Order::Lsf);
    let recovered_msg = str::from_utf8(&recovered_digits).unwrap();

    println!("recovered message {:?}", recovered_msg);
}

fn print_binary(x: &Integer, indent: u32) {
    let mut output = String::new();
    for _ in 0..indent {
        print!(" ")
    }
    let mut x = x.clone();
    while x > 0 {
        output.insert_str(0, &format!("{:04b} ", x.clone() & 0xf));
        x >>= 4;
    }
    print!("{}", output);
}
