// #![allow(dead_code, non_snake_case)]
// use ndarray::{Array2, Array3, ArrayView3};
use rug::integer::{IsPrime, Order};
// use rug::ops::Pow;
use rug::rand::RandState;
use rug::{Assign, Integer, Rational};
// use std::num::IntErrorKind;
use std::str::FromStr;
use std::time::Instant;
use std::{str, vec};

mod cryptlib;

fn main() {
    test_rsa();
    test_crt();
    test_gsp();
    test_gsp_equivalence();
    test_lll();
    test_coppersmith();
    test_hastad_broadcast();
    test_div_poly_zn();
    test_poly_euclid();
    test_franklin_reiter();
    test_determinant();
    test_short_pad();
}

fn test_short_pad() {
    let n_bits = 512;
    // let e = Integer::from(3);
    let e = Integer::from(3);

    let mut p = Integer::new();
    let mut q = Integer::new();

    let mut rand = RandState::new();
    rand.seed(&Integer::from(2));
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

    let m1 = (m.clone() << 32) + Integer::from(91418461);
    let m2 = (m.clone() << 32) + Integer::from(94392911);
    let diff = Integer::from(&m2 - &m1);

    let c1 = cryptlib::fast_power(&m1, &e, &n);
    let c2 = cryptlib::fast_power(&m2, &e, &n);

    let mut g1:Vec<Vec<Integer>> = (0..4).map(|_x| (0..4).map(|_y| Integer::from(0)).collect()).collect();
    g1[0][0] -= &c1; // x^3 - c1
    g1[3][0] += 1;

    let mut g2:Vec<Vec<Integer>> = (0..4).map(|_x| (0..4).map(|_y| Integer::from(0)).collect()).collect();
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
    let delta = cryptlib::coppersmith(&inverted, &n, 1, 18);
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
    let values = [[[5,0],[3,1],[2,1]],[[1,1],[4,0],[1,0]],[[3,0],[2,0],[3,1]]];
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
    let value = cryptlib::coppersmith(&f, &n, m, epsilon_denom);

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

        let config = HastadRSAConfig {
            n: n,
            f: f,
        };

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
    let x_0 = cryptlib::coppersmith(&g, &n, m, epsilon_denom);
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
