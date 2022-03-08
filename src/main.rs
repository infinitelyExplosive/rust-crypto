#![allow(dead_code, non_snake_case)]
use rug::integer::{IsPrime, Order};
use rug::rand::RandState;
use rug::{Assign, Integer, Rational};
use std::str::FromStr;
use std::{str, vec};
use std::time::Instant;

mod cryptlib;

fn main() {
    // test_rsa();
    // test_crt();
    // test_gsp();
    // test_gsp_equivalence();
    // test_lll();
    test_coppersmith();
    // test_hastad_broadcast();
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

    for line in other_data{
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
    for line in &b_star_2{
        print!(" ");
        for val in line {
            print!("{:+.3},  ", val.to_f32());
        }
        println!("")
    }

    for (line1, line2) in b_star_1.iter().zip(b_star_2.iter()) {
        for (val1, val2) in line1.iter().zip(line2){
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
    for line in &b_star_2{
        print!(" ");
        for val in line {
            print!("{:+.3},  ", val.to_f32());
        }
        println!("")
    }

    for (line1, line2) in b_star_1.iter().zip(b_star_2.iter()) {
        for (val1, val2) in line1.iter().zip(line2){
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
    let COEFFS: Vec<&str> = vec!["1942528644709637042", "1234567890123456789", "987654321987654321", "1"];
    let m = 2;
    let epsilon_denom = 10;

    for value in COEFFS {
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
    //TODO: FINISH ME
    struct HastadRSAConfig {
        n: Integer,
        e: Integer,
        d: Integer,
        f: Vec<Integer>,
    }

    let N_BITS = 32;
    let CONFIGS = 3; // degree (x+c)^3 = 3
    let e = Integer::from(2);

    let mut configs: Vec<HastadRSAConfig> = Vec::new();

    for i in 0..CONFIGS {
        let mut p = Integer::new();
        let mut q = Integer::new();

        let mut rand = RandState::new();
        while p.is_probably_prime(40) == IsPrime::No
            || Integer::from(&p % &e) == 0
            || configs.iter().any(|config| &config.n % p.clone() == 0)
        {
            p.assign(Integer::random_bits((N_BITS / 2)-4, &mut rand));
        }

        while q.is_probably_prime(40) == IsPrime::No
            || Integer::from(&q % &e) == 0
            || configs.iter().any(|config| &config.n % q.clone() == 0)
        {
            q.assign(Integer::random_bits((N_BITS / 2)-4, &mut rand));
        }
        println!("p:{} q:{}\nn:{}\n", p, q, Integer::from(&p*&q));

        let phi_n = Integer::from(&p - 1) * Integer::from(&q - 1);
        let d = cryptlib::find_inverse(&phi_n, &e);

        let mut f = Vec::new();

        f.push(Integer::from(32*i));
        f.push(Integer::from(1));

        let config = HastadRSAConfig {
            n: Integer::from(&p * &q),
            e: e.clone(),
            d: d,
            f: f,
        };

        configs.push(config);
    }

    // let msg = Integer::from_digits("a".as_bytes(), Order::Lsf);
    // let cs: Vec<Integer> = configs
    //     .iter()
    //     .map(|config| {
    //         let f_x = cryptlib::eval_poly(&msg, &config.f, &config.n);
    //         cryptlib::fast_power(&f_x, &e, &config.n)
    //     })
    //     .collect();

    // let gs: Vec<Vec<Integer>> = configs
    //     .iter()
    //     .zip(cs)
    //     .map(|(config, c)| {
    //         let mut g = cryptlib::exp_poly(&config.f, &e);
    //         g[0] -= &c;
    //         g
    //     })
    //     .collect();

    let n1s: Vec<Integer> = (0..CONFIGS)
        .map(|i| {
            (&configs[0..i])
                .iter()
                .chain(&configs[i + 1..CONFIGS])
                .fold(Integer::from(1), |acc, config| acc * &config.n)
        })
        .collect();

    let ts: Vec<Integer> = configs
        .iter()
        .zip(n1s)
        .map(|(config, n1)| {
            let (m1, _m2) = cryptlib::bezout(&n1, &config.n);
            Integer::from(&n1 * m1)
        })
        .collect();

    let mut g = Vec::new();
    g.push(Integer::from(0));
    g.push(Integer::from(0));
    g.push(Integer::from(0));

    for (config, t) in configs.iter().zip(ts.iter()) {
        for i in 0..config.f.len() {
            g[i] += Integer::from(&config.f[i] * t);
        }
    }

    let n = Integer::from(&ts[0] * &configs[0].n);
    println!("overall:");
    println!("n:{}", n);
    println!("g: {:?}", g);
    let x_0 = cryptlib::coppersmith(&g, &n, 10, 20);
    println!("{}", x_0);
}

fn test_lll() {
    let mut basis = Vec::new();

    // let data = vec![vec![1, 1, 1], vec![-1, 0, 2], vec![3, 5, 6]];
    let m = 10001;
    let x = 10;
    let data = vec![vec![m, 0, 0, 0], vec![0, m*x, 0, 0], vec![0, 0, m*x, 0], vec![-222, 5000*x, 10*x, x]];
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
    let N_BITS: u32 = 512;
    let e: Integer = Integer::from(65537);

    cryptlib::extended_euclidian(&Integer::from(240), &Integer::from(47));
    let mut p = Integer::new();
    let mut q = Integer::new();

    let mut rand = RandState::new();
    while p.is_probably_prime(40) == IsPrime::No || Integer::from(&p % &e) == 0 {
        p.assign(Integer::random_bits(N_BITS / 2, &mut rand));
    }

    while q.is_probably_prime(40) == IsPrime::No || Integer::from(&q % &e) == 0 {
        q.assign(Integer::random_bits(N_BITS / 2, &mut rand));
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
