#![allow(dead_code)]
use rug::integer::{IsPrime, Order};
use rug::rand::RandState;
use rug::{Assign, Integer, Rational};
use std::{str, vec};

mod cryptlib;

fn main() {
    // test_rsa();
    // test_crt();
    // test_gsp();
    test_lll();
}

fn test_lll() {
    let mut basis = Vec::new();

    let data = vec![vec![1,1,1],vec![-1,0,2],vec![3,5,6]];
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
    let data = vec![vec![3,5,-2,3,3], vec![1,2,3,4,5], vec![-5,5,-5,5,0], vec![3,0,2,0,1], vec![1,1,1,-1,1]];
    for line in data {
        let mut row = Vec::new();
        for val in line {
            row.push(Rational::from(val));
        }
        basis.push(row);
    }

    println!("{:?}", basis);
    let reduced= cryptlib::gsp(&basis);
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
