use num_bigint::{BigInt, Sign};
use rug::integer::{IsPrime, Order};
use rug::rand::RandState;
use rug::{Assign, Integer};
use std::str;

fn main() {
    let mut i1 = BigInt::from(0);
    let i2 = BigInt::from_bytes_be(Sign::Plus, b"100000000");
    let i3 = BigInt::from_bytes_be(Sign::Plus, b"10000000000000000");

    i1 = &i1 + &i2 + &i3;
    println!("{} {} {}", i1, i2, i3);

    let N_BITS: u32 = 512;
    let e: Integer = Integer::from(65537);

    extended_euclidian(&Integer::from(240), &Integer::from(47));
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

    let d = find_inverse(&e, &phi_n);

    println!("e:{}\nd:{}", e, d);

    let product = Integer::from(&e * &d);
    let remainder = product % &n;
    println!("{}", remainder);

    let msg = "test message";
    let digits = msg.as_bytes();
    let msg_int = Integer::from_digits(digits, Order::Lsf);

    println!("{}", msg_int);

    let c = fast_power(&msg_int, &e, &n);

    println!("ciphertext: {}", c);
    
    let d = fast_power(&c, &d, &n);

    println!("recovered: {}", d);

    let recovered_digits = d.to_digits::<u8>(Order::Lsf);
    let recovered_msg = str::from_utf8(&recovered_digits).unwrap();

    println!("recovered message {:?}", recovered_msg);


}

fn find_inverse(e: &Integer, n: &Integer) -> Integer {
    let (_r, _s, t) = extended_euclidian(n, e);

    if t < 0 {
        let t = t + n;
        return t;
    } else {
        return t;
    }
}

fn extended_euclidian(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
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

fn fast_power(x: &Integer, e: &Integer, n: &Integer) -> Integer {
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
