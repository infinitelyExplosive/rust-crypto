use rug::Integer;

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
