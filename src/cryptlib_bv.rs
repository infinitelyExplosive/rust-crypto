use crate::cryptlib;
use rug::{integer::IsPrime, ops::Pow, Assign, Float, Integer, Rational};

pub fn eval_poly_bv(f: &Vec<Vec<Integer>>, x: &Integer, y: &Integer, n: &Integer) -> Integer {
    let mut sum = Integer::from(0);
    let mut x_power = Integer::from(1);
    for row in f {
        let mut y_power = Integer::from(1);
        // println!("row {:?}", row);
        for term in row {
            let part = Integer::from(&x_power * &y_power) * term;
            // println!(" {} (x{} y{} term{})", part, x_power, y_power, term);
            sum += part;
            y_power *= y;
            if *n > 0 {
                sum %= n;
                y_power %= n;
            }
        }
        x_power *= x;
        if *n > 0 {
            x_power %= n;
        }
    }
    return sum;
}
/// Returns the x-first leading coefficient as (coefficient, x power, y power)
fn lead_bv(f: &Vec<Vec<Integer>>) -> (Integer, usize, usize) {
    for (i, row) in f.iter().enumerate().rev() {
        for (j, val) in row.iter().enumerate().rev() {
            if *val != 0 {
                return (val.clone(), i, j);
            }
        }
    }
    return (Integer::from(0), 0, 0);
}

/// Returns the degree of the polynomial as (x degree, y degree)
fn degree_bv(f: &Vec<Vec<Integer>>) -> (usize, usize) {
    let mut xmax = 0;
    let mut ymax = 0;
    for (i, row) in f.iter().enumerate() {
        for (j, val) in row.iter().enumerate() {
            if *val != 0 {
                if i > xmax {
                    xmax = i;
                }
                if j > ymax {
                    ymax = j;
                }
            }
        }
    }
    return (xmax, ymax);
}

/// Returns the result of multiplying two polynomials mod n, or over integers if n=-1
pub fn multiply_poly_bv(
    f: &Vec<Vec<Integer>>,
    g: &Vec<Vec<Integer>>,
    n: &Integer,
) -> Vec<Vec<Integer>> {
    let debug = false;
    if debug {
        print!(" mult {:?}*{:?}", f, g);
    }
    let mut prod = Vec::new();
    let (f_dx, f_dy) = degree_bv(f);
    let (g_dx, g_dy) = degree_bv(g);
    for _ in 0..(f_dx + g_dx + 1) {
        let mut row = Vec::new();
        for _ in 0..(f_dy + g_dy + 1) {
            row.push(Integer::from(0));
        }
        prod.push(row);
    }

    for i in 0..(f_dx + 1) {
        for j in 0..(g_dx + 1) {
            for u in 0..(f_dy + 1) {
                for v in 0..(g_dy + 1) {
                    if *n < 0 {
                        prod[i + j][u + v] += Integer::from(&f[i][u] * &g[j][v]);
                    } else {
                        prod[i + j][u + v] += Integer::from(&f[i][u] * &g[j][v]) % n;
                    }
                }
            }
        }
    }
    if debug {
        println!("\t|\t   result {:?}", prod);
    }
    return prod;
}

/// Returns f to the e mod n, or over integers if n=-1
pub fn exp_poly_bv(f: &Vec<Vec<Integer>>, e: &Integer, n: &Integer) -> Vec<Vec<Integer>> {
    if *e == 0 {
        let (dx, dy) = degree_bv(f);
        let mut result: Vec<Vec<Integer>> = (0..=dx)
            .map(|_i| (0..=dy).map(|_j| Integer::from(0)).collect())
            .collect();
        result[0][0].assign(1);
        return result;
    }
    let mut result = f.clone();

    for _ in 0..Integer::from(e - 1).to_i32().unwrap() {
        result = multiply_poly_bv(&result, f, n);
    }
    return result;
}

/// Finds a small root (x0, y0) of f where x0 < X, y0 < Y
/// larger values of k allow for larger values of X, Y
pub fn coppersmith_bv(
    f: &Vec<Vec<Integer>>,
    cap_x: &Integer,
    cap_y: &Integer,
    k: usize,
) -> Option<(Integer, Integer)> {
    let debug = false;

    //remove gcd if present
    let mut gcd = Integer::from(0);
    for row in f {
        for elem in row {
            gcd = cryptlib::gcd(&gcd, elem);
        }
    }
    let mut f = f.clone();
    if gcd > 1 {
        if debug {
            println!("changing f from\n{:?}", f);
        }
        for row in f.iter_mut() {
            for elem in row.iter_mut() {
                *elem /= &gcd;
            }
        }
    }
    if debug {
        println!("to\n{:?}", f);
    }

    let mut cap_x = cap_x.clone();
    let mut cap_y = cap_y.clone();

    // adjust X,Y if gcd(p_0,0, XY) != 1
    if cryptlib::gcd(&f[0][0], &Integer::from(&cap_x * &cap_y)) != 1 {
        let mut alt_x: Integer = cap_x.clone() + 1;
        while f[0][0].clone() % &alt_x == 0 || alt_x.is_probably_prime(40) == IsPrime::No {
            alt_x += 1;
        }
        let mut alt_y: Integer = cap_y.clone() + 1;
        while f[0][0].clone() % &alt_y == 0
            || alt_y == alt_x
            || alt_y.is_probably_prime(40) == IsPrime::No
        {
            alt_y += 1;
        }
        cap_x = alt_x;
        cap_y = alt_y;
        assert!(
            cryptlib::gcd(&f[0][0], &Integer::from(&cap_x * &cap_y)) == 1,
            "failed to find new X, Y (this should never happen)"
        );
    }

    // unsure
    let epsilon_recip = 10;

    let (dx, dy) = degree_bv(&f);
    let delta = std::cmp::max(dx, dy);

    // W = max a_ij X^i Y^j
    let cap_w = f
        .iter()
        .enumerate()
        .map(|(i, row)| {
            row.iter()
                .enumerate()
                .map(|(j, elem)| {
                    let val =
                        (cap_x.clone().pow(i as u32) * cap_y.clone().pow(j as u32) * elem).abs();
                    if debug {
                        println!("W{},{} = {} from {}", i, j, val, elem);
                    }
                    val
                })
                .max()
                .unwrap()
        })
        .max()
        .unwrap();

    let omega = (delta + k + 1).pow(2);

    let u = ((Integer::from(1 - &cap_w) % f[0][0].clone().abs()) + f[0][0].clone().abs())
        % f[0][0].clone().abs()
        + &cap_w;

    let lhs = Float::with_val(512, omega).sqrt() / Integer::from(2).pow(omega as u32) * &cap_w;
    let rhs = cap_w.clone() * 2;

    assert!(lhs <= u, "sqrtw * 2^(-w) * W > u");
    assert!(u < rhs, "u >= 2W");

    assert!(cryptlib::gcd(&u, &f[0][0]) == 1, "gcd u,p00 != 1");

    let n = Integer::from(&cap_x * &cap_y).pow(k as u32) * &u;

    assert!(cryptlib::gcd(&n, &f[0][0]) == 1, "gcd n, p00 != 1");

    // check sqrt(ω) * 2^(-ω) * (XY)^k * W <= n < 2 * (KY)^k * W
    let lhs = Float::with_val(512, omega).root(2) / Integer::from(2).pow(omega as u32)
        * (cap_x.clone() * cap_y.clone()).pow(k as u32)
        * cap_w.clone();
    let rhs = 2 * (cap_x.clone() * cap_y.clone()).pow(k as u32) * cap_w.clone();
    assert!(lhs <= n, "sqrtw * 2^(-w)(XY)^k * W > n");
    assert!(n < rhs, "n >= 2(XY)^k * W");

    if debug {
        print!("f");
        print_poly_bv(&f);
        println!(
            "X = {}\nY = {}\nk = {}\nW = {}\n𝛿 = {}\nω = {}\nu = {}\nn = {}",
            cap_x, cap_y, k, cap_w, delta, omega, u, n
        );
        println!();
    }

    let lhs = Float::with_val(1024, cap_x.clone() * &cap_y);
    let w_pow = (2 * epsilon_recip - 3 * delta) as u32;
    let w_root = (3 * epsilon_recip * delta) as u32;
    let rhs = Float::with_val(1024, cap_w.clone()).pow(w_pow).root(w_root);

    println!(" condition {}: {:.4}", lhs < rhs, (lhs / rhs).log10());

    assert!(f[0][0] != 0, "p_0,0 = 0");

    let inverse = cryptlib::find_inverse(&f[0][0], &n);
    let test = (((f[0][0].clone() * &inverse) % &n) + &n) % &n;
    // println!("inverse {} test {}", inverse, test);
    assert!(test == 1, "f_0,0^-1 not inverse");

    let q: Vec<Vec<Integer>> = f
        .iter()
        .map(|row| {
            row.iter()
                .map(|elem| (((elem.clone() * &inverse) % &n) + &n) % &n)
                .collect()
        })
        .collect();

    assert!(q[0][0] == 1, "q_0,0 != 1");

    if debug {
        print!("q");
        print_poly_bv(&q);
    }

    // create q
    let mut qs: Vec<Vec<Vec<Integer>>> = (0..=(k + delta))
        .map(|i| {
            (0..=(k + delta))
                .map(|j| {
                    let mut new_q: Vec<Vec<Integer>> = (0..=(dx + k))
                        .map(|_| (0..=(dy + k)).map(|_| Integer::from(0)).collect())
                        .collect();
                    if i <= k && j <= k {
                        let cap_x_pow = cap_x.clone().pow((k - i) as u32);
                        let cap_y_pow = cap_y.clone().pow((k - j) as u32);
                        if debug {
                            print!("q_{},{} with X^: {} y^: {}  [", i, j, cap_x_pow, cap_y_pow);
                        }
                        for (x, row) in q.iter().enumerate() {
                            for (y, elem) in row.iter().enumerate() {
                                if *elem != 0 {
                                    new_q[x + i][y + j].assign(elem);
                                    new_q[x + i][y + j] *= &cap_x_pow;
                                    new_q[x + i][y + j] *= &cap_y_pow;
                                    if debug {
                                        print!("{}.x{}.y{}, ", new_q[x + i][y + j], x + i, y + j);
                                    }
                                }
                            }
                        }
                        if debug {
                            println!("]");
                        }
                    } else {
                        new_q[i][j].assign(&n);
                    }
                    let new_q = new_q.concat();
                    new_q
                })
                .collect()
        })
        .collect();

    // create q~
    for i in 0..=(k + delta) {
        for j in 0..=(k + delta) {
            let func_idx = i * (k + delta + 1) + j;
            let cap_x_pow = cap_x.clone().pow(i as u32);
            let cap_y_pow = cap_y.clone().pow(j as u32);
            for (_x_idx, row) in qs.iter_mut().enumerate() {
                for (_y_idx, func) in row.iter_mut().enumerate() {
                    func[func_idx] *= &cap_x_pow;
                    func[func_idx] *= &cap_y_pow;
                }
            }
        }
    }

    let lattice = qs.concat();

    if debug {
        println!("L:");
        for v in &lattice {
            println!(" {:?}", v);
        }
    }

    let (reduced, _min_idx) = cryptlib::lll(&lattice);

    if debug {
        println!("reduced:");
        for v in &reduced {
            println!(" {:?}", v);
        }
    }

    let mut normalized_h = reduced[0].clone();
    for i in 0..(dx + k) {
        for j in 0..(dy + k) {
            normalized_h[j * (dx + k) + i] /= cap_x.clone().pow(i as u32);
            normalized_h[j * (dx + k) + i] /= cap_y.clone().pow(j as u32);
        }
    }

    let reversed_f: Vec<Vec<Integer>> = (0..=dy)
        .map(|j| (0..=dx).map(|i| f[i][j].clone()).collect())
        .collect();
    let reversed_h: Vec<Vec<Integer>> = (0..=(dy + k))
        .map(|j| {
            (0..=(dx + k))
                .map(|i| {
                    reduced[0][i * (dx + k + 1) + j].clone()
                        / cap_x.clone().pow(i as u32)
                        / cap_y.clone().pow(j as u32)
                })
                .collect()
        })
        .collect();
    if debug {
        println!("h^T {:?}", reversed_h);
        println!("f^T {:?}", reversed_f);
    }
    let mut cap_q = cryptlib::resultant(&reversed_h, &reversed_f, &Integer::from(-1));
    let mut gcd = Integer::from(0);
    for coef in &cap_q {
        gcd = cryptlib::gcd(coef, &gcd);
    }
    if debug {
        println!("gcd: {}", gcd);
    }

    for (i, coef) in cap_q.iter_mut().enumerate() {
        *coef /= &gcd;
        *coef *= cap_x.clone().pow(i as u32);
    }

    if debug {
        println!("Q {:?}", cap_q);
    }

    let x_candidates = cryptlib::approximate_zero(&cap_q, &cap_x);

    for (i, coef) in cap_q.iter_mut().enumerate() {
        *coef /= cap_x.clone().pow(i as u32);
    }

    let search_range = 5;
    let mut x_val;
    for candidate in &x_candidates {
        for i in -search_range..=search_range {
            let x = Integer::from(candidate + i);
            let f_of_x = cryptlib::eval_poly(&x, &cap_q, &Integer::from(-1));
            if f_of_x == 0 {
                if debug {
                    println!("f({}) = 0", x);
                }
                x_val = x;

                let mut y_poly: Vec<Integer> = reversed_f
                    .iter()
                    .map(|term| cryptlib::eval_poly(&x_val, term, &Integer::from(-1)))
                    .collect();

                if debug {
                    println!("poly for y: {:?}", y_poly);
                }

                for (i, coef) in y_poly.iter_mut().enumerate() {
                    *coef *= cap_y.clone().pow(i as u32);
                }

                let y_candidates = cryptlib::approximate_zero(&y_poly, &cap_y);

                for (i, coef) in y_poly.iter_mut().enumerate() {
                    *coef /= cap_y.clone().pow(i as u32);
                }

                let y_val;
                for candidate in &y_candidates {
                    // println!("trying {}", candidate);
                    for i in -search_range..=search_range {
                        let y = Integer::from(candidate + i);
                        let f_of_y = cryptlib::eval_poly(&y, &y_poly, &Integer::from(-1));
                        if f_of_y == 0 {
                            if debug {
                                println!("f({}) = 0", y);
                            }
                            y_val = y;
                            return Some((x_val, y_val));
                        }
                    }
                }
            }
        }
    }

    return None;
}

fn print_basis_bv(basis: &Vec<Vec<Vec<Rational>>>, indent: i32) {
    let indent_str: String = (0..indent).map(|_| " ").collect();
    for vec in basis {
        print!("{}", indent_str);
        for (i, row) in vec.iter().enumerate() {
            print!("{:?}x{}, ", row, i);
        }
        println!("]");
    }
}

fn print_poly_bv(f: &Vec<Vec<Integer>>) {
    print!(": [");
    for (i, row) in f.iter().enumerate() {
        print!("{:?}x^{}, ", row, i);
    }
    println!("]");
}
