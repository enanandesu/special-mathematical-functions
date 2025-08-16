///|
pub fn gamma(x : Double) -> Double {
  @kmath.gamma(x)
}

///|
pub fn gammaln(x : Double) -> Double {
  @math.ln(gamma(x).abs())
}

///|
pub fn gammasgn(x : Double) -> Double {
  if gamma(x) > 0 {
    1.0
  } else {
    -1.0
  }
}

///|
pub fn gammainc(a : Double, x : Double) -> Double {
  if a <= 0 {
    abort("a must be positive!")
  }
  if x < 0 {
    abort("x must be in [0,+inf)")
  }
  if x == 0 {
    return 0.0
  }
  let f = t => @math.pow(t, a - 1.0) * @math.exp(-t)
  let (result, _, _) = @itgt.adap_quad_gk(
    f, @itgt.kronrod_r31, 1.0e-50, x, 1.0e-12, 1.0e-12, 200,
  ).unwrap()
  result / gamma(a)
}

///|
pub fn gammaincinv(a : Double, y : Double) -> Double {
  if a <= 0 {
    abort("a must be positive!")
  }
  if y < 0 || y > 1 {
    abort("y must be in [0,1]")
  }
  let f = gammainc(a, _)
  let ans = bisection_guess(y, f, 0.0, 50.0, "I")
  let result = match ans {
    None => abort("No solution!")
    Some(val) => val
  }
  result
}

///|
pub fn beta(a : Double, b : Double) -> Double {
  gamma(a) * gamma(b) / gamma(a + b)
}

///|
pub fn betaln(a : Double, b : Double) -> Double {
  @math.ln(beta(a, b).abs())
}

///|
pub fn betainc(a : Double, b : Double, x : Double) -> Double {
  if x < 0 || x > 1 {
    abort("x must be in [0,1]")
  }
  if x == 0 || x == 1 {
    return x
  }
  let f = t => @math.pow(t, a - 1.0) * @math.pow(1.0 - t, b - 1.0)
  let (result, _, _) = @itgt.adap_quad_gk(
    f, @itgt.kronrod_r31, 1.0e-50, x, 1.0e-12, 1.0e-12, 200,
  ).unwrap()
  result / beta(a, b)
}

///|
pub fn betaincc(a : Double, b : Double, x : Double) -> Double {
  1 - betainc(a, b, x)
}

///|
pub fn betaincinv(a : Double, b : Double, y : Double) -> Double {
  if y < 0 || y > 1 {
    abort("y must be in [0,1]")
  }
  if y == 0 || y == 1 {
    return y
  }
  let f = betainc(a, b, _)
  let ans = bisection_guess(y, f, 0.0, 1.0, "I")
  let result = match ans {
    None => abort("No solution!")
    Some(val) => val
  }
  result
}

///|
pub fn betainccinv(a : Double, b : Double, y : Double) -> Double {
  betaincinv(a, b, 1.0 - y)
}

///|
pub fn rgamma(x : Double) -> Double {
  1.0 / gamma(x)
}

///|
pub fn digamma(x : Double) -> Double {
  let (d, _) = @deriv.deriv_central(gamma, x, 1.0e-3)
  d / gamma(x)
}

///|
pub fn polygamma(n : Int, x : Double) -> Double {
  if n < 0 {
    abort("n must be non-negative")
  }
  let sgn = if n % 2 == 1 { 1.0 } else { -1.0 }
  let nfact = factorial(n).to_double()
  let mut k = 0.0
  let mut tmp = @math.pow(x + k, (n + 1).to_double())
  let mut sum = 0.0
  while nfact / tmp > 1.0e-10 {
    sum = sum + nfact / tmp
    k = k + 1
    tmp = @math.pow(x + k, (n + 1).to_double())
  }
  sum * sgn
}

///|
fn factorial(n : Int) -> Int {
  if n == 0 || n == 1 {
    1
  } else {
    let mut result = 1
    for i = 2; i <= n; i = i + 1 {
      result = result * i
    }
    result
  }
}

///|
pub fn poch(z : Double, m : Double) -> Double {
  gamma(z + m) / gamma(z)
}

///|
test {
  inspect(gamma(1.73), content="0.9146653711862724")
  inspect(gammaln(1.73), content="-0.08919699514701278")
  inspect(gammasgn(1.73), content="1")
  inspect(gammainc(0.71, 1.73), content="0.8946863363707583")
  inspect(gammaincinv(0.71, 0.26), content="0.14269189008369043")
  inspect(beta(1.73, 2.89), content="0.12035540077280849")
  inspect(betaln(1.73, 2.89), content="-2.1173062402068927")
  inspect(betainc(1.73, 2.89, 0.7), content="0.926323292643955")
  inspect(betaincc(1.73, 2.89, 0.7), content="0.07367670735604503")
  inspect(betaincinv(1.73, 2.89, 0.7), content="0.48237044609452767")
  inspect(betainccinv(1.73, 2.89, 0.7), content="0.24058042039824823")
  inspect(rgamma(1.73), content="1.0932960091219515")
  inspect(digamma(1.73), content="0.23207695930131053")
  inspect(poch(1.73, 2.89), content="15.049901208898861")
  inspect(polygamma(3, 1.73), content="0.834568977525701")
}

///|
test {
  inspect(gammainc(0.11, 2.0), content="0.9936669667532537")
  inspect(gammaincinv(0.71, 0.88), content="1.6137426489557072")
  inspect(gammaincinv(0.37, 0.8), content="0.5910570267289916")
  inspect(gammaincinv(0.11, 0.23), content="9.643432719030898e-7")
  inspect(gammaincinv(0.9, 0.5), content="0.5967430489562275")
}

///|
test {
  inspect(betaincinv(0.5, 0.5, 0.7), content="0.7938926261462598")
  inspect(betaincinv(0.1, 0.1, 0.9), content="0.9999999113521802")
  inspect(betaincinv(0.9, 0.9, 0.1), content="0.0859772604256932")
  inspect(betaincinv(0.1, 0.9, 0.1), content="1.1806022826021945e-10")
}
