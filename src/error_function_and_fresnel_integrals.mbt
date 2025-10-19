///|
/// Returns the error function of argument.
pub fn erf(x : Double) -> Double {
  let f = t => @math.exp(-t * t)
  let (result, _, _) = @itgt.adap_quad_gk(
    f, @itgt.kronrod_r21, 0.0, x, 1.0e-10, 1.0e-10, 100,
  ).unwrap()
  2.0 / @math.PI.sqrt() * result
}

///|
test "erf" {
  inspect(erf(0.0), content="0")
  inspect(erf(1.0), content="0.8427007929497149")
  inspect(erf(-1.0), content="-0.8427007929497149")
  inspect(erf(2.0), content="0.9953222650189526")
  inspect(erf(-2.0), content="-0.9953222650189526")
  inspect(erf(0.5), content="0.5204998778130465")
  inspect(erf(-0.5), content="-0.5204998778130465")
  inspect(erf(3.0), content="0.9999779095030016")
  inspect(erf(-3.0), content="-0.9999779095030016")
  inspect(erf(5.0), content="0.9999999999984626")
  inspect(erf(-5.0), content="-0.9999999999984626")
  inspect(erf(10.0), content="1")
  inspect(erf(-10.0), content="-1")
}

///|
/// Complementary error function, 1 - erf(x).
pub fn erfc(x : Double) -> Double {
  1 - erf(x)
}

///|
/// Scaled complementary error function, exp(x**2) * erfc(x).
pub fn erfcx(x : Double) -> Double {
  @math.exp(x * x) * erfc(x)
}

///|
/// Inverse of the error function.
pub fn erfinv(x : Double) -> Double {
  guard -1.0 <= x && x <= 1.0 else { abort("Range Error!") }
  let df = t => 2.0 / @math.PI.sqrt() * @math.exp(-t * t)
  let ans = func_inv(x, erf, df, @math.PI.sqrt() / 2.0 * x)
  let result = match ans {
    None => abort("No solution!")
    Some(val) => val
  }
  result
}

///|
test "erfinv" {
  inspect(erfinv(0.5), content="0.4769362762044699")
  inspect(erfinv(0.8), content="0.9061938024368235")
  inspect(erfinv(0.9999), content="2.7510639057122948")
  inspect(erfinv(-0.8), content="-0.9061938024368235")
  inspect(erfinv(-0.9999), content="-2.7510639057122948")
}

///|
/// Inverse of the complementary error function.
pub fn erfcinv(y : Double) -> Double {
  erfinv(1.0 - y)
}

///|
/// Dawson's integral.
pub fn dawsn(x : Double) -> Double {
  let f = t => @math.exp(t * t)
  let (result, _, _) = @itgt.adap_quad_gk(
    f, @itgt.kronrod_r21, 0.0, x, 1.0e-10, 1.0e-10, 100,
  ).unwrap()
  @math.exp(x * x) * result
}

///|
/// Fresnel integrals.
pub fn fresnel(x : Double) -> (Double, Double) {
  let fs = t => @math.sin(@math.PI * 0.5 * t * t)
  let fc = t => @math.cos(@math.PI * 0.5 * t * t)
  let (s, _, _) = @itgt.adap_quad_gk(
    fs, @itgt.kronrod_r21, 0.0, x, 1.0e-10, 1.0e-10, 100,
  ).unwrap()
  let (c, _, _) = @itgt.adap_quad_gk(
    fc, @itgt.kronrod_r21, 0.0, x, 1.0e-10, 1.0e-10, 100,
  ).unwrap()
  (s, c)
}
