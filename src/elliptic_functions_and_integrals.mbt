///|
/// Incomplete elliptic integral of the first kind F(phi, m)
pub fn ellipkinc(phi : Double, m : Double) -> Double {
  if m >= 1.0 {
    abort("Parameter m must be less than 1")
  }
  let f = t => @math.pow(1.0 - m * @math.sin(t) * @math.sin(t), -0.5)
  let (result, _, _) = @itgt.adap_quad_gk(
    f, @itgt.kronrod_r21, 0.0, phi, 1.0e-10, 1.0e-10, 100,
  ).unwrap()
  result
}

///|
/// Complete elliptic integral of the first kind K(m)
pub fn ellipk(m : Double) -> Double {
  ellipkinc(@math.PI / 2.0, m)
}

///|
/// Incomplete elliptic integral of the second kind F(phi, m)
pub fn ellipeinc(phi : Double, m : Double) -> Double {
  if m >= 1.0 {
    abort("Parameter m must be less than 1")
  }
  let f = t => (1.0 - m * @math.sin(t) * @math.sin(t)).sqrt()
  let (result, _, _) = @itgt.adap_quad_gk(
    f, @itgt.kronrod_r21, 0.0, phi, 1.0e-10, 1.0e-10, 100,
  ).unwrap()
  result
}

///|
/// Complete elliptic integral of the second kind E(m)
pub fn ellipe(m : Double) -> Double {
  ellipeinc(@math.PI / 2.0, m)
}

///|
/// Jacobian elliptic functions sn, cn, dn and amplitude am
/// Similar to SciPy's ellipj function
/// 
/// Parameters:
/// - u: argument
/// - m: parameter (0 <= m < 1)
/// 
/// Returns: (sn, cn, dn, phi) where phi is the amplitude
pub fn ellipj(u : Double, m : Double) -> (Double, Double, Double, Double) {
  if u == 0.0 {
    return (0.0, 1.0, 1.0, 0.0)
  }
  if m == 0.0 {
    let sn = @math.sin(u)
    let cn = @math.cos(u)
    let dn = 1.0
    let phi = u
    return (sn, cn, dn, phi)
  }
  let k = ellipk(m)
  let u_reduced = u % (4.0 * k)
  let mut u_work = u_reduced
  let mut sign_sn = 1.0
  let mut sign_cn = 1.0
  if u_work > 2.0 * k {
    u_work = u_work - 2.0 * k
    sign_sn = -sign_sn
    sign_cn = -sign_cn
  }
  if u_work > k {
    u_work = 2.0 * k - u_work
    sign_cn = -sign_cn
  }
  let ellipkinc_m = ellipkinc(_, m)
  let df = t => {
    let s = @math.sin(t)
    1.0 / (1.0 - m * s * s).sqrt()
  }
  let phi = match
    func_inv(u_work, ellipkinc_m, df, u, max_iterations=200, tol=1.0e-12) {
    None => abort("Failed to converge")
    Some(val) => val
  }
  let sin_phi = @math.sin(phi)
  let cos_phi = @math.cos(phi)
  let sn = sign_sn * sin_phi
  let cn = sign_cn * cos_phi
  let dn = (1.0 - m * sin_phi * sin_phi).sqrt()
  (sn, cn, dn, phi)
}

///|
test "ellipj basic cases" {
  let (sn, cn, dn, phi) = ellipj(0.0, 0.5)
  inspect(sn, content="0")
  inspect(cn, content="1")
  inspect(dn, content="1")
  inspect(phi, content="0")
  let (sn2, cn2, dn2, _) = ellipj(@math.PI / 4.0, 0.0)
  let expected_sn = @math.sin(@math.PI / 4.0)
  let expected_cn = @math.cos(@math.PI / 4.0)
  inspect((sn2 - expected_sn).abs() < 1.0e-10, content="true")
  inspect((cn2 - expected_cn).abs() < 1.0e-10, content="true")
  inspect(dn2, content="1")
}

///|
test "ellipj numerical values" {
  let (sn, cn, dn, phi) = ellipj(1.0, 0.5)

  // sn^2 + cn^2 = 1
  let identity1 = sn * sn + cn * cn
  inspect((identity1 - 1.0).abs() < 1.0e-10, content="true")

  // dn^2 + m*sn^2 = 1
  let identity2 = dn * dn + 0.5 * sn * sn
  inspect((identity2 - 1.0).abs() < 1.0e-10, content="true")
  inspect(
    (sn, cn, dn, phi),
    content="(0.8030018248956439, 0.5959765676721407, 0.8231610016315963, 0.9323150798838539)",
  )
}
