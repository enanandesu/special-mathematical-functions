///|
/// Calculates the inverse of a function
pub fn func_inv(
  target : Double,
  f : (Double) -> Double,
  d : (Double) -> Double,
  initial_guess : Double,
  max_iterations~ : Int = 100,
  tol~ : Double = 1.0e-12,
) -> Double? {
  let mut x = initial_guess
  for iter in 0..<max_iterations {
    let fx = f(x)
    let dfx = d(x)
    if dfx.abs() < 1.0e-15 {
      return None
    }
    let new_x = x - (fx - target) / dfx
    if (new_x - x).abs() < tol {
      return Some(new_x)
    }
    x = new_x
  }
  None
}
