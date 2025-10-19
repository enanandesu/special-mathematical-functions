///|
pub(all) enum Monotonicity {
  Inc
  Dec
}

///|
/// Calculates the inverse of a function
pub fn func_inv(
  target : Double,
  f : (Double) -> Double,
  d : (Double) -> Double,
  initial_guess : Double,
  max_iterations? : Int = 100,
  tol? : Double = 1.0e-12,
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

///|
pub fn bisection_guess(
  target : Double,
  f : (Double) -> Double,
  min_range : Double,
  max_range : Double,
  max_iterations? : Int = 200,
  tol? : Double = 1.0e-14,
  monotonicity : Monotonicity,
) -> Double? {
  let mut left = min_range
  let mut right = max_range
  let f_left = f(left)
  let f_right = f(right)
  let target_in_range = match monotonicity {
    Inc => f_left <= target && target <= f_right
    Dec => f_right <= target && target <= f_left
  }
  if not(target_in_range) {
    return None
  }
  for iter in 0..<max_iterations {
    let mid = (left + right) / 2.0
    let f_mid = f(mid)
    if (f_mid - target).abs() < tol {
      return Some(mid)
    }
    match monotonicity {
      Inc => if f_mid < target { left = mid } else { right = mid }
      Dec => if f_mid > target { left = mid } else { right = mid }
    }
    if (right - left).abs() < tol {
      return Some((left + right) / 2.0)
    }
  }
  None
}
