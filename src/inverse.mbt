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
  monotonicity : String,
) -> Double? {
  if monotonicity != "I" && monotonicity != "D" {
    return None
  }
  let mut left = min_range
  let mut right = max_range
  let f_left = f(left)
  let f_right = f(right)
  let target_in_range = match monotonicity {
    "I" => f_left <= target && target <= f_right
    "D" => f_right <= target && target <= f_left
    _ => false
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
      "I" => if f_mid < target { left = mid } else { right = mid }
      "D" => if f_mid > target { left = mid } else { right = mid }
      _ => return None
    }
    if (right - left).abs() < tol {
      return Some((left + right) / 2.0)
    }
  }
  None
}
