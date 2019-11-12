//// Pairing library from zokrates verifier contract

use rstd::prelude::*;
use primitives::{U256};



/// Add two twist points
/// 
/// # Arguments
/// 
/// * pt1xx Coefficient 1 of x on point 1
/// * pt1xy Coefficient 2 of x on point 1
/// * pt1yx Coefficient 1 of y on point 1
/// * pt1yy Coefficient 2 of y on point 1
/// * pt2xx Coefficient 1 of x on point 2
/// * pt2xy Coefficient 2 of x on point 2
/// * pt2yx Coefficient 1 of y on point 2
/// * pt2yy Coefficient 2 of y on point 2
/// * returns (pt3xx, pt3xy, pt3yx, pt3yy)
pub fn ECTwistAdd(pt1xx: U256, pt1xy: U256, pt1yx: U256, pt1yy: U256, pt2xx: U256, pt2xy:U256, pt2yx:U256, pt2yy: U256) -> (U256, U256, U256, U256) {
    // Get constants
    let (FIELD_MODULUS, TWISTBX, TWISTBY) = get_constants();

    if (pt1xx == U256::from(0) && pt1xy == U256::from(0) && pt1yx == U256::from(0) && pt1yy == U256::from(0)) {
        if(!(pt2xx == U256::from(0) && pt2xy == U256::from(0) && pt2yx == U256::from(0) && pt2yy == U256::from(0))) {
            if(!_isOnCurve(pt2xx, pt2xy, pt2yx, pt2yy, FIELD_MODULUS, TWISTBX, TWISTBY)) {
                // replace solidity's assert with returning error value
                return(U256::from(1), U256::from(1), U256::from(1), U256::from(1));    
            }
            return (pt2xx, pt2xy, pt2yx, pt2yy);
        }
    } else if (pt2xx== U256::from(0) && pt2xy == U256::from(0) && pt2yx == U256::from(0) && pt2yy == U256::from(0)) {
        if(!_isOnCurve(pt1xx, pt1xy, pt1yx, pt1yy, FIELD_MODULUS, TWISTBX, TWISTBY)) {
            // replace solidity's assert with returning error value
            return(U256::from(1), U256::from(1), U256::from(1), U256::from(1));
        }
        return (pt1xx, pt1xy, pt1yx, pt1yy);
    } 

    if(!_isOnCurve(pt1xx, pt1xy, pt1yx, pt1yy, FIELD_MODULUS, TWISTBX, TWISTBY)) {
        // replace solidity's assert with returning error value
        return(U256::from(1), U256::from(1), U256::from(1), U256::from(1));
    }

    if(!_isOnCurve(pt2xx, pt2xy, pt2yx, pt2yy, FIELD_MODULUS, TWISTBX, TWISTBY)) {
        // replace solidity's assert with returning error value
        return(U256::from(1), U256::from(1), U256::from(1), U256::from(1));
    }

    let mut pt3: Vec<U256> = _ECTwistAddJacobian(pt1xx, pt1xy, pt1yx, pt1yy, U256::from(1), U256::from(0), pt2xx, pt2xy, pt2yx, pt2yy, U256::from(1), U256::from(0), FIELD_MODULUS);
        

    return _fromJacobian(pt3[PTXX], pt3[PTXY], pt3[PTYX], pt3[PTYY], pt3[PTZX], pt3[PTZY], FIELD_MODULUS);
}

/// Multiply a twist point by a scalar
/// 
/// # Arguments
/// 
/// * s Scalar to multiply by
/// * pt1xx Coefficient 1 of x on point 1
/// * pt1xy Coefficient 2 of x on point 1
/// * pt1yx Coefficient 1 of y on point 1
/// * pt1yy Coefficient 2 of y on point 1
/// * returns (pt2xx, pt2xy, pt2yx, pt2yy)
pub fn ECTwistMul(s: U256, mut pt1xx: U256, pt1xy: U256, mut pt1yx: U256, pt1yy: U256) -> (U256, U256, U256, U256) {
    // Get constants
    let (FIELD_MODULUS, TWISTBX, TWISTBY) = get_constants();
    let mut pt1zx: U256 = U256::from(1);

    if (pt1xx == U256::from(0) && pt1xy == U256::from(0) && pt1yx == U256::from(0) && pt1yy == U256::from(0)) {
        pt1xx = U256::from(1);
        pt1yx = U256::from(1);
        pt1zx = U256::from(0);
    } else {
        if(!_isOnCurve(pt1xx, pt1xy, pt1yx, pt1yy, FIELD_MODULUS, TWISTBX, TWISTBY)) {
            // replace solidity's assert with returning error value
            return(U256::from(1), U256::from(1), U256::from(1), U256::from(1));
        }
    }

    let pt2 = _ECTwistMulJacobian(s, pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, U256::from(0), FIELD_MODULUS);

    return _fromJacobian(pt2[PTXX], pt2[PTXY], pt2[PTYX], pt2[PTYY], pt2[PTZX], pt2[PTZY], FIELD_MODULUS);
}
// helper functions


const PTXX: usize = 0;
const PTXY: usize = 1;
const PTYX: usize = 2;
const PTYY: usize = 3;
const PTZX: usize = 4;
const PTZY: usize = 5;


/// get all constants for calculation
fn get_constants() -> (U256, U256, U256) {
    let FIELD_MODULUS: U256 = U256::from_dec_str("30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47").expect("p to be the field modulus");
    let TWISTBX: U256 = U256::from_dec_str("2b149d40ceb8aaae81be18991be06ac3b5b4c5e559dbefa33267e6dc24a138e5").expect("p to be TWISTBX");
    let TWISTBY: U256 = U256::from_dec_str("9713b03af0fed4cd2cafadeed8fdf4a74fa084e52d1852e4a2bd0685c315d2").expect("p to be TWISTBY");

    return (FIELD_MODULUS, TWISTBX, TWISTBY);
}

/// Compute addition in finite system
fn addmod(a: U256, b: U256, n: U256) -> U256 {
    return a + b % n;
}

/// Compute substraction in finite system
fn submod(a: U256, b: U256, n: U256) -> U256 {
    return addmod(a, n - b, n);
}

/// Compute multiplication in finite system
fn mulmod(a: U256, b: U256, n: U256) -> U256 {
    return a * b % n;
}


fn _fq2_mul(xx: U256, xy: U256, yx: U256, yy: U256, field_modulus: U256) -> (U256, U256) {
    return (
        submod(mulmod(xx, yx, field_modulus), mulmod(xx, yx, field_modulus), field_modulus),
        addmod(mulmod(xx, yy, field_modulus), mulmod(xy, yx, field_modulus), field_modulus)
    );
}

fn _fq2_muc(xx: U256, xy: U256, c: u32, field_modulus: U256) -> (U256, U256) {
    return (
        mulmod(xx, U256::from(c), field_modulus),
        mulmod(xy, U256::from(c), field_modulus)
    );
}

fn _fq2_add(xx:U256, xy:U256, yx: U256, yy: U256, field_modulus: U256) -> (U256, U256) {
    return (
        addmod(xx, yx, field_modulus),
        addmod(xy, yy, field_modulus)
    );
}

fn _fq2_sub(xx:U256, xy:U256, yx: U256, yy: U256, field_modulus: U256) -> (U256, U256) {
    return (
        submod(xx, yx, field_modulus),
        submod(xy, yy, field_modulus)
    );
}

/// Extended Eulidean algorithm(rucursive)
fn egcd(a: U256, b: U256) -> (U256, U256, U256) {
    assert!(a < b);
    if a == U256::zero() {
        return (b, U256::zero(), U256::one());
    }
    else {
        let (g, x, y) = egcd(b % a, a);
        return (g, y - (b / a) * x, x);
    }
}


fn _modInv(a: U256, n: U256) -> Option<U256> {
    let (g, x, _) = egcd(a, n);
    if g != U256::one() {
        return None;
    }
    else {
        return Some(x % n);
    }
}

fn _fq2_inv(x: U256, y: U256, field_modulus: U256) -> (U256, U256) {
    let inv: U256 = _modInv(addmod(mulmod(y, y, field_modulus), mulmod(x, x, field_modulus), field_modulus), field_modulus).unwrap();
    return (
        mulmod(x, inv, field_modulus),
        field_modulus - mulmod(y, inv, field_modulus)
    );
}

fn _isOnCurve(xx: U256, xy: U256, yx: U256, yy: U256, field_modulus: U256, TWISTBX: U256, TWISTBY: U256) -> bool {
    let (yyx, yyy) = _fq2_mul(yx, yy, yx, yy, field_modulus);
    let (xxxx, xxxy) = _fq2_mul(xx, xy, xx, xy, field_modulus);
    let (xxxx_1, xxxy_1) = _fq2_mul(xxxx, xxxy, xx, yy, field_modulus);
    let (yyx_1, yyy_1) = _fq2_sub(yyx, yyy, xxxx, xxxy, field_modulus);
    let (yyx_2, yyy_2) = _fq2_sub(yyx_1.clone(), yyy_1.clone(), TWISTBX, TWISTBY, field_modulus);
    return yyx_2 == U256::from(0) && yyy_2 == U256::from(0);
}




fn _fromJacobian(pt1xx: U256, pt1xy: U256, pt1yx: U256, pt1yy: U256, pt1zx: U256, pt1zy: U256, field_modulus: U256) -> (U256, U256, U256, U256) {
    let (invzx, invzy) = _fq2_inv(pt1zx, pt1zy, field_modulus);
    let (pt2xx, pt2xy) = _fq2_mul(pt1xx, pt1xy, invzx, invzy, field_modulus);
    let (pt2yx, pt2yy) = _fq2_mul(pt1yx, pt1yy, invzx, invzy, field_modulus);
    return (pt2xx, pt2xy, pt2yx, pt2yy);
}


fn _ECTwistAddJacobian(
    mut pt1xx: U256, mut pt1xy: U256, mut pt1yx: U256, mut pt1yy: U256, mut pt1zx: U256, mut pt1zy: U256,
    mut pt2xx: U256, mut pt2xy: U256, mut pt2yx: U256, mut pt2yy: U256, mut pt2zx: U256, mut pt2zy: U256,
    field_modulus: U256
    ) -> Vec<U256> {
        let mut pt3: Vec<U256> = vec!{U256::from(0),U256::from(0), U256::from(0), U256::from(0), U256::from(0), U256::from(0)};
        if (pt1zx == U256::from(0) && pt1zy == U256::from(0)) || (pt2zx == U256::from(0) && pt2zy == U256::from(0)) {
            pt3[PTXX] = pt1xx;
            pt3[PTXY] = pt1xy;
            pt3[PTYX] = pt1yx;
            pt3[PTYY] = pt1yy;
            pt3[PTZX] = pt1zx;
            pt3[PTZY] = pt1zy;
            return pt3;
        }

        let U1 = _fq2_mul(pt2yx, pt2yy, pt1zx, pt1zy, field_modulus); // U1 = y2 * z1
        pt2yx = U1.0;
        pt2yx = U1.1;
        let U2 = _fq2_mul(pt1yx, pt1yy, pt2zx, pt2zy, field_modulus); // U2 = y1 * z2
        pt3[PTYX] = U2.0;
        pt3[PTYY] = U2.1;
        let V1 = _fq2_mul(pt2xx, pt2xy, pt1zx, pt1zy, field_modulus); // V1 = x2 * z1
        pt2xx = V1.0;
        pt2xy = V1.1;
        let V2 = _fq2_mul(pt1xx, pt1xy, pt2zx, pt2zy, field_modulus); // V2 = x1 * z2
        pt3[PTZX] = V2.0;
        pt3[PTZY] = V2.1;

        if pt2xx == pt3[PTZX] && pt2xy == pt3[PTZY] {
            if pt2yx == pt3[PTYX] && pt2yy == pt3[PTYY] {
                let result  = _ECTwistDoubleJacobian(
                    pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, pt1zy,
                    field_modulus);
                pt3[PTXX] = result.0;
                pt3[PTXY] = result.1;
                pt3[PTYX] = result.2;
                pt3[PTYY] = result.3;
                pt3[PTZX] = result.4;
                pt3[PTZY] = result.5;
                return pt3;
            } 
            pt3[PTXX] = U256::from(1);
            pt3[PTXY] = U256::from(0);
            pt3[PTYX] = U256::from(1);
            pt3[PTYY] = U256::from(0);
            pt3[PTZX] = U256::from(0);
            pt3[PTZY] = U256::from(0);
            return pt3;
        }

        let W = _fq2_mul(pt1zx, pt1zy, pt2zx, pt2zy, field_modulus);     // W = z1 * z2
        pt2zx = W.0;
        pt2zy = W.1;
        let U = _fq2_sub(pt2yx, pt2yy, pt3[PTYX], pt3[PTYY], field_modulus); // U = U1 - U2
        pt1xx = U.0;
        pt1xy = U.1;
        let V = _fq2_sub(pt2xx, pt2xy, pt3[PTZX], pt3[PTZY], field_modulus); // V = V1 - V2
        pt1yx = V.0;
        pt1yy = V.1;
        let V_squared = _fq2_mul(pt1yx, pt1yy, pt1yx, pt1yy, field_modulus);     // V_squared = V * V
        pt1zx = V_squared.0;
        pt1zy = V_squared.1;
        let V_squared_times_V2 = _fq2_mul(pt1zx, pt1zy, pt3[PTZX], pt3[PTZY], field_modulus); // V_squared_times_V2 = V_squared * V2
        pt2yx = V_squared_times_V2.0;
        pt2yy = V_squared_times_V2.1;
        let V_cubed  = _fq2_mul(pt1zx, pt1zy, pt1yx, pt1yy, field_modulus);     // V_cubed = V * V_squared
        pt1zx = V_cubed.0;
        pt1zy = V_cubed.1;
        let newz = _fq2_mul(pt1zx, pt1zy, pt2zx, pt2zy, field_modulus);     // newz = V_cubed * W
        pt3[PTZX] = newz.0;
        pt3[PTZY] = newz.1;     
        let UxU     = _fq2_mul(pt1xx, pt1xy, pt1xx, pt1xy, field_modulus);     // U * U
        pt2xx = UxU.0;     
        pt2xy = UxU.1;
        let UxUxW     = _fq2_mul(pt2xx, pt2xy, pt2zx, pt2zy, field_modulus);     // U * U * W
        pt2xx = UxUxW.0;     
        pt2xy = UxUxW.1;
        let UxUxW_minus_V_cubed     = _fq2_sub(pt2xx, pt2xy, pt1zx, pt1zy, field_modulus);     // U * U * W - V_cubed
        pt2xx = UxUxW_minus_V_cubed.0;
        pt2xy = UxUxW_minus_V_cubed.1;
        let twoxV_squared_times_V2     = _fq2_muc(pt2yx, pt2yy, 2, field_modulus);                    // 2 * V_squared_times_V2
        pt2zx = twoxV_squared_times_V2.0;
        pt2zy = twoxV_squared_times_V2.1;
        let  A  = _fq2_sub(pt2xx, pt2xy, pt2zx, pt2zy,  field_modulus);     // A = U * U * W - V_cubed - 2 * V_squared_times_V2
        pt2xx = A.0;
        pt2xy = A.1;
        let newx = _fq2_mul(pt1yx, pt1yy, pt2xx, pt2xy, field_modulus);     // newx = V * A
        pt3[PTXX] = newx.0;
        pt3[PTXY] = newx.1;
        let V_squared_times_V2_minus_A     = _fq2_sub(pt2yx, pt2yy, pt2xx, pt2xy, field_modulus);     // V_squared_times_V2 - A
        pt1yx = V_squared_times_V2_minus_A.0;
        pt1yy = V_squared_times_V2_minus_A.1;
        let  UxV_squared_times_V2_minus_A    = _fq2_mul(pt1xx, pt1xy, pt1yx, pt1yy, field_modulus);     // U * (V_squared_times_V2 - A)
        pt1yx = UxV_squared_times_V2_minus_A.0;
        pt1yy = UxV_squared_times_V2_minus_A.1;
        let  V_cubedxU2 = _fq2_mul(pt1zx, pt1zy, pt3[PTYX], pt3[PTYY], field_modulus); // V_cubed * U2
        pt1xx = V_cubedxU2.0;     
        pt1xy = V_cubedxU2.1;
        let newy = _fq2_sub(pt1yx, pt1yy, pt1xx, pt1xy, field_modulus);     // newy = U * (V_squared_times_V2 - A) - V_cubed * U2
        pt3[PTYX] = newy.0;
        pt3[PTYY] = newy.1;
        return pt3;
}



fn _ECTwistDoubleJacobian(
    mut pt1xx: U256, mut pt1xy: U256, mut pt1yx: U256, mut pt1yy: U256, mut pt1zx: U256, mut pt1zy: U256,
    field_modulus: U256
    ) -> (U256, U256, U256, U256, U256, U256) {
        /// Rust does not have multiple variable assignment yet. fuck.
        let three_x = _fq2_muc(pt1xx, pt1xy, 3, field_modulus);            // 3 * x
        let mut pt2xx = three_x.0;  
        let mut pt2xy = three_x.1;
        let W  = _fq2_mul(pt2xx, pt2xy, pt1xx, pt1xy, field_modulus); // W = 3 * x * x
        pt2xx = W.0; 
        pt2xy = W.1;
        let S = _fq2_mul(pt1yx, pt1yy, pt1zx, pt1zy, field_modulus); // S = y * z
        pt1zx = S.0;
        pt1zy = S.1;
        let xy = _fq2_mul(pt1xx, pt1xy, pt1yx, pt1yy, field_modulus); // x * y
        let mut pt2yx = xy.0;
        let mut pt2yy = xy.1;
        let B = _fq2_mul(pt2yx, pt2yy, pt1zx, pt1zy, field_modulus); // B = x * y * S
        pt2yx = B.0;
        pt2yy = B.1;
        let W_squared = _fq2_mul(pt2xx, pt2xy, pt2xx, pt2xy, field_modulus); // W * W
        pt1xx = W_squared.0;
        pt1xy = W_squared.1;
        let eight_B = _fq2_muc(pt2yx, pt2yy, 8, field_modulus);            // 8 * B
        let mut pt2zx = eight_B.0;
        let mut pt2zy = eight_B.1;
        let H = _fq2_sub(pt1xx, pt1xy, pt2zx, pt2zy, field_modulus); // H = W * W - 8 * B
        pt1xx = H.0;
        pt1xy = H.1;
        let S_squared = _fq2_mul(pt1zx, pt1zy, pt1zx, pt1zy, field_modulus); // S_squared = S * S
        pt2zx = S_squared.0;
        pt2zy = S_squared.1;
        let four_B = _fq2_muc(pt2yx, pt2yy, 4, field_modulus);            // 4 * B
        pt2yx = four_B.0;
        pt2yy = four_B.1;
        let four_B_minus_H = _fq2_sub(pt2yx, pt2yy, pt1xx, pt1xy, field_modulus); // 4 * B - H
        pt2yx = four_B_minus_H.0;
        pt2yy = four_B_minus_H.1;
        let Wxfour_b_minus_H = _fq2_mul(pt2yx, pt2yy, pt2xx, pt2xy, field_modulus); // W * (4 * B - H)
        pt2yx = Wxfour_b_minus_H.0;
        pt2yy = Wxfour_b_minus_H.1;
        let eight_y = _fq2_muc(pt1yx, pt1yy, 8, field_modulus);            // 8 * y
        pt2xx = eight_y.0;
        pt2xy = eight_y.1;
        let eight_y_squared = _fq2_mul(pt2xx, pt2xy, pt1yx, pt1yy, field_modulus); // 8 * y * y
        pt2xx = eight_y_squared.0;
        pt2xy = eight_y_squared.1;
        let eight_y_squared_S_squared = _fq2_mul(pt2xx, pt2xy, pt2zx, pt2zy, field_modulus); // 8 * y * y * S_squared
        pt2xx = eight_y_squared_S_squared.0;
        pt2xy = eight_y_squared_S_squared.1;
        let newy = _fq2_sub(pt2yx, pt2yy, pt2xx, pt2xy, field_modulus); // newy = W * (4 * B - H) - 8 * y * y * S_squared
        pt2yx = newy.0;
        pt2yy = newy.1;
        let two_H = _fq2_muc(pt1xx, pt1xy, 2, field_modulus);            // 2 * H
        pt2xx = two_H.0;
        pt2xy = two_H.1;
        let newx = _fq2_mul(pt2xx, pt2xy, pt1zx, pt1zy, field_modulus); // newx = 2 * H * S
        pt2xx = newx.0;
        pt2xy = newx.1;
        let S_S_squared = _fq2_mul(pt1zx, pt1zy, pt2zx, pt2zy, field_modulus); // S * S_squared
        pt2zx = S_S_squared.0;
        pt2zy = S_S_squared.1;
        let newz = _fq2_muc(pt2zx, pt2zy, 8, field_modulus);            // newz = 8 * S * S_squared
        pt2zx = newz.0;
        pt2zy = newz.1;
    return (pt2xx, pt2xy, pt2yx, pt2yy, pt2zx, pt2zy);
}


fn _ECTwistMulJacobian(
    mut d: U256,
    mut pt1xx: U256, mut pt1xy: U256, mut pt1yx: U256, mut pt1yy: U256, mut pt1zx: U256, mut pt1zy: U256,
    field_modulus: U256
) -> Vec<U256> {
    let mut pt2: Vec<U256> = vec!{U256::from(0),U256::from(0), U256::from(0), U256::from(0), U256::from(0), U256::from(0)};
    while d != U256::from(0) {
        if((d & U256::from(1)) != U256::from(0)) {
            pt2 = _ECTwistAddJacobian(
                pt2[PTXX], pt2[PTXY], pt2[PTYX], pt2[PTYY], pt2[PTZX], pt2[PTZY],
                pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, pt1zy,
                field_modulus
                );
        }
        
        let epoch = _ECTwistDoubleJacobian(
            pt1xx, pt1xy, pt1yx, pt1yy, pt1zx, pt1zy,
            field_modulus
            );

        pt1xx = epoch.0;
        pt1xy = epoch.1;
        pt1yx = epoch.2;
        pt1yy = epoch.3;
        pt1zx = epoch.4;
        pt1zy = epoch.5;

        d = d / 2;
    }

    return pt2;
}

