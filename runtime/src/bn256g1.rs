/// G1 point arithmetic library
use rstd::prelude::*;
use primitives::{U256};

/// Galois field point
type gfp = Vec<U256>;

pub fn new_gfp(x: U256) -> Vec<U256>{
    if x >= U256::from(0) {
        let out = vec!{[x,x,x,x]};
        return out;
    } else {
        let out = vec!{[-x,-x,-x,-x]}
        return gfp_neg(out, out);
    }
}

pub fn set_gfp()





