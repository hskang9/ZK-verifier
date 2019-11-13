use primitives::{U256};

use crate::bn256g1::*;
use crate::bn256g2::*;

pub struct G1Point {
    X: U256,
    Y: U256
}

pub struct G2Point {
    X: [U256; 2],
    Y: [U256; 2]
}

pub fn P1() -> G1Point {
    return G1Point{
        X: U256::from(1), 
        Y: U256::from(2)
    };
}

pub fn P2() -> G2Point {
    return G2Point{
        X: [U256::from_dec_str("11559732032986387107991004021392285783925812861821192530917403151452391805634").expect("p is xx"),
             U256::from_dec_str("10857046999023057135944570762232829481370756359578518086990519993285655852781").expect("p is xy")],
        Y: [U256::from_dec_str("4082367875863433681332203403145435568316851327593401208105741076214120093531").expect("p is yx"),
             U256::from_dec_str("8495653923123431417604973247489272438418190587263600148770280649306958101930").expect("p is yy")]
    };
}

pub fn negate(p: G1Point) -> G1Point {
    let q: U256 = U256::from_dec_str("21888242871839275222246405745257275088696311157297823662689037894645226208583").expect("p is q");
    if p.X == U256::from(0) && p.Y == U256::from(0) {
        return G1Point(U256::from(0), U256::from(0));
    }
    return G1Point{p.X, q - (p.Y % q)};
}

/// EC point addition
/// curvePoint implements the elliptic curve y²=x³+3. Points are kept in Jacobian
/// form and t=z² when valid. G₁ is the set of points of this curve on GF(p).
/// Check geth's https://github.com/ethereum/go-ethereum/blob/22e3bbbf0aa711fcab3876ac1a34b8cdff8a0ccd/core/vm/contracts.go#L297 runBn256Add.
/// for Adding the point https://github.com/ethereum/go-ethereum/blob/22e3bbbf0aa711fcab3876ac1a34b8cdff8a0ccd/crypto/bn256/cloudflare/bn256.go#L72
/// gfp https://github.com/ethereum/go-ethereum/blob/22e3bbbf0aa711fcab3876ac1a34b8cdff8a0ccd/crypto/bn256/cloudflare/gfp.go#L10
pub fn addition(p1: G1Point, p2: G1Point) {
    let input: Vec<U256> = vec!{[U256::from(0), U256::from(0), U256::from(0), U256::from(0)]};
    input[0] = p1.X;
    input[1] = p1.Y;
    input[2] = p2.X;
    input[3] = p2.Y;
    // Run bn256 Add between input's p1 and p2 
    //bn256g1::runBn256Add()
}