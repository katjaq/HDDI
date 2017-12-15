const fs = require('fs');
const code = fs.readFileSync('linalg.js').toString();
eval(code);

const SVD = require('./svd.js').svd;

/*
Example 1, from http://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm
Results should be:
U=[[0.82, -0.58, 0,0],
[0.58, 0.82, 0,0],
[0,0,1,0],
[0,0,0,1]]

S=[5.47, 0.37]

V=[[0.40, -0.91], [0.91, 0.40]]
*/

let A = [[2,4], [1,3], [0,0], [0,0]];

console.log("emscripten transcompiled SVD");
console.log(SVD(A));

console.log("manually translated from C code");
console.log(HDDILinAlg.svd(A));

/*
// Example 2, from http://stitchpanorama.sourceforge.net/Python/svd.py
*/
let B = [[22.,10., 2.,  3., 7.],
      [14., 7.,10.,  0., 8.],
      [-1.,13.,-1.,-11., 3.],
      [-3.,-2.,13., -2., 4.],
      [ 9., 8., 1., -2., 4.],
      [ 9., 1.,-7.,  5.,-1.],
      [ 2.,-6., 6.,  5., 1.],
      [ 4., 5., 0., -2., 2.]];

console.log("emscripten transcompiled SVD");
console.log(SVD(B));

console.log("manually translated from C code");
console.log(HDDILinAlg.svd(B));

console.log("\nTest matrix inversion");
inv = HDDILinAlg.svdinv([
    [0,1,2],
    [3,4,5],
    [6,7,9]
]);
console.log(inv);

console.log("\nTest multiplication and transpose");
let x=[[1,2],[3,4]];
let y=[[1,1,1,1],[2,2,2,2]];
let z=[[1,2,3,4],[5,6,7,8]];
console.log("mult",HDDILinAlg.mult(x,y));
console.log("transpose",HDDILinAlg.transpose(z));

console.log("\nTest diag");
const v = [1,2,3,4];
console.log("diag",HDDILinAlg.diag(v));
