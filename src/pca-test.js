const fs = require('fs');
const linalgCode = fs.readFileSync('linalg.js').toString();
const pcaCode = fs.readFileSync('pca.js').toString();
eval(linalgCode);
eval(pcaCode);
var HDDI = function HDDI() {};
props=Object.keys(HDDILinAlg);
for(ind in props) {
    prop = props[ind];
    HDDI.prototype[prop] = HDDILinAlg[prop];
}
props=Object.keys(HDDIPCA);
for(ind in props) {
    prop = props[ind];
    HDDI.prototype[prop] = HDDIPCA[prop];
}
const h = new HDDI();


/**
  * @desc mvnrnd([0 0 0],[1.5 0 0;0 1 0;0 0 0.5],10)
  */

const M = [
    {x: -1.612941, y: 0.872697, z: 0.103421},
    {x: 0.564237, y: 0.042549, z: 0.410649},
    {x: 1.629862, y: -0.456570, z: 0.369193},
    {x: 1.990129, y: -0.124546, z: -1.170756},
    {x: -1.278364, y: -0.421444, z: -0.477042},
    {x: -0.912089, y: 0.595086, z: -0.158703},
    {x: -0.649690, y: -2.229267, z: -1.167557},
    {x: 1.920755, y: 0.664606, z: 0.712530},
    {x: 3.290676, y: 0.191062, z: 0.023735},
    {x: -0.407322, y: 1.194171, z: -0.837170}
];

const pc1 = h.firstPC(M);
const pc2 = h.nextPC(M, pc1);

console.log(pc1);
console.log(pc2);

