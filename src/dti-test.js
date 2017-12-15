const fs = require('fs');
const linalgCode = fs.readFileSync('linalg.js').toString();
const dtiCode = fs.readFileSync('dti.js').toString();
eval(linalgCode);
eval(dtiCode);
var HDDI = function HDDI() {};
props=Object.keys(HDDILinAlg);
for(ind in props) {
    prop = props[ind];
    HDDI.prototype[prop] = HDDILinAlg[prop];
}
props=Object.keys(HDDIDTI);
for(ind in props) {
    prop = props[ind];
    HDDI.prototype[prop] = HDDIDTI[prop];
}
const h = new HDDI();

// gradient table
const G = [
    {x: 0.817324, y: -0.49673, z: -0.29196},
    {x: 0.465087, y: -0.03533, z: 0.88456},
    {x: 0.820439, y: -0.31517, z: 0.477018},
    {x: -0.80334, y: 0.593293, z: -0.05141},
    {x: -0.15636, y: 0.788990, z: -0.59418},
    {x: -0.11253, y: -0.34483, z: -0.93189}
];
//const signal = [312, 206, 129, 174, 233, 154, 187]; // 1st value is b0
const signal = [0, 206, 129, 174, 233, 154, 187]; // 1st value is b0
signal[0] = Math.max(...signal);
console.log(signal);

console.log("\nDTI using algebra");
console.time('svd1');
const iH = h.invHfromG(G);
const dif = h.diffusionMatrix6(signal, iH);
const eva = h.eigenvalues(dif);
console.log("eigenvalues:", eva);
const eve = h.eigenvectors(eva, dif);
console.log("eigenvectors:", eve);
const mea = h.md(eva);
console.log("MD:", mea);
const fra = h.fa(eva);
console.log("FA:", fra);
console.timeEnd('svd1');

console.log("\nDTI using SVD");
console.time('svd2');
//console.log("svd");
console.log(h.svd([
    [dif[0],dif[3],dif[4]],
    [dif[3],dif[1],dif[5]],
    [dif[4],dif[5],dif[2]]
]));
const {U,S,V} = h.svd([
    [dif[0],dif[3],dif[4]],
    [dif[3],dif[1],dif[5]],
    [dif[4],dif[5],dif[2]]
]);
console.log("eigenvalues:", S);
console.log("eigenvectors:", V);
console.log("MD:",h.md({l1:S[0],l2:S[1],l3:S[2]}));
console.log("FA:",h.fa({l1:S[0],l2:S[1],l3:S[2]}));
console.timeEnd('svd2');

const vec = {x:1000,y:0,z:0};
let i, j;
let dot = [];
/*
console.log("\nDTI on a volume");
let dwi = [];
const dim = [1,1,1];
for(i=0;i<6;i++) {
    dot[i] = Math.abs(h.dot(vec,G[i]));
}
console.log(dot);
for(i=0;i<6;i++) {
    dwi[i]=[];
    for(j=0;j<dim[0]*dim[1]*dim[2];j++) {
        dwi[i][j]=dot[i];
    }
}
console.log(JSON.stringify(h.dti(dwi, dim, G),2," "));
*/

console.log("\nTest");
dot = [];
for(i=0;i<6;i++) {
    let v = h.scale(G[i], h.dot(vec,G[i]));
    dot[i]=[v.x,v.y,v.z];
}
console.log(dot);
console.log(h.svd(dot));

