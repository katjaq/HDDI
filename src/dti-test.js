const fs = require('fs');
const code = fs.readFileSync('dti.js').toString();
eval(code);

const dir = [
    {x: 0.817324, y: -0.49673, z: -0.29196},
    {x: 0.465087, y: -0.03533, z: 0.88456},
    {x: 0.820439, y: -0.31517, z: 0.477018},
    {x: -0.80334, y: 0.593293, z: -0.05141},
    {x: -0.15636, y: 0.788990, z: -0.59418},
    {x: -0.11253, y: -0.34483, z: -0.93189}
];
const S = [312, 206, 129, 174, 233, 154, 187]; // 1st value is b0

const iH = invHfromG(dir);
const dif = (diffusionMatrix6(S, iH))._data;
console.log("dif:", dif);
const eva = eigenvalues(dif);
console.log("eigenvalues:", eva);
const eve = eigenvectors(eva, dif);
console.log("eigenvectors:", eve);
const mea = md(eva);
console.log("MD:", mea);
const fra = fa(eva);
console.log("FA:", fra);
