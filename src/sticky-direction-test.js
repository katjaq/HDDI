const fs = require('fs');
const linalgCode = fs.readFileSync('linalg.js').toString();
const rndCode = fs.readFileSync('random.js').toString();
const simCode = fs.readFileSync('sim.js').toString();
eval(linalgCode);
eval(rndCode);
eval(simCode);

var HDDI = function HDDI() {};

props=Object.keys(HDDILinAlg);
for(ind in props) {
    prop = props[ind];
    HDDI.prototype[prop] = HDDILinAlg[prop];
}

props=Object.keys(HDDIRandom);
for(ind in props) {
    prop = props[ind];
    HDDI.prototype[prop] = HDDIRandom[prop];
}

props=Object.keys(HDDISim);
for(ind in props) {
    prop = props[ind];
    HDDI.prototype[prop] = HDDISim[prop];
}
const h = new HDDI();

let mdir = {x:1, y:0, z:0};
let anisotropy = 0;
let g = 0.9;
let wd = 0.999;
let a, dot, i, j, r, v;
for(j=0;j<100;j++) {
    g = j/100;
    for(i=0;i<100000;i++) {
        r = h.randomDirection();
        v = {x:1.5, y:1.2, z:0.5};
        v = {
            x: v.x*r.x,
            y: v.y*r.y,
            z: v.z*r.z
        };
        dot = h.dot(mdir, v);
        if(dot<0) {
            v = h.scale(v, -1);
        }
        mdir = {
            x: wd * mdir.x + (1-wd)*v.x,
            y: wd * mdir.y + (1-wd)*v.y,
            z: wd * mdir.z + (1-wd)*v.z
        }
        a = h.norm(h.sub(v, h.scale(mdir, Math.abs(h.dot(v,mdir))/h.dot(mdir,mdir))));
        anisotropy = wd*anisotropy + (1-wd)*a;
    }
    console.log(mdir.x,mdir.y,mdir.z,1-anisotropy);
}