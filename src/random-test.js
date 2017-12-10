const fs = require('fs');
const code = fs.readFileSync('random.js').toString();
eval(code);

let hist = [];
const nbins = 30;
const nsamples = 10000;
let i, x;
let max, min;

console.log("nrandom, "+nsamples+" samples.");
min = -4;
max = 4;
for(i=0;i<nbins;i++) {
    hist[i] = 0;
}
for(i=0;i<nsamples;i++) {
    x = HDDIRandom.nrandom();
    if(x >= min && x < max) {
        hist[parseInt(nbins*(x-min)/(max-min))] += 1;
    }
}
for(i=0;i<nbins;i++) {
    console.log('|'+(new Array(parseInt(hist[i]/25))).join('*'));
}

console.log("randomDirection, Marsaglia's algorithm, X coordinate, "+nsamples+" samples.");
min = -1;
max = 1;
for(i=0;i<nbins;i++) {
    hist[i] = 0;
}
for(i=0;i<nsamples;i++) {
    x = HDDIRandom.randomDirection();
    hist[parseInt(nbins*(x.x-min)/(max-min))] += 1;
}
for(i=0;i<nbins;i++) {
    console.log('|'+(new Array(parseInt(hist[i]/25))).join('*'));
}
