/*
Hypothesis-driven Diffusion Imaging.
A script to generate fibers based on the anatomy.

This simulation tries to see how much of a real DWI-based tractography we are able to recover by applying anatomical constraints to randomly drifting particles within the mask.
These constraints could be implemented based on the hypotheses that
    fibres run perpendicular to gyral crowns,
    follow sulcal fundi,
    have a homogeneous density throughout the white matter
    and are sticky, forming bundles of similar orientation.

Randomly, a surface voxel will be chosen and start propagating a particle into a random direction until it hits a surface. This process is repeated a huge number of times.
Work in progress: Direction values will be stored throughout all iterations, while counting voxel visits and repelling fibers based on density...

2016, katja & roberto
*/

const progress = require('cli-progress');
const sys = require('util')
const fsl = require('./interfaces/fsl');
const mrtrix = require('./interfaces/mrtrix');

var HDDIApp = {
    debug: 0,
    boundaryValue: 2,        // border value
    boundary: null,     // border voxels
    vol: null,      // new volum array (storing direction value sum of all passed fibers)
    rho: null,    // fibre density volume
    dir: null     // fibre direction volume
}

var HDDI = function HDDI() {};
let ind, prop;

console.log('\nExtending HDDI from HDDIApp');
//keys=Object.keys(HDDIApp);for(p in keys){HDDIApp[keys[p]]);HDDI.prototype[keys[p]] = HDDIApp[keys[p]]}
props=Object.keys(HDDIApp);
for(ind in props) {
    prop = props[ind];
    if(typeof prop !== 'undefined') {
        HDDI.prototype[prop] = HDDIApp[prop];
    }
}

console.log('\nExtending HDDI from HDDISim');
props=Object.keys(HDDISim);
for(ind in props) {
    prop = props[ind];
    if(typeof prop !== 'undefined') {
        HDDI.prototype[prop] = HDDISim[prop];
    }
}

console.log('\nExtending HDDI from HDDISobel3');
props=Object.keys(HDDISobel3);
for(ind in props) {
    prop = props[ind];
    if(typeof prop !== 'undefined') {
        HDDI.prototype[prop] = HDDISobel3[prop];
    }
}

console.log('\nExtending HDDI from HDDIRandom');
props=Object.keys(HDDIRandom);
for(ind in props) {
    prop = props[ind];
    if(typeof prop !== 'undefined') {
        HDDI.prototype[prop] = HDDIRandom[prop];
    }
}

/*
    To generate gradient tables:
    http://www.emmanuelcaruyer.com/WebApp/q-space-sampling.php?nbPoints=6&nbShells=1&alpha=1
*/