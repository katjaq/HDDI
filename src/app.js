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
    nifti1: null, // mask or T1 in format .nii
    nifti2: null, // identify surface voxels (set to bv) and inside volume (set 1 if > 0) and outside (0),
    nifti3: null, // contains value 1 for every visited voxel
    nifti4: null, // stores direction values from one visit
    bv: 2,        // border value
    bvox: [],     // border voxels
    vol: [],      // new volum array (storing direction value sum of all passed fibers)
    vvC: [],      // voxelVisitCount: array (blocksize) containing integer: how many times has voxel been passed by a fiber
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
        console.log("["+prop+"]");
        HDDI.prototype[prop] = HDDIApp[prop];
    }
}

console.log('\nExtending HDDI from HDDISim');
props=Object.keys(HDDISim);
for(ind in props) {
    prop = props[ind];
    if(typeof prop !== 'undefined') {
        console.log("["+prop+"]");
        HDDI.prototype[prop] = HDDISim[prop];
    }
}

console.log('\nExtending HDDI from HDDISobel3');
props=Object.keys(HDDISobel3);
for(ind in props) {
    prop = props[ind];
    if(typeof prop !== 'undefined') {
        console.log("["+prop+"]");
        HDDI.prototype[prop] = HDDISobel3[prop];
    }
}
//HDDI.prototype = extend(HDDI.prototype, HDDISim.prototype);
//HDDI.prototype = extend(HDDI.prototype, HDDISobel3.prototype);
console.log(HDDI.prototype);

/*
    To generate gradient tables:
    http://www.emmanuelcaruyer.com/WebApp/q-space-sampling.php?nbPoints=6&nbShells=1&alpha=1
*/