const exec = require('child_process').execSync;
var exports = module.exports = {};

const mrconvert = function(args) {
    const arr = ['mrconvert ', ...args];
    exec( arr.join(' '), (res) => {console.log(res);});
}

const dwidenoise = function(args) {
    const arr = ['dwidenoise', ...args];
    exec( arr.join(' '), (res) => {console.log(res);});
}

const dwi2tensor = function(args) {
    const arr = ['dwi2tensor', ...args];
    exec( arr.join(' '), (res) => {console.log(res);});
}

const tensor2metric = function(args) {
    const arr = ['tensor2metric', ...args];
    exec( arr.join(' '), (res) => {console.log(res);});
}

const tckgen = function(args) {
    const arr = ['tckgen', ...args];
    exec( arr.join(' '), (res) => {console.log(res);});
}

exports.mrconvert = mrconvert;
exports.dwidenoise = dwidenoise;
exports.dwi2tensor = dwi2tensor;
exports.tensor2metric = tensor2metric;
exports.tckgen = tckgen;