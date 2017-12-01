const exec = require('child_process').execSync;
var exports = module.exports = {};

const merge = function(args) {
    const arr = ['fslmerge ', ...args];
    exec( arr.join(' '), (res) => {console.log(res);});
}

const maths = function(args) {
    const arr = ['fslmaths', ...args];
    exec( arr.join(' '), (res) => {console.log(res);});
}

exports.merge = merge;
exports.maths = maths;