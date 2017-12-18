const gulp = require('gulp');
const concat = require('gulp-concat');
const print = require('gulp-print');

const sourceFiles = ['tract.js', 'nii.js', 'linalg.js', 'dti.js', 'sim.js', 'sobel3.js', 'random.js', 'app.js'];
const experimentFiles = [
/*
    '01-ellipsoid/script.js',
    '02-rectangle/script.js',
    '03-sphere/script.js',
*/
    '04-homogeneous-sphere/script.js'
/*
    '05-twofolds/script.js'
    '06-twofolds-different-gradients/script.js',
    '07-ferret-p4/script.js',
    '08-ferret-p4-longer-fibres/script.js',
    '09-sphere-sticky/script.js',
    '10-twofolds-different-gradients-sticky/script.js'
*/
];

const src = 'src/';
const exp = 'experiments/';
const dest = './';

gulp.task('default', ['pack-hddi']);

gulp.task('pack-hddi', function () {
    return gulp.src([
        ...(sourceFiles.map(o => src + o)),
        ...(experimentFiles.map(o => exp + o))
        ])
        .pipe(print())
        .pipe(concat('hddi.js'))
        .pipe(gulp.dest(dest));
});
