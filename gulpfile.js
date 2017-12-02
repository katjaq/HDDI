const gulp = require('gulp');
const concat = require('gulp-concat');
const print = require('gulp-print');

const sourceFiles = ['extend.js', 'mii.js', 'dti.js', 'sim.js', 'app.js'];
const experimentFiles = ['01-ellipsoid/script.js', '02-rectangle/script.js']
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
