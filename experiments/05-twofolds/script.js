(function () {
    "use strict";

    console.log("Sphere test");

    let hddi = new HDDI();

    const params = {
        w: 0.9, // stiffness parameter
        step: 0.5, // step size
        minFibLength: 10,
        ns: 1e+7, // number of streamlines to throw
        dir: [      // array containing the diffusion directions
            {x: 0.817324, y: -0.49673, z: -0.29196},
            {x: 0.465087, y: -0.03533, z: 0.88456},
            {x: 0.820439, y: -0.31517, z: 0.477018},
            {x: -0.80334, y: 0.593293, z: -0.05141},
            {x: -0.15636, y: 0.788990, z: -0.59418},
            {x: -0.11253, y: -0.34483, z: -0.93189}
        ]
    };
    hddi.params = params;

    const wdir = 'experiments/05-twofolds/results/';
    const exec = require('child_process').execSync;

    // generate an ellipsoid with two folds
    let vol = [];
    const dim = [80, 80, 80];
    let i, j, k;
    let f, r, theta;
    for(i=0;i<dim[0];i++) {
        for(j=0;j<dim[1];j++) {
            for(k=0;k<dim[2];k++) {
                r = Math.sqrt(Math.pow(i-dim[0]/2,2) + Math.pow(j-dim[1]/2,2) + Math.pow(k-dim[2]/2,2));
                theta = Math.atan2(k-dim[2]/2, i-dim[0]/2);
                f = 30
                   *(1+Math.cos(2*theta)/10)
                   *(1-0.5*Math.exp(-Math.pow((theta-Math.PI/2)/0.2, 2)))
                   *(1-0.5*Math.exp(-Math.pow((theta+Math.PI/2)/0.2, 2)));
                if(r<f) {
                    hddi.setValue(vol, dim, i, j, k, 1);
                } else {
                    hddi.setValue(vol, dim, i, j, k, 0);
                }
            }
        }
    }

    // identify surface
    const ident = hddi.identifyVoxels(vol, dim);
    saveNiftiData(ident, dim, wdir + 'mask.nii.gz');
    
    // generate streamlines
    hddi.initialise(dim);
    hddi.genStreamlines(ident, dim);
    saveNiftiData(hddi.rho, dim, wdir + 'rho.nii.gz');

    // create a b0 image as max(dir)
    const b0 = hddi.computeB0(dim);

/*
    // compute and save first diffusion directions
    const frst = hddi.firstDirection(res.dir, dim);
    Promise.all([
        saveNiftiData(frst[0],dim, wdir + 'red.nii.gz'),
        saveNiftiData(frst[1],dim, wdir + 'green.nii.gz'),
        saveNiftiData(frst[2],dim, wdir + 'blue.nii.gz'),
    ]).then(() => {
        fsl.merge([
            wdir + 'rgb.nii.gz',
            wdir + 'red.nii.gz',
            wdir + 'green.nii.gz',
            wdir + 'blue.nii.gz'
        ]);
        exec(['rm', wdir + 'red.nii.gz', wdir + 'green.nii.gz', wdir + 'blue.nii.gz'].join(' '));
    });
*/

    // save results
    const dir = [];
    const dirs = [];
    for(i=0; i<hddi.params.dir.length; i++) {
        dir.push(wdir + 'dir{num}.nii.gz'.replace('{num}', i));
        dirs.push(wdir + 'dir{num}s.nii.gz'.replace('{num}', i));
    }
    let arr = [saveNiftiData(b0,dim, wdir + 'b0.nii.gz')];
    for(i=0; i<hddi.params.dir.length; i++) {
        arr.push(saveNiftiData(hddi.dir[i],dim,dir[i]));
    }
    Promise
    .all(arr)
    .then(() => {
        // gaussian smooth
        for(i=0; i<hddi.params.dir.length; i++) {
            fsl.maths( [dir[i], '-kernel gauss 3 -fmean', dirs[i]] );
        }
        // merge
        fsl.merge(['-t', wdir + 'dwi.nii.gz', wdir + 'b0.nii.gz', ...dir]);
        // remove intermediate files
        exec(['rm', ...dir, ...dirs, wdir + 'b0.nii.gz'].join(' '));

        // do tractography
        // transform bvecs: negative x-axis
        let vecs = hddi.params.dir.map(o=>{o.x = -o.x; return o});
        // save bvals & bvecs
        fs.writeFileSync(wdir + 'bvals.txt', [0, ...[...Array(vecs.length)].map(o=>1000)].join(' '));
        fs.writeFileSync(wdir + 'bvecs.txt', [0, ...vecs.map(o=>o.x), '\n'].join(' '));
        fs.appendFileSync(wdir + 'bvecs.txt', [0, ...vecs.map(o=>o.y), '\n'].join(' '));
        fs.appendFileSync(wdir + 'bvecs.txt', [0, ...vecs.map(o=>o.z), '\n'].join(' '));

        mrtrix.mrconvert([wdir + 'dwi.nii.gz','-fslgrad', wdir + 'bvecs.txt', wdir + 'bvals.txt', wdir + 'dwi.mif', '-force']);
        mrtrix.dwi2tensor([wdir + 'dwi.mif', wdir + 'dt.mif', '-force']);
        mrtrix.tensor2metric(['-fa', wdir + 'fa.nii.gz', '-adc', wdir + 'adc.nii.gz', '-num', 1 ,'-vector', wdir + 'v1.nii.gz', wdir + 'dt.mif ', '-force']);
        mrtrix.dwi2tensor([wdir + 'dwi.mif', wdir + 'dt.mif', '-force']);
        mrtrix.tckgen([wdir + 'dwi.mif', wdir + 'streamlines.50k-det.tck', '-algorithm', 'Tensor_Det', '-seed_image', wdir + 'mask.nii.gz', '-select', 50000, '-force']);
    });
} ());