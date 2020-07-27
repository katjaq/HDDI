(function () {
    "use strict";

    console.log("Two folds test, with different gradients, sticky fibres");

    let hddi = new HDDI();

    const params = {
        w: 0.9, // stiffness parameter
        step: 1, // step size
        minFibLength: 15,
        dmass: 0.999,
        ns: 1e+5, // number of streamlines to throw
        dir: [      // array containing the diffusion directions
            { x:-0.049, y:-0.996, z:-0.074 },
            { x:-0.996, y:0.043, z:0.080 },
            { x:0.078, y:-0.078, z:0.994 },
            { x:0.498, y:-0.556, z:-0.665 },
            { x:-0.524, y:-0.663, z:0.535 },
            { x:0.708, y:0.344, z:0.617 }
        ]
    };
    hddi.params = params;

    const wdir = 'experiments/10-twofolds-different-gradients-sticky/results/';
    const exec = require('child_process').execSync;
    const imgdir = 'experiments/images/';
    const imgprefix0 = '10-twofolds-different-gradients-sticky_sag_'
    const imgprefix1 = '10-twofolds-different-gradients-sticky_cor_'
    const imgprefix2 = '10-twofolds-different-gradients-sticky_axi_'
    const imgzoom = 20;

    try {
        fs.statSync(wdir);
    } catch(e) {
        fs.mkdirSync(wdir);
    }
    try {
        fs.statSync(imgdir);
    } catch(e) {
        fs.mkdirSync(imgdir);
    }

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
    let mask = ident.map((v) => (v === 1));
    saveNiftiData(mask, dim, wdir + 'mask.nii.gz');

    // generate streamlines
    hddi.initialise(dim);
    hddi.genStickyStreamlines(ident, dim);
    saveNiftiData(hddi.rho, dim, wdir + 'rho.nii.gz');
    saveNiftiData(hddi.params.anis, dim, wdir + 'anisotropy.nii.gz');

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
        mrtrix.tckgen([wdir + 'dwi.mif', wdir + 'streamlines.50k-det.tck', '-algorithm', 'Tensor_Det', '-seed_image', wdir + 'mask.nii.gz', '-mask', wdir + 'mask.nii.gz', '-select', 50000, '-force']);
    });
    const fafile = wdir + 'fa.nii.gz';
    const tckfile = wdir + 'streamlines.50k-det.tck';
    //size width, height in pixel units
    //fov field of view in mm
    const imgsize0 = dim[1] * imgzoom + ',' + dim[2] * imgzoom;
    const imgfov0 = Math.max(dim[1], dim[2]);
    const imgsize1 = dim[0] * imgzoom + ',' + dim[2] * imgzoom;
    const imgfov1 = Math.max(dim[0], dim[2]);
    const imgsize2 = dim[0] * imgzoom + ',' + dim[1] * imgzoom;
    const imgfov2 = Math.max(dim[0], dim[1]) + imgzoom;
    Promise
    .all(tckfile)
    .then(() => {
        console.log(fafile);
        console.log(imgdir);
        //-plane index --> Set the viewing plane, according to the mappping 0: sagittal; 1: coronal; 2: axial.
        //-tractography.thickness
        mrtrix.mrview([`-load ${fafile}`, '-noannotations', '-plane 0', `-size ${imgsize0}`, `-fov ${imgfov0}`, `-tractography.load ${tckfile}`, '-tractography.opacity 0.3', '-imagevisible 0', `-capture.folder ${imgdir}`, `-capture.prefix ${imgprefix0}`, '-capture.grab', '-force', '-exit']);
        mrtrix.mrview([`-load ${fafile}`, '-noannotations', '-plane 1', `-size ${imgsize1}`, `-fov ${imgfov1}`, `-tractography.load ${tckfile}`, '-tractography.opacity 0.3', '-imagevisible 0', `-capture.folder ${imgdir}`, `-capture.prefix ${imgprefix1}`, '-capture.grab', '-force', '-exit']);
        mrtrix.mrview([`-load ${fafile}`, '-noannotations', '-plane 2', `-size ${imgsize2}`, `-fov ${imgfov2}`, `-tractography.load ${tckfile}`, '-tractography.opacity 0.3', '-imagevisible 0', `-capture.folder ${imgdir}`, `-capture.prefix ${imgprefix2}`, '-capture.grab', '-force', '-exit']);
    });
} ());