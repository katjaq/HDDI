(function () {
    "use strict";

    console.log("Two folds test, with different gradients");

    let hddi = new HDDI();

    const params = {
        w: 0.95, // stiffness parameter
        step: 0.5, // step size
        minFibLength: 80,
        ns: 1e+5, // number of streamlines to throw
        dir: [      // array containing the diffusion directions
            {x:-0.049, y:-0.996, z:-0.074},
            {x:-0.996, y:0.043, z:0.080},
            {x:0.078, y:-0.078, z:0.994},
            {x:0.498, y:-0.556, z:-0.665},
            {x:-0.524, y:-0.663, z:0.535},
            {x:0.708, y:0.344, z:0.617}
        ]
    };
    hddi.params = params;

    const wdir = 'experiments/08-ferret-p4-longer-fibres/results/';
    const exec = require('child_process').execSync;
    const imgdir = 'experiments/images/';
    const imgprefix0 = '08-ferret-p4-longer-fibres_sag_'
    const imgprefix1 = '08-ferret-p4-longer-fibres_cor_'
    const imgprefix2 = '08-ferret-p4-longer-fibres_axi_'
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

    // load nifti1 mask
    loadNiftiData('experiments/08-ferret-p4-longer-fibres/F08_P4.t2_to_ref.sel.nii.gz')
    .then((res) => {
        let {vol, dim} = res;

        // identify surface
        const ident = hddi.identifyVoxels(vol, dim);
        let mask = ident.map((v) => (v === 1));
        saveNiftiData(mask, dim, wdir + 'mask.nii.gz');
    
        // generate streamlines
        hddi.initialise(dim);
        hddi.genStreamlines(ident, dim);
        saveNiftiData(hddi.rho, dim, wdir + 'rho.nii.gz');

        // create a b0 image as max(dir)
        const b0 = hddi.computeB0(dim);

        // save results
        const dir = [];
        const dirs = [];
        let i;
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
    })
    .catch((err) => console.log(err));
} ());