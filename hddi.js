const fs = require('fs');
const Struct = require('Struct');
const zlib = require('zlib');

const NiiHdr = new Struct()
        .word32Sle('sizeof_hdr')        // Size of the header. Must be 348 (bytes)
        .chars('data_type', 10)         // Not used; compatibility with analyze.
        .chars('db_name', 18)           // Not used; compatibility with analyze.
        .word32Sle('extents')           // Not used; compatibility with analyze.
        .word16Sle('session_error')     // Not used; compatibility with analyze.
        .word8('regular')               // Not used; compatibility with analyze.
        .word8('dim_info')              // Encoding directions (phase, frequency, slice).
        .array('dim', 8, 'word16Sle')   // Data array dimensions.
        .floatle('intent_p1')           // 1st intent parameter.
        .floatle('intent_p2')           // 2nd intent parameter.
        .floatle('intent_p3')           // 3rd intent parameter.
        .word16Sle('intent_code')       // nifti intent.
        .word16Sle('datatype')          // Data type.
        .word16Sle('bitpix')            // Number of bits per voxel.
        .word16Sle('slice_start')       // First slice index.
        .array('pixdim', 8, 'floatle')  // Grid spacings (unit per dimension).
        .floatle('vox_offset')          // Offset into a .nii file.
        .floatle('scl_slope')           // Data scaling, slope.
        .floatle('scl_inter')           // Data scaling, offset.
        .word16Sle('slice_end')         // Last slice index.
        .word8('slice_code')            // Slice timing order.
        .word8('xyzt_units')            // Units of pixdim[1..4].
        .floatle('cal_max')             // Maximum display intensity.
        .floatle('cal_min')             // Minimum display intensity.
        .floatle('slice_duration')      // Time for one slice.
        .floatle('toffset')             // Time axis shift.
        .word32Sle('glmax')             // Not used; compatibility with analyze.
        .word32Sle('glmin')             // Not used; compatibility with analyze.
        .chars('descrip', 80)           // Any text.
        .chars('aux_file', 24)          // Auxiliary filename.
        .word16Sle('qform_code')        // Use the quaternion fields.
        .word16Sle('sform_code')        // Use of the affine fields.
        .floatle('quatern_b')           // Quaternion b parameter.
        .floatle('quatern_c')           // Quaternion c parameter.
        .floatle('quatern_d')           // Quaternion d parameter.
        .floatle('qoffset_x')           // Quaternion x shift.
        .floatle('qoffset_y')           // Quaternion y shift.
        .floatle('qoffset_z')           // Quaternion z shift.
        .array('srow_x', 4, 'floatle')  // 1st row affine transform
        .array('srow_y', 4, 'floatle')  // 2nd row affine transform.
        .array('srow_z', 4, 'floatle')  // 3rd row affine transform.
        .chars('intent_name', 16)       // Name or meaning of the data.
        .chars('magic', 4);             // Magic string.

function createNiftiHeader(dim, pixdim, dir) {
/*eslint-disable camelcase*/
    var datatype = 16; // float, 2; // uchar
    var i;
    var newHdr = {
        sizeof_hdr: 348,
        data_type: '',
        db_name: '',
        extents: 0,
        session_error: 0,
        regular: 0,
        dim_info: 0,
        dim: [3, dim[0], dim[1], dim[2], 1, 1, 1, 1],
        intent_p1: 0,
        intent_p2: 0,
        intent_p3: 0,
        intent_code: 0,
        datatype: datatype,
        bitpix: 32,
        slice_start: 0,
        pixdim: [-1, pixdim[0], pixdim[1], pixdim[2], 0, 1, 1, 1],
        vox_offset: 352,
        scl_slope: 1,
        scl_inter: 0,
        slice_end: 0,
        slice_code: 0,
        xyzt_units: 10,
        cal_max: 0,
        cal_min: 0,
        slice_duration: 0,
        toffset: 0,
        glmax: 0,
        glmin: 0,
        descrip: 'HDDI, 16 November 2017',
        aux_file: '',
        qform_code: 0,
        sform_code: 1,
        quatern_b: 0,
        quatern_c: 0,
        quatern_d: 0,
        qoffset_x: 0,
        qoffset_y: 0,
        qoffset_z: 0,
        srow_x: [dir[0][0], dir[0][1], dir[0][2], dir[0][3]],
        srow_y: [dir[1][0], dir[1][1], dir[1][2], dir[1][3]],
        srow_z: [dir[2][0], dir[2][1], dir[2][2], dir[2][3]],
        intent_name: '',
        magic: 'n+1'
    };
/*eslint-enable camelcase*/

    NiiHdr.allocate();
    let niihdr = NiiHdr.buffer();
    for(i in newHdr) {
        NiiHdr.fields[i] = newHdr[i];
    }

    return niihdr;
}
function saveNiftiData(vol, dim, path) {
    const sz = dim[0]*dim[1]*dim[2]+1;
    let mri = null;
    let hdr = createNiftiHeader(dim, [1,1,1], [[1,0,0,0],[0,1,0,0],[0,0,1,0]]);
    let img = new Float32Array(sz); //Buffer(sz);
    let buffer = new Buffer(img.buffer);
    let i;
    for(i = 0; i<sz; i++) {
        img[i] = vol[i];
    }
    mri = new Buffer(buffer.length + hdr.length);
    hdr.copy(mri);
    buffer.copy(mri, hdr.length);
    var pr = new Promise(function (resolve, reject) {
        zlib.gzip(mri, function (err, mrigz) {
            if(err) {
                reject(err);

                return;
            }
            fs.writeFileSync(path, mrigz);
            resolve("Nifti file saved");
        });
    });

    // fs.writeFileSync(path, hdr);
    // fs.appendFileSync(path,new Buffer(img.buffer));

    return pr;
}

/*
http://www.diffusion-imaging.com/2014/04/from-diffusion-weighted-images-to.html
*/

const math = require('mathjs');

function normVec(v) {
    let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
    return {
        x: v.x / lv,
        y: v.y / lv,
        z: v.z / lv
    };
}

/**
 * @func eigenvalues
 * @desc Returns the eigenvalues of a symmetric 3x3 matrix
 * @param array d A symmetric 3x3 diffusion matrix
 */
function eigenvalues(d) {
    const [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = d;
    const I1 = Dxx + Dyy + Dzz;
    const I2 = Dxx*Dyy + Dyy*Dzz + Dzz*Dxx - (Dxy*Dxy + Dxz*Dxz + Dyz*Dyz);
    const I3 = Dxx*Dyy*Dzz + 2*Dxy*Dxz*Dyz - (Dzz*Dxy*Dxy + Dyy*Dxz*Dxz + Dxx*Dyz*Dyz);
    const v = Math.pow(I1/3, 2) - I2/3;
    const s = Math.pow(I1/3, 3) - I1*I2/6 + I3/2;
    const phi = Math.acos(s/Math.pow(v, 3/2))/3;
    const l1 = Math.abs(I1/3 + 2*Math.sqrt(v)*Math.cos(phi));
    const l2 = Math.abs(I1/3 - 2*Math.sqrt(v)*Math.cos(Math.PI/3+phi));
    const l3 = Math.abs(I1/3 - 2*Math.sqrt(v)*Math.cos(Math.PI/3-phi));

    return {l1: l1, l2: l2, l3: l3};
}

/**
 * @func eigenvectors
 * @desc Returns the eigenvectors of a symmetric 3x3 matrix
 * @param array l A vector with 3 eigenvalues
 * @param array d A symmetric 3x3 diffusion matrix
 * @return object An object with three eigenvectors
 */
function eigenvectors(l, d) {
    const [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = d;
    //console.log("d:",Dxx, Dyy, Dzz, Dxy, Dxz, Dyz);
    const {l1, l2, l3} = l;
    //console.log("l:",l1,l2,l3);
    const A1 = Dxx - l1;
    const B1 = Dyy - l1;
    const C1 = Dzz - l1;
    const A2 = Dxx - l2;
    const B2 = Dyy - l2;
    const C2 = Dzz - l2;
    const A3 = Dxx - l3;
    const B3 = Dyy - l3;
    const C3 = Dzz - l3;
    const e1 = {
        x: (Dxy*Dyz - B1*Dxz) * (Dxz*Dyz - C1*Dxy),
        y: (Dxz*Dyz - C1*Dxy) * (Dxy*Dxz - A1*Dyz),
        z: (Dxy*Dxz - A1*Dyz) * (Dxy*Dyz - B1*Dxz)
    }
    const e2 = {
        x: (Dxy*Dyz - B2*Dxz) * (Dxz*Dyz - C2*Dxy),
        y: (Dxz*Dyz - C2*Dxy) * (Dxy*Dxz - A2*Dyz),
        z: (Dxy*Dxz - A2*Dyz) * (Dxy*Dyz - B2*Dxz)
    }
    const e3 = {
        x: (Dxy*Dyz - B3*Dxz) * (Dxz*Dyz - C3*Dxy),
        y: (Dxz*Dyz - C3*Dxy) * (Dxy*Dxz - A3*Dyz),
        z: (Dxy*Dxz - A3*Dyz) * (Dxy*Dyz - B3*Dxz)
    }

    return {e1: normVec(e1), e2: normVec(e2), e3: normVec(e3)};
}

function invHfromG(G) {
    const invH = math.inv(math.matrix([
        [G[0].x*G[0].x, G[0].y*G[0].y, G[0].z*G[0].z, 2*G[0].x*G[0].y, 2*G[0].x*G[0].z, 2*G[0].y*G[0].z],
        [G[1].x*G[1].x, G[1].y*G[1].y, G[1].z*G[1].z, 2*G[1].x*G[1].y, 2*G[1].x*G[1].z, 2*G[1].y*G[1].z],
        [G[2].x*G[2].x, G[2].y*G[2].y, G[2].z*G[2].z, 2*G[2].x*G[2].y, 2*G[2].x*G[2].z, 2*G[2].y*G[2].z],
        [G[3].x*G[3].x, G[3].y*G[3].y, G[3].z*G[3].z, 2*G[3].x*G[3].y, 2*G[3].x*G[3].z, 2*G[3].y*G[3].z],
        [G[4].x*G[4].x, G[4].y*G[4].y, G[4].z*G[4].z, 2*G[4].x*G[4].y, 2*G[4].x*G[4].z, 2*G[4].y*G[4].z],
        [G[5].x*G[5].x, G[5].y*G[5].y, G[5].z*G[5].z, 2*G[5].x*G[5].y, 2*G[5].x*G[5].z, 2*G[5].y*G[5].z]
    ]));

    return invH;
}

/**
 * @func diffusionMatrix6
 * @desc Computes a diffusion matrix based on signals measured in 6 different axes
 * @param array S An array of 7 values, first the signal without gradients, followed by signals in 6 directions
 * @param array G An array of 6 directions, each with x, y and z components
 * @return array A diffusion matrix expressed as a vector, d = [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]
 */
function diffusionMatrix6(S, invH) {
    let b = 1000; // s/mm2 (changing b only changes the eigenvalues, and thus MD, but not the eigenvectors or the FA)
    let S0 = S[0]; // larger than any Sk
    const Y = math.matrix([
        Math.log(S0/S[1])/b, 
        Math.log(S0/S[2])/b, 
        Math.log(S0/S[3])/b, 
        Math.log(S0/S[4])/b, 
        Math.log(S0/S[5])/b, 
        Math.log(S0/S[6])/b
    ]);

    return math.multiply(invH, Y)._data;
}

/**
 * @desc Mean diffusivity
 */
function md(evals) {
    return (evals.l1 + evals.l2 + evals.l3)/3;
}

/**
 * @desc Fractional anisotropy
 */
function fa(evals) {
    const {l1, l2, l3} = evals;
    const MD = (l1 + l2 + l3)/3;
    return Math.sqrt(
        (3*Math.pow(l1 - MD, 2)
        + Math.pow(l2 - MD, 2)
        + Math.pow(l3 - MD, 2))
        /(2*(l1*l1 + l2*l2 + l3*l3)));
}

var HDDISim = {
    getValue: function getValue(vol, dim, x, y, z) {
        return vol[z*dim[1]*dim[0] + y*dim[0] + x];
    },

    getValueN: function getValueN(vol, dim, x, y, z) {
        let i, val = [];
        for(i=0; i<this.params.dir.length; i++) {
            val.push(vol[i][z*dim[1]*dim[0] + y*dim[0] + x]);
        }

        return val;
    },

    setValue: function setValue(vol,dim, x, y, z, val) {
        vol[z*dim[1]*dim[0] + y*dim[0] + x] = val;
    },

    setValueN: function setValueN(vol,dim, x, y, z, val) {
        let i;
        for(i=0; i<this.params.dir.length; i++) {
            vol[i][z*dim[1]*dim[0] + y*dim[0] + x] = Math.abs(this.dot(this.params.dir[i], val));
        }
    },

    addValue: function addValue(vol,dim, x, y, z, val) {
        vol[z*dim[1]*dim[0] + y*dim[0] + x] += val;
    },

    addValueN: function addValueN(vol,dim, x, y, z, val) {
        let i;
        for(i=0; i<this.params.dir.length; i++) {
            vol[i][z*dim[1]*dim[0] + y*dim[0] + x] += Math.abs(this.dot(this.params.dir[i], val));
        }
    },

    normalise: function normalise(v, l) {
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: l * v.x / lv,
            y: l * v.y / lv,
            z: l * v.z / lv
        };
    },

    scale: function scale(v, l) {
        return {
            x: l * v.x,
            y: l * v.y,
            z: l * v.z
        };
    },

    /**
      * @func identifyVoxels
      * @desc Identify voxels as background (0), surface (boundaryValue) or core (1)
      * @param array vol An array containing voxels
      * @param object dim An array containing [nx, ny, nz], the dimensions of vol
      */
    identifyVoxels: function identifyVoxels(vol, dim) {
        console.log( '%cidentify mask voxels', 'color: green; ' );
        let idvol = [];
        let i, j, k, l;
        let x, y, z;
        let val;
        let con26 = [];

        // initialise 26-connected neighbourhood
        for( i = -1; i <= 1; ++i ) {
            for( j = -1; j<= 1; ++j ) {
                for( k= -1; k<= 1; ++k ) {
                    con26.push([i, j, k]);
                }
            }
        }

        // initialise volume to zero
        for( i = 0; i< dim[0]*dim[1]*dim[2]; i++) {
            idvol[i] = 0;
        }

        // mark "core" values
        for( x = 0; x < dim[0]; ++x ) {
            for( y = 0; y < dim[1]; ++y ) {
                for( z = 0; z < dim[2]; ++z) {
                    val = this.getValue( vol, dim, x , y, z );
                    if( val > 0 ) {
                        this.setValue( idvol, dim, x, y, z, 1 );
                    }
                    else {
                        this.setValue( idvol, dim, x, y, z, 0 );
                    }
                }
            }
        }

        // mark "surface" values
        for( x = 1; x < dim[0]-1; ++x ) {
            for( y = 1; y < dim[1]-1; ++y ) {
                for( z = 1; z < dim[2]-1; ++z) {
                    val = this.getValue( vol, dim, x , y, z );
                    if( val === 0 ) {
                        continue;
                    }
                    for( l=0; l<con26.length; l++ ) {
                        [i, j, k] = con26[l];
                        if( this.getValue( vol, dim, x+i, y+j, z+k ) === 0 ) {
                            this.setValue( idvol, dim, x, y, z, this.boundaryValue );
                            break;
                        }
                    }
                }
            }
        }

        // store coordinates of border voxels
        this.boundary = [];
        for( x = 0; x < dim[0]; ++x ) {
            for( y = 0; y < dim[1]; ++y ) {
                for( z = 0; z < dim[2]; ++z) {
                    if(this.getValue(idvol, dim, x, y, z) === this.boundaryValue ) {
                        this.boundary.push({ x: x, y: y, z: z });
                    }
                }
            }
        }

        if( this.debug > 2 ) {
            console.log( '%cborder voxels (index x, y, z): ', 'color: orange; ' );
            console.log( this.boundary );
        }

        return idvol;
    },

    /**
      * @func initialise
      * @desc Initialise the density (rho) and directions (dir) volumes
      */
    initialise: function initialise(dim) {
        console.log( '%cinitialise', 'color: green; ' );

        const sz = dim[0]*dim[1]*dim[2];
        let i, j;

        this.rho = new Float32Array(sz);
        this.dir = [];

        // initialise rho and dir volumes to 0
        for(i=0;i<sz;i++) {
            this.rho[i] = 0;
        }
        for(i=0;i<this.params.dir.length;i++) {
            this.dir[i] = new Float32Array(sz);
            for(j=0;j<sz;j++) {
                this.dir[i][j] = 0;
            }
        }
    },

    /**
      * @func genStreamlines
      * @desc Generates random fibres
      */
    genStreamlines: function genStreamlines(vol, dim) {
        console.log( '%cgenerate fibres', 'color: green; ' );

        let i;
        let count;
        let pos, rindex;
        let ix, iy, iz;
        let length, max, min;
        let v, v2;
        let bar = new progress.Bar({}, progress.Presets.shades_classic);

        bar.start(this.params.ns, 0);

        // generate fibres
        count = 0;
        while( count < this.params.ns) {
            let fib = [];

            if( count%1e+3 === 0 ) {
                bar.update(count);
            }

            //take random surface voxel to start fiber
            rindex = parseInt( Math.random() * ( this.boundary.length - 1 ) );
            pos = {
                x: this.boundary[rindex].x + 0.5,
                y: this.boundary[rindex].y + 0.5,
                z: this.boundary[rindex].z + 0.5
            };
            v = this.scale(this.randomDirection(), this.params.step);
            fib.push([pos, v]);
            length = 0;
            while( this.getValue( vol, dim, parseInt(pos.x), parseInt(pos.y), parseInt(pos.z) ) > 0 ) {
                // random direction
                v2 = this.scale(this.randomDirection(), this.params.step);
                //new direction = mix of old plus random
                v = this.normalise({
                    x: this.params.w * v.x + (1-this.params.w)* v2.x ,
                    y: this.params.w * v.y + (1-this.params.w)* v2.y,
                    z: this.params.w * v.z + (1-this.params.w)* v2.z
                }, this.params.step);

                // advance streamline
                pos = {
                    x: pos.x + v.x,
                    y: pos.y + v.y,
                    z: pos.z + v.z
                };
                fib.push([pos, v]);
                // compute length
                length += Math.sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
            }

            // filter out short fibres
            if(length < this.params.minFibLength ) {
                continue;
            }

            // compute min and max fibre length
            if(count === 0) {
                min = length;
                max = length;
            } else {
                if(length < min) {
                    min = length;
                }
                if(length > max) {
                    max = length;
                }
            }

            for(i = 0; i<fib.length; i += 1) {
                [pos, v] = fib[i];
                //set every visited voxel to 1
                ix = parseInt(pos.x);
                iy = parseInt(pos.y);
                iz = parseInt(pos.z);
                this.addValue( this.rho, dim, ix, iy, iz, 1 );
                this.addValueN( this.dir, dim, ix, iy, iz, v );
            }

            count += 1;
        }
        bar.update(this.params.ns);
        bar.stop();
        if( this.debug > 2 ) {
            console.log( '%cgenerate voxels finito', 'color: light-green; ' );
        }
        console.log("Min and Max fibre length:", min, max);
    },

    /**
      * @func genHomogeneousStreamlines
      * @desc Generates random fibres which avoid regions of high fibre density
      */
    genHomogeneousStreamlines: function genHomogeneousStreamlines(vol, dim) {
        console.log( '%cgenerate fibres', 'color: green; ' );

        let i;
        let count;
        let pos, rindex;
        let ix, iy, iz;
        let length, max, min;
        let v, v2, v3;
        let bar = new progress.Bar({}, progress.Presets.shades_classic);

        bar.start(this.params.ns, 0);

        // compute gradient of density image
        const dx = this.xder(this.rho, dim);
        const dy = this.yder(this.rho, dim);
        const dz = this.zder(this.rho, dim);
        for(i=0;i<dim[0]*dim[1]*dim[2];i++) {
            this.rho[i] = 0;
        }

        // generate fibres
        count = 0;
        while( count < this.params.ns) {
            let fib = [];

            if( count%1e+3 === 0 ) {
                bar.update(count);

            }

            // take random surface voxel to start fiber
            rindex = parseInt( Math.random() * ( this.boundary.length - 1 ) );
            pos = {
                x: this.boundary[rindex].x + 0.5,
                y: this.boundary[rindex].y + 0.5,
                z: this.boundary[rindex].z + 0.5
            };
            v = this.scale(this.randomDirection(), this.params.step);
            fib.push([pos, v]);
            length = 0;
            while( this.getValue( vol, dim, parseInt(pos.x), parseInt(pos.y), parseInt(pos.z) ) > 0 ) {
                // combine with random direction
                v2 = this.scale(this.randomDirection(), this.params.step);
                v = this.normalise({
                    x: this.params.w * v.x + (1-this.params.w)* v2.x ,
                    y: this.params.w * v.y + (1-this.params.w)* v2.y,
                    z: this.params.w * v.z + (1-this.params.w)* v2.z
                }, this.params.step);

                // combine with density-decreasing direction
                v3 = {
                    x: -this.getValue(dx, dim, parseInt(pos.x), parseInt(pos.y), parseInt(pos.z)),
                    y: -this.getValue(dy, dim, parseInt(pos.x), parseInt(pos.y), parseInt(pos.z)),
                    z: -this.getValue(dz, dim, parseInt(pos.x), parseInt(pos.y), parseInt(pos.z))
                };
                v = {
                    x: this.params.wh * v.x + (1-this.params.wh)* v3.x,
                    y: this.params.wh * v.y + (1-this.params.wh)* v3.y,
                    z: this.params.wh * v.z + (1-this.params.wh)* v3.z
                };

                // advance streamline
                pos = {
                    x: pos.x + v.x,
                    y: pos.y + v.y,
                    z: pos.z + v.z
                };
                fib.push([pos, v]);
                // compute length
                length += Math.sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
            }

            // filter out short fibres
            if(length < this.params.minFibLength ) {
                continue;
            }

            // fs.writeFileSync('fib.json', JSON.stringify(fib.map((o) => [o[0].x, o[0].y, o[0].z].join(' '))));

            // compute min and max fibre length
            if(count === 0) {
                min = length;
                max = length;
            } else {
                if(length < min) {
                    min = length;
                }
                if(length > max) {
                    max = length;
                }
            }

            for(i = 0; i<fib.length; i += 1) {
                [pos, v] = fib[i];
                //set every visited voxel to 1
                ix = parseInt(pos.x);
                iy = parseInt(pos.y);
                iz = parseInt(pos.z);
                this.addValue( this.rho, dim, ix, iy, iz, 1 );
                this.addValueN( this.dir, dim, ix, iy, iz, v );
            }

            count += 1;
        }
        bar.update(this.params.ns);
        bar.stop();
        if( this.debug > 2 ) {
            console.log( '%cgenerate voxels finito', 'color: light-green; ' );
        }
        console.log("Min and Max fibre length:", min, max);
            },
    
    /**
      * @func normaliseRho
      * @desc Normalise the density volume (rho) so that the maximum value is 1 if
      *       all streamlines pass through that coordinate. 
      */
    normaliseRho: function normaliseRho(dim) {
        let i, sz = dim[0]*dim[1]*dim[2];
        for(i=0; i<sz; i++) {
            this.rho /= this.params.ns;
        }
    },

    computeB0: function computeB0(dim) {
        let i, j, k;
        let b0 = new Float32Array(dim[0] * dim[1] *dim[2]);
        for(i=0;i<dim[0];i++) {
            for(j=0;j<dim[1];j++) {
                for(k=0;k<dim[2];k++) {
                    let sig = this.getValueN(this.dir, dim, i, j, k);
                    this.setValue(b0, dim, i, j, k, Math.max(...sig));
                }
            }
        }

        return b0;
    },

    firstDirection: function firstDirection(dir, dim) {
        let i, j, sz = dim[0] * dim[1] * dim[2];
        let max;
        let dif, evals, evecs;
        let first = [new Float32Array(sz), new Float32Array(sz), new Float32Array(sz)];
        const iH = invHfromG(this.params.dir);
        for(i=0; i<sz; i++) {
            let S = [0];
            for(j=0; j<this.params.dir.length; j++) {
                S[1+j] = dir[j][i];
            }
            // take S0 (signal without gradient) as the maximum value present
            max = Math.max(...S);
            if(max < 1) {
                continue;
            }
            S[0] = max;
            // compute diffusion matrix
            dif = diffusionMatrix6(S, iH);
            // get eigenvalues
            evals = eigenvalues(dif);
            // get eigenvectors
            evecs = eigenvectors(evals, dif);
            // store first eigenvector
            first[0][i] = evecs.e1.x;
            first[1][i] = evecs.e1.y;
            first[2][i] = evecs.e1.z;
        }

        return first;
    },

    dot: function dot( a, b ) {
        return( a.x * b.x + a.y * b.y + a.z * b.z )
    },

    getBvecsIn: function getBvecsIn() {
        console.log( '%cget bvecs in', 'color: green; ' );

        this.dir.push( { x:-0.26901271050312, y:0.426217878501591, z:-0.864094805089212 } );
        this.dir.push( { x:0.389308175473584, y:-0.464435056102534, z:-0.795882261217864 } );
        this.dir.push( { x:0.445081230516361, y:-0.1976420477346, z:0.873801848108661 } );
        // @todo read file from disc, genre { bvec[i].x, bvec[i].y, bvec[i].z };

        if( this.debug > 2 ) {
            console.log( '%cdirections dir (x, y, z): ', 'color: orange; ' );
            console.log( 'directions: ', this.dir );
        }
    }
}

var HDDISobel3 = {
    xder: function xder(vol, dim) {
        let i, j, k;
        let xres = new Float32Array(dim[0]*dim[1]*dim[2]);
        let yres = new Float32Array(dim[0]*dim[1]*dim[2]);
        let zres = new Float32Array(dim[0]*dim[1]*dim[2]);

        // h'(x)
        for(i=1; i<dim[0]-1; i++) {
            for(j=0; j<dim[1]; j++) {
                for(k=0; k<dim[2]; k++) {
                    xres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        - 1*vol[k*dim[1]*dim[0] + j*dim[0] + (i-1)]
                        + 1*vol[k*dim[1]*dim[0] + j*dim[0] + (i+1)];
                }
            }
        }

        // h(y)
        for(i=0; i<dim[0]; i++) {
            for(j=1; j<dim[1]-1; j++) {
                for(k=0; k<dim[2]; k++) {
                    yres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        1*xres[k*dim[1]*dim[0] + (j-1)*dim[0] + i]
                        + 2*xres[k*dim[1]*dim[0] + j*dim[0] + i]
                        + 1*xres[k*dim[1]*dim[0] + (j+1)*dim[0] + i];
                }
            }
        }

        // h(z)
        for(i=0; i<dim[0]; i++) {
            for(j=0; j<dim[1]; j++) {
                for(k=1; k<dim[2]-1; k++) {
                    zres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        1*yres[(k-1)*dim[1]*dim[0] + j*dim[0] + i]
                        + 2*yres[k*dim[1]*dim[0] + j*dim[0] + i]
                        + 1*yres[(k+1)*dim[1]*dim[0] + j*dim[0] + i];
                }
            }
        }

        return zres;
    },

    yder: function yder(vol, dim) {
        let i, j, k;
        let xres = new Float32Array(dim[0]*dim[1]*dim[2]);
        let yres = new Float32Array(dim[0]*dim[1]*dim[2]);
        let zres = new Float32Array(dim[0]*dim[1]*dim[2]);

        // h(x)
        for(i=1; i<dim[0]-1; i++) {
            for(j=0; j<dim[1]; j++) {
                for(k=0; k<dim[2]; k++) {
                    xres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        1*vol[k*dim[1]*dim[0] + j*dim[0] + (i-1)]
                        + 2*vol[k*dim[1]*dim[0] + j*dim[0] + i]
                        + 1*vol[k*dim[1]*dim[0] + j*dim[0] + (i+1)];
                }
            }
        }

        // h'(y)
        for(i=0; i<dim[0]; i++) {
            for(j=1; j<dim[1]-1; j++) {
                for(k=0; k<dim[2]; k++) {
                    yres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        - 1*xres[k*dim[1]*dim[0] + (j-1)*dim[0] + i]
                        + 1*xres[k*dim[1]*dim[0] + (j+1)*dim[0] + i];
                }
            }
        }

        // h(z)
        for(i=0; i<dim[0]; i++) {
            for(j=0; j<dim[1]; j++) {
                for(k=1; k<dim[2]-1; k++) {
                    zres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        1*yres[(k-1)*dim[1]*dim[0] + j*dim[0] + i]
                        + 2*yres[k*dim[1]*dim[0] + j*dim[0] + i]
                        + 1*yres[(k+1)*dim[1]*dim[0] + j*dim[0] + i];
                }
            }
        }

        return zres;
    },

    zder: function zder(vol, dim) {
        let i, j, k;
        let xres = new Float32Array(dim[0]*dim[1]*dim[2]);
        let yres = new Float32Array(dim[0]*dim[1]*dim[2]);
        let zres = new Float32Array(dim[0]*dim[1]*dim[2]);

        // h(x)
        for(i=1; i<dim[0]-1; i++) {
            for(j=0; j<dim[1]; j++) {
                for(k=0; k<dim[2]; k++) {
                    xres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        1*vol[k*dim[1]*dim[0] + j*dim[0] + (i-1)]
                        + 2*vol[k*dim[1]*dim[0] + j*dim[0] + i]
                        + 1*vol[k*dim[1]*dim[0] + j*dim[0] + (i+1)];
                }
            }
        }

        // h(y)
        for(i=0; i<dim[0]; i++) {
            for(j=1; j<dim[1]-1; j++) {
                for(k=0; k<dim[2]; k++) {
                    yres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        1*xres[k*dim[1]*dim[0] + (j-1)*dim[0] + i]
                        + 2*xres[k*dim[1]*dim[0] + j*dim[0] + i]
                        + 1*xres[k*dim[1]*dim[0] + (j+1)*dim[0] + i];
                }
            }
        }

        // h(z)
        for(i=0; i<dim[0]; i++) {
            for(j=0; j<dim[1]; j++) {
                for(k=1; k<dim[2]-1; k++) {
                    zres[k*dim[1]*dim[0] + j*dim[0] + i] =
                        - 1*yres[(k-1)*dim[1]*dim[0] + j*dim[0] + i]
                        + 1*yres[(k+1)*dim[1]*dim[0] + j*dim[0] + i];
                }
            }
        }

        return zres;
    }
}

var HDDIRandom = {
    /**
     * @desc https://en.wikipedia.org/wiki/Marsaglia_polar_method
     */
    nrandomPrecomputed: 0,
    nrandomIsPrecomputedReady: false,
    nrandom: function nrandom() {
        if(this.nrandomIsPrecomputedReady) {
            this.nrandomIsPrecomputedReady = false;

            return this.nrandomPrecomputed;
        } else {
            let x1, x2, s, y1;
            do {
                x1 = 2*Math.random() - 1;
                x2 = 2*Math.random() - 1;
                s = x1 * x1 + x2 * x2;
            } while ( s >= 1 || s === 0 );
            s = Math.sqrt( (-2*Math.log(s) ) / s );
            y1 = x1 * s;
            this.nrandomPrecomputed = x2 * s;
            this.nrandomIsPrecomputedReady = true;

            return y1;
        }
    },

    randomDirectionNaive: function randomDirection() {
        let v = {
            x: Math.random() - 0.5,
            y: Math.random() - 0.5,
            z: Math.random() - 0.5
        };
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: v.x/lv,
            y: v.y/lv,
            z: v.z/lv
        };
    },

    /**
      * @desc Muller's method from http://mathworld.wolfram.com/SpherePointPicking.html
      */
    randomDirectionMuller: function randomDirection() {
        let v = {
            x: this.nrandom(),
            y: this.nrandom(),
            z: this.nrandom()
        };
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: v.x/lv,
            y: v.y/lv,
            z: v.z/lv
        };
    },

    /**
      * @desc Marsaglia's method from http://mathworld.wolfram.com/SpherePointPicking.html
      */
    randomDirection: function randomDirection() {
        let x1, x2, s;
        do {
            x1 = 2*Math.random() - 1;
            x2 = 2*Math.random() - 1;
            s = x1 * x1 + x2 * x2;
        } while ( s >= 1 || s === 0);
        return v = {
            x: 2 * x1 * Math.sqrt(1 - s),
            y: 2 * x2 * Math.sqrt(1 - s),
            z: 1 - 2 * s
        };
    }
}

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
(function () {
    "use strict";

    console.log("Two folds test, with different gradients");

    let hddi = new HDDI();

    const params = {
        w: 0.9, // stiffness parameter
        step: 0.5, // step size
        minFibLength: 10,
        ns: 1e+7, // number of streamlines to throw
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

    const wdir = 'experiments/06-twofolds-different-gradients/results/';
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
    let mask = ident.map((v) => (v === 1));
    saveNiftiData(mask, dim, wdir + 'mask.nii.gz');
    
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
        mrtrix.tckgen([wdir + 'dwi.mif', wdir + 'streamlines.50k-det.tck', '-algorithm', 'Tensor_Det', '-seed_image', wdir + 'mask.nii.gz', '-mask', wdir + 'mask.nii.gz', '-select', 50000, '-force']);
    });
} ());