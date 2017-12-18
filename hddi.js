if(typeof fs === 'undefined') {
    var fs = require('fs');
}
if(typeof Struct === 'undefined') {
    var Struct = require('struct');
}

function readTck(path) {
    data = fs.readFileSync(path);
    let i, j, str;
    let hdr = {};
    let entry;
    j=0;
    for(i=0;i<data.length;i++) {
        if(data[i] === 10) {
            str = data.slice(j,i).toString();
            entry = str.split(/:[ ]*/);
            hdr[entry[0]] = entry[1];
            j=i+1;
            if( str === 'END') {
                break;
            }
        }
    }
    console.log(hdr);
    let count = parseInt(hdr.total_count);
    let offset = parseInt(hdr.file.split(' ')[1]);

    const s = data.slice(offset);

    let x, y, z;
    let strm = [];
    let trk = [];
    for(i=0;i<count;i++) {
        x = s.readFloatLE((3*i+0)*4);
        y = s.readFloatLE((3*i+1)*4);
        z = s.readFloatLE((3*i+2)*4);
        if(isNaN(x)) {
            trk.push(strm);
            strm = [];
            continue;
        }
        strm.push([x, y, z]);
    }

    return {hdr: hdr, tck: trk};
}

function writeTck(streamlines, path) {
    const fd = fs.openSync(path, 'w');
    let i, j, nvertices = 0;
    for(i=0;i<streamlines.length;i++) {
        nvertices += streamlines[i].length;
    }
    let hdrSize = 700;
    let hdr = [
        "mrtrix tracks",
        "init_threshold: 0.1",
        "max_angle: 9",
        "max_dist: 24",
        "max_num_seeds: 10000000",
        "max_num_tracks: " + 50000,
        "max_seed_attempts: 1000",
        "method: TensorDet",
        "min_dist: 1.2",
        "mrtrix_version: f835a76b",
        "rk4: 0",
        "source: simulation.mif",
        "step_size: 0.02",
        "stop_on_all_include: 0",
        "threshold: 0.1",
        "timestamp: 1506408029.8440241814",
        "unidirectional: 0",
        "roi: seed simulation.nii.gz",
        "roi: mask simulation.nii.gz",
        "datatype: Float32LE",
        "file: . " + hdrSize,
        "count: " + streamlines.length,
        "total_count: " + nvertices,
        "END"
    ].join("\n");
    fs.writeSync(fd, hdr);

    let offset = new Uint8Array(hdrSize - hdr.length).fill(0);
    fs.writeSync(fd, offset);

    const nanVal = Math.sqrt(-1);
    const nullVal = null;
    let floatVal = new Float32Array(3);
    for(i=0;i<streamlines.length;i++) {
        for(j=0;j<streamlines[i].length;j++) {
            floatVal[0] = streamlines[i][j][0];
            floatVal[1] = streamlines[i][j][1];
            floatVal[2] = streamlines[i][j][2];
            fs.writeSync(fd, new Buffer(floatVal.buffer));
        }
        floatVal[0] = NaN;
        floatVal[1] = NaN;
        floatVal[2] = NaN;
        fs.writeSync(fd, new Buffer(floatVal.buffer));
    }
    fs.writeSync(fd, new Buffer((new Float32Array(1).fill(null)).buffer));
    fs.closeSync(fd);
}

const TrkHdr = new Struct()
    .chars('id_string', 6)
    .array('dim', 3, 'word16Sle')
    .array('voxel_size', 3, 'floatle')
    .array('origin', 3, 'floatle')
    .word16Sle('n_scalars')
    .chars('scalar_name', 200)
    .word16Sle('n_properties')
    .chars('property_name', 200)
    .array('vox_to_ras', 16, 'floatle')
    .chars('reserved', 444)
    .chars('voxel_order', 4)
    .chars('pad2', 4)
    .array('image_orientation_patient', 6, 'floatle')
    .chars('pad1', 2)
    .word8('invert_x')
    .word8('invert_y')
    .word8('invert_z')
    .word8('invert_xy')
    .word8('invert_yz')
    .word8('invert_zx')
    .word32Sle('n_count')
    .word32Sle('version')
    .word32Sle('hdr_size');

function readTrk(path) {
    data = fs.readFileSync(path);
    TrkHdr.allocate();
    TrkHdr._setBuff(data);
    var hdr = JSON.parse(JSON.stringify(TrkHdr.fields));
    console.log(hdr);

    // read tractography
    s = data.slice(hdr.hdr_size);
    console.log(`n_count=${hdr.n_count}`);
    let i, j, k, l, n;
    i = 0;
    let trk = [];
    for(j=0;j<hdr.n_count;j++) {
        let pr = [];
        let strm=[];

        // read number of vertices in streamline
        n = s.readUInt32LE((i++)*4);

        // read streamline
        for(k=0;k<n;k++) {
            // read x, y and z coordinates of vertices
            strm.push([
                s.readFloatLE((i++)*4),
                s.readFloatLE((i++)*4),
                s.readFloatLE((i++)*4)
            ]);
            // read vertex properties
            for(l=0;l<hdr.n_properties; l++) {
                pr.push(s.readFloatLE((i++)*4));
            }
        }

        trk.push(strm);
    }
    return {hdr: hdr, trk: trk};
}
console.log(typeof fs);
console.log(typeof Struct);

if(typeof fs === 'undefined') {
    var fs = require('fs');
}
if(typeof Struct === 'undefined') {
    var Struct = require('struct');
}
var zlib = require('zlib');

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

function loadNiftiData(path) {
    var niigz = fs.readFileSync(path);

    var pr = new Promise(function (resolve, reject) {
        zlib.gunzip(niigz, function (err, nii) {
            if(err) {
                reject(err);

                return;
            }

            NiiHdr._setBuff(nii);
            var h = JSON.parse(JSON.stringify(NiiHdr.fields));

            if (h.sizeof_hdr != 348) {
                reject(new Error("Can't read Big Endian data"));

                return;
            }

            var vox_offset = h.vox_offset;
            var sizeof_hdr = h.sizeof_hdr;

            const datatype = h.datatype;
            const dim = [h.dim[1],h.dim[2],h.dim[3]];
            const pixdim = [h.pixdim[1],h.pixdim[2],h.pixdim[3]];
            let data, j, tmp;

            switch(datatype) {
                case 2: // UCHAR
                    data = nii.slice(vox_offset);
                    break;
                case 4: // SHORT
                    tmp = nii.slice(vox_offset);
                    data = new Int16Array(dim[0]*dim[1]*dim[2]);
                    for(j = 0; j<dim[0]*dim[1]*dim[2]; j += 1) {
                        data[j] = tmp.readInt16LE(j*2);
                    }
                    break;
                case 8: // INT
                    tmp = nii.slice(vox_offset);
                    data = new Uint32Array(dim[0]*dim[1]*dim[2]);
                    for(j = 0; j<dim[0]*dim[1]*dim[2]; j += 1) {
                        data[j] = tmp.readUInt32LE(j*4);
                    }
                    break;
                case 16: // FLOAT
                    tmp = nii.slice(vox_offset);
                    data = new Float32Array(dim[0]*dim[1]*dim[2]);
                    for(j = 0; j<dim[0]*dim[1]*dim[2]; j += 1) {
                        data[j] = tmp.readFloatLE(j*4);
                    }
                    break;
                case 256: // INT8
                    tmp = nii.slice(vox_offset);
                    data = new Int8Array(dim[0]*dim[1]*dim[2]);
                    for(j = 0; j<dim[0]*dim[1]*dim[2]; j += 1) {
                        data[j] = tmp.readInt8(j);
                    }
                    break;
                case 512: // UINT16
                    tmp = nii.slice(vox_offset);
                    data = new Uint16Array(dim[0]*dim[1]*dim[2]);
                    for(j = 0; j<dim[0]*dim[1]*dim[2]; j += 1) {
                        data[j] = tmp.readUInt16LE(j*2);
                    }
                    break;
                default: {
                    reject(new Error("ERROR: Unknown dataType: " + datatype));

                    return;
                }
            }

            resolve({vol: data, dim: dim});
        });
    });

    return pr;
}


var HDDILinAlg = {
    /**
      * @desc 3D vector addition
      */
    add: function add( a, b ) {
        return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z };
    },

    /**
      * @desc 3D vector subtraction
      */
    sub: function sub( a, b ) {
        return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
    },

    /**
      * @desc 3D vector dot product
      */
    dot: function dot( a, b ) {
        return( a.x * b.x + a.y * b.y + a.z * b.z )
    },

    normalise: function normalise(v, l) {
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: l * v.x / lv,
            y: l * v.y / lv,
            z: l * v.z / lv
        };
    },

    normVec: function normVec(v) {
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: v.x / lv,
            y: v.y / lv,
            z: v.z / lv
        };
    },

    norm: function norm(v) {
        return Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
    },

    scale: function scale(v, l) {
        return {
            x: l * v.x,
            y: l * v.y,
            z: l * v.z
        };
    },

    /**
      * @desc multiplication of matrices a and b
      */
    mult: function mult(a, b) {
        let r,c,i;
        let res = new Array(a.length).fill(0);
        res=res.map(o=>new Array(b[0].length).fill(0));
        for(r=0;r<a.length;r++) {
            for(c=0;c<b[0].length;c++) {
                for(i=0;i<a[0].length;i++) {
                    res[r][c] += a[r][i]*b[i][c];
                }
            }
        }

        return res;
    },

    /**
      * @desc transpose of matrix a
      */
    transpose: function transpose(a) {
        let i,j,r = new Array(a[0].length).fill(0);
        r=r.map(o=>new Array(a.length).fill(0));
        for(i=0;i<a.length;i++)
        for(j=0;j<a[0].length;j++) {
            r[j][i] = a[i][j];
        }
        return r;
    },

    /**
      * @desc make a diagonal matrix from vector v
      */
    diag: function diag(v) {
        let m = new Array(v.length).fill(0);
        m = m.map((o) => new Array(v.length).fill(0));
        v.map((o,i)=>m[i][i]=o);

        return m;
    },

    /**
      * @desc This function comes from the package numeric.js 
      *       Shanti Rao sent me this routine by private email. I had to modify it
      *       slightly to work on Arrays instead of using a Matrix object.
      *       It is apparently translated from http://stitchpanorama.sourceforge.net/Python/svd.py
      */

    svd: function svd(A) {
        //Compute the thin SVD from G. H. Golub and C. Reinsch, Numer. Math. 14, 403-420 (1970)
        let temp;
        let prec= 2.220446049250313e-16; //numeric.epsilon; //Math.pow(2,-52) // assumes double prec
        let tolerance= 1.e-64/prec;
        let itmax= 50;
        let c=0;
        let i=0;
        let j=0;
        let k=0;
        let l=0;

        let u= JSON.parse(JSON.stringify(A)); // numeric.clone(A);
        let m= u.length;
        let n= u[0].length;

        if (m < n) {
            throw "Need more rows than columns";
        }

        let e = new Array(n);
        let q = new Array(n);
        for (i=0; i<n; i++) {
            e[i] = 0.0;
            q[i] = 0.0;
        }
        let v = new Array(n).fill(0);
        v = v.map((o) => new Array(n).fill(0));

        function pythag(a,b) {
            a = Math.abs(a)
            b = Math.abs(b)
            if (a > b)
                return a*Math.sqrt(1.0+(b*b/a/a))
            else if (b == 0.0)
                return a
            return b*Math.sqrt(1.0+(a*a/b/b))
        }

        //Householder's reduction to bidiagonal form

        let f= 0.0;
        let g= 0.0;
        let h= 0.0;
        let x= 0.0;
        let y= 0.0;
        let z= 0.0;
        let s= 0.0;

        for (i=0; i < n; i++) {
            e[i]= g;
            s= 0.0;
            l= i+1;
            for (j=i; j < m; j++)
                s += (u[j][i]*u[j][i]);
            if (s <= tolerance) {
                g= 0.0;
            } else {
                f= u[i][i];
                g= Math.sqrt(s);
                if (f >= 0.0) g= -g;
                h= f*g-s
                u[i][i]=f-g;
                for (j=l; j < n; j++)
                {
                    s= 0.0
                    for (k=i; k < m; k++)
                        s += u[k][i]*u[k][j]
                    f= s/h
                    for (k=i; k < m; k++)
                        u[k][j]+=f*u[k][i]
                }
            }
            q[i]= g
            s= 0.0
            for (j=l; j < n; j++) {
                s= s + u[i][j]*u[i][j]
            }
            if (s <= tolerance) {
                g= 0.0;
            } else {
                f= u[i][i+1]
                g= Math.sqrt(s)
                if (f >= 0.0) g= -g
                h= f*g - s
                u[i][i+1] = f-g;
                for (j=l; j < n; j++) {
                    e[j]= u[i][j]/h;
                }
                for (j=l; j < m; j++) {
                    s=0.0
                    for (k=l; k < n; k++) {
                        s += (u[j][k]*u[i][k])
                    }
                    for (k=l; k < n; k++) {
                        u[j][k]+=s*e[k]
                    }
                }
            }
            y= Math.abs(q[i])+Math.abs(e[i])
            if (y>x)
                x=y
        }

        // accumulation of right hand transformations
        for (i=n-1; i != -1; i+= -1) {
            if (g != 0.0) {
                h= g*u[i][i+1];
                for (j=l; j < n; j++) {
                    v[j][i]=u[i][j]/h;
                }
                for (j=l; j < n; j++) {
                    s=0.0
                    for (k=l; k < n; k++) {
                        s += u[i][k]*v[k][j];
                    }
                    for (k=l; k < n; k++) {
                        v[k][j]+=(s*v[k][i]);
                    }
                }
            }
            for (j=l; j < n; j++) {
                v[i][j] = 0;
                v[j][i] = 0;
            }
            v[i][i] = 1;
            g= e[i]
            l= i
        }
        // accumulation of left hand transformations
        for (i=n-1; i != -1; i+= -1) {
            l= i+1
            g= q[i]
            for (j=l; j < n; j++)
                u[i][j] = 0;
            if (g != 0.0) {
                h= u[i][i]*g
                for (j=l; j < n; j++) {
                    s=0.0
                    for (k=l; k < m; k++) s += u[k][i]*u[k][j];
                    f= s/h
                    for (k=i; k < m; k++) u[k][j]+=f*u[k][i];
                }
                for (j=i; j < m; j++) u[j][i] = u[j][i]/g;
            }
            else
                for (j=i; j < m; j++) u[j][i] = 0;
            u[i][i] += 1;
        }

        // diagonalization of the bidiagonal form
        prec= prec*x
        for (k=n-1; k != -1; k+= -1) {
            for (let iteration=0; iteration < itmax; iteration++)
            { // test f splitting
                let test_convergence = false
                for (l=k; l != -1; l+= -1)
                {
                    if (Math.abs(e[l]) <= prec) {
                        test_convergence= true;
                        break;
                    }
                    if (Math.abs(q[l-1]) <= prec) {
                        break
                    }
                }
                if (!test_convergence)
                { // cancellation of e[l] if l>0
                    c= 0.0
                    s= 1.0
                    let l1= l-1
                    for (i =l; i<k+1; i++)
                    {
                        f= s*e[i]
                        e[i]= c*e[i]
                        if (Math.abs(f) <= prec)
                            break
                        g= q[i]
                        h= pythag(f,g)
                        q[i]= h
                        c= g/h
                        s= -f/h
                        for (j=0; j < m; j++)
                        {
                            y= u[j][l1]
                            z= u[j][i]
                            u[j][l1] =  y*c+(z*s)
                            u[j][i] = -y*s+(z*c)
                        }
                    }
                }
                // test f convergence
                z= q[k]
                if (l == k)
                { //convergence
                    if (z<0.0)
                    { //q[k] is made non-negative
                        q[k]= -z
                        for (j=0; j < n; j++)
                            v[j][k] = -v[j][k]
                    }
                    break  //break out of iteration loop and move on to next k value
                }
                if (iteration >= itmax-1) {
                    throw 'Error: no convergence.'
                }
                // shift from bottom 2x2 minor
                x= q[l]
                y= q[k-1]
                g= e[k-1]
                h= e[k]
                f= ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
                g= pythag(f,1.0)
                if (f < 0.0)
                    f= ((x-z)*(x+z)+h*(y/(f-g)-h))/x
                else
                    f= ((x-z)*(x+z)+h*(y/(f+g)-h))/x
                // next QR transformation
                c= 1.0
                s= 1.0
                for (i=l+1; i< k+1; i++) {
                    g= e[i]
                    y= q[i]
                    h= s*g
                    g= c*g
                    z= pythag(f,h)
                    e[i-1]= z
                    c= f/z
                    s= h/z
                    f= x*c+g*s
                    g= -x*s+g*c
                    h= y*s
                    y= y*c
                    for (j=0; j < n; j++)
                    {
                        x= v[j][i-1]
                        z= v[j][i]
                        v[j][i-1] = x*c+z*s
                        v[j][i] = -x*s+z*c
                    }
                    z= pythag(f,h)
                    q[i-1]= z
                    c= f/z
                    s= h/z
                    f= c*g+s*y
                    x= -s*g+c*y
                    for (j=0; j < m; j++) {
                        y= u[j][i-1]
                        z= u[j][i]
                        u[j][i-1] = y*c+z*s
                        u[j][i] = -y*s+z*c
                    }
                }
                e[l]= 0.0
                e[k]= f
                q[k]= x
            }
        }
        //vt= transpose(v)
        //return (u,q,vt)
        for (i=0;i<q.length; i++)
          if (q[i] < prec) q[i] = 0

        //sort eigenvalues
        for (i=0; i< n; i++) {
        //writeln(q)
         for (j=i-1; j >= 0; j--) {
          if (q[j] < q[i]) {
        //  writeln(i,'-',j)
           c = q[j]
           q[j] = q[i]
           q[i] = c
           for(k=0;k<u.length;k++) {
             temp = u[k][i];
             u[k][i] = u[k][j];
             u[k][j] = temp;
           }
           for(k=0;k<v.length;k++) {
             temp = v[k][i];
             v[k][i] = v[k][j];
             v[k][j] = temp;
           }
        //   u.swapCols(i,j)
        //    v.swapCols(i,j)
           i = j
          }
         }
        }

        return {U:u,S:q,V:v}
    },

    /**
      * @desc The inverse of matrix a, where a=USV', is a^(-1) = VSU': https://www.cse.unr.edu/~bebis/CS791E/Notes/SVD.pdf
      */
    svdinv: function svdinv(a) {
        const {U:u,S:s,V:v} = this.svd(a);
        let i, j;
        let d = [];
        for(i=0;i<s.length;i++) {
            d[i] = v[i].map((o,j)=>o/s[j]);
        }
        return this.mult(d,this.transpose(u));
    }
}

/*
http://www.diffusion-imaging.com/2014/04/from-diffusion-weighted-images-to.html
*/

// const math = require('mathjs');

var HDDIDTI = {
    /**
     * @func eigenvalues
     * @desc Returns the eigenvalues of a symmetric 3x3 matrix
     * @param array d A symmetric 3x3 diffusion matrix
     */
    eigenvalues: function eigenvalues(d) {
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
    },

    /**
     * @func eigenvectors
     * @desc Returns the eigenvectors of a symmetric 3x3 matrix
     * @param array l A vector with 3 eigenvalues
     * @param array d A symmetric 3x3 diffusion matrix
     * @return object An object with three eigenvectors
     */
    eigenvectors: function eigenvectors(l, d) {
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

        return {e1: this.normVec(e1), e2: this.normVec(e2), e3: this.normVec(e3)};
    },

    invHfromG: function invHfromG(G) {
        const invH = this.svdinv([
            [G[0].x*G[0].x, G[0].y*G[0].y, G[0].z*G[0].z, 2*G[0].x*G[0].y, 2*G[0].x*G[0].z, 2*G[0].y*G[0].z],
            [G[1].x*G[1].x, G[1].y*G[1].y, G[1].z*G[1].z, 2*G[1].x*G[1].y, 2*G[1].x*G[1].z, 2*G[1].y*G[1].z],
            [G[2].x*G[2].x, G[2].y*G[2].y, G[2].z*G[2].z, 2*G[2].x*G[2].y, 2*G[2].x*G[2].z, 2*G[2].y*G[2].z],
            [G[3].x*G[3].x, G[3].y*G[3].y, G[3].z*G[3].z, 2*G[3].x*G[3].y, 2*G[3].x*G[3].z, 2*G[3].y*G[3].z],
            [G[4].x*G[4].x, G[4].y*G[4].y, G[4].z*G[4].z, 2*G[4].x*G[4].y, 2*G[4].x*G[4].z, 2*G[4].y*G[4].z],
            [G[5].x*G[5].x, G[5].y*G[5].y, G[5].z*G[5].z, 2*G[5].x*G[5].y, 2*G[5].x*G[5].z, 2*G[5].y*G[5].z]
        ]);

        return invH;
    },

    /**
     * @func diffusionMatrix6
     * @desc Computes a diffusion matrix based on signals measured in 6 different axes
     * @param array S An array of 7 values, first the signal without gradients, followed by signals in 6 directions
     * @param array G An array of 6 directions, each with x, y and z components
     * @return array A diffusion matrix expressed as a vector, d = [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]
     */
    diffusionMatrix6: function diffusionMatrix6(S, invH) {
        let b = 1000; // s/mm2 (changing b only changes the eigenvalues, and thus MD, but not the eigenvectors or the FA)
        let S0 = S[0]; // larger than any Sk
        const Y = [
            Math.log(S0/S[1])/b,
            Math.log(S0/S[2])/b,
            Math.log(S0/S[3])/b,
            Math.log(S0/S[4])/b,
            Math.log(S0/S[5])/b,
            Math.log(S0/S[6])/b
        ];
        return this.mult(invH, h.transpose([Y])).map(o=>o[0]);
    },

    /**
      * @desc fit diffusion tensors to a diffusion weighted volume
      */
    dti: function dti(dwvol, dim, G) {
        let i;
        let dt = [];
        const invH = this.invHfromG(G);
        for(i=0; i<dim[0]*dim[1]*dim[2]; i++) {
            let signal = [];
            for(j=0;j<dwvol.length;j++) {
                signal.push(dwvol[j][i]);
            }
            signal = [Math.max(...signal), ...signal];
            const dif = this.diffusionMatrix6(signal, invH);
            console.log("dif",dif);
            /*
            const {U,S,V} = h.svd([
                [dif[0],dif[3],dif[4]],
                [dif[3],dif[1],dif[5]],
                [dif[4],dif[5],dif[2]]
            ]);
            */
            //dt[i] = this.mult(U, this.diag(S));
            dt[i] = h.svd([
                [dif[0],dif[3],dif[4]],
                [dif[3],dif[1],dif[5]],
                [dif[4],dif[5],dif[2]]
            ]);
        }

        return dt;
    },

    /**
     * @desc Mean diffusivity (or apparent diffusion coefficient)
     */
    md: function md(evals) {
        return (evals.l1 + evals.l2 + evals.l3)/3;
    },

    /**
     * @desc Fractional anisotropy
     */
    fa: function fa(evals) {
        const {l1, l2, l3} = evals;
        const MD = (l1 + l2 + l3)/3;
        return Math.sqrt(
            (3*Math.pow(l1 - MD, 2)
            + Math.pow(l2 - MD, 2)
            + Math.pow(l3 - MD, 2))
            /(2*(l1*l1 + l2*l2 + l3*l3)));
    }
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
        let fibres = [];

        bar.start(this.params.ns, 0);

        // generate fibres
        count = 0;
        while( count < this.params.ns) {
            let fib = [];

            if( count%1e+3 === 0 ) {
                bar.update(count);
            }

            // take random boundary voxel to start fiber
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

                // new direction = mix of old plus random
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

            // stroke the streamline in the density and direction volumes
            for(i = 0; i < fib.length; i++) {
                [pos, v] = fib[i];
                //set every visited voxel to 1
                ix = parseInt(pos.x);
                iy = parseInt(pos.y);
                iz = parseInt(pos.z);
                this.addValue( this.rho, dim, ix, iy, iz, 1 );
                this.addValueN( this.dir, dim, ix, iy, iz, v );
            }

            // if part of the last 5k fibres, store it
            if(this.params.ns - count <5000) {
                let fib2 = [];
                for(i=0;i<fib.length;i+=Math.min(1, fib.length/10)) {
                    fib2.push([fib[i][0].x,fib[i][0].y,fib[i][0].z]);
                }
                fibres.push(fib2);
            }

            count += 1;
        }
        bar.update(this.params.ns);
        bar.stop();
        if( this.debug > 2 ) {
            console.log( '%cgenerate voxels finito', 'color: light-green; ' );
        }
        console.log("Min and Max fibre length:", min, max);
        writeTck(fibres, "test.tck");
    },

    /**
      * @func genStickyStreamlines
      * @desc Generates random fibres which tend to follow the previous directions
      */
    genStickyStreamlines: function genStickyStreamlines(vol, dim) {
        console.log( '%cgenerate fibres', 'color: green; ' );

        let i;
        let count;
        let pos, rindex;
        let ix, iy, iz;
        let length, max, min;
        let v, v2, v3;
        let dot;
        let bar = new progress.Bar({}, progress.Presets.shades_classic);
        let fibres = [];

        bar.start(this.params.ns, 0);

        // initialise substrate orientation to random
        let md, mdir = [];
        let a;
        this.params.anis = [];
        let anis = this.params.anis;
        for(ix=0;ix<dim[0];ix++) {
            for(iy=0;iy<dim[1];iy++) {
                for(iz=0;iz<dim[2];iz++) {
                    if( this.getValue( vol, dim, ix, iy, iz ) > 0 ) {
                        mdir[iz*dim[1]*dim[0]+iy*dim[0]+ix] = this.randomDirection();//{x:ix/dim[0],y:0,z:1}; // this.randomDirection();
                        this.params.anis[iz*dim[1]*dim[0]+iy*dim[0]+ix] = 0.2;
                    } else {
                        mdir[iz*dim[1]*dim[0]+iy*dim[0]+ix] = {x:0,y:0,z:0};
                        this.params.anis[iz*dim[1]*dim[0]+iy*dim[0]+ix] = 0.2;
                    }
                }
            }
        }

        // generate fibres
        count = 0;
        while( count < this.params.ns) {
            let fib = [];

            if( count%1e+3 === 0 ) {
                bar.update(count);
            }

            // take random boundary voxel to start fiber
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

                // local main orientation
                v3 = this.scale(mdir[parseInt(pos.z)*dim[1]*dim[0] + parseInt(pos.y)*dim[0] + parseInt(pos.x)], this.params.step);

                //  make consistent with fibre direction
                dot = this.dot(v3, v);
                if(dot<0) {
                    v3 = this.scale(v3, -1);
                }

                // substrate anisotropy
                a = Math.max(0, this.params.anis[parseInt(pos.z)*dim[1]*dim[0] + parseInt(pos.y)*dim[0] + parseInt(pos.x)] - 0.2)/0.8;

                // combine own direction, substrate orientation and random orientation
                if(count<this.params.ns-1e+6) {
                    a = 0;
                }

                v = this.normalise({
                    x: a*v3.x + (1-a)*(this.params.w*v.x + (1-this.params.w)*v2.x) ,
                    y: a*v3.y + (1-a)*(this.params.w*v.y + (1-this.params.w)*v2.y) ,
                    z: a*v3.z + (1-a)*(this.params.w*v.z + (1-this.params.w)*v2.z) ,
                }, this.params.step);

                // advance streamline
                pos = {
                    x: pos.x + v.x,
                    y: pos.y + v.y,
                    z: pos.z + v.z
                };

                // store new vertex
                fib.push([pos, v]);

                // update fibre length
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

            // stroke the streamline in the density and direction volumes
            // update main direction and anisotropy
            for(i = 0; i < fib.length; i++) {
                [pos, v] = fib[i];

                // voxel coordinates
                ix = parseInt(pos.x);
                iy = parseInt(pos.y);
                iz = parseInt(pos.z);

                //set every visited voxel to 1
                this.addValue( this.rho, dim, ix, iy, iz, 1 );

                // update the directions
                this.addValueN( this.dir, dim, ix, iy, iz, v );

                // update substrate orientation
                md = mdir[iz*dim[1]*dim[0] + iy*dim[0] + ix];
                a = this.params.anis[iz*dim[1]*dim[0] + iy*dim[0] + ix];
                dot = this.dot(md, v);
                if(dot<0) {
                    v = this.scale(v, -1);
                }
                md = {
                    x: this.params.dmass * md.x + (1-this.params.dmass)*v.x,
                    y: this.params.dmass * md.y + (1-this.params.dmass)*v.y,
                    z: this.params.dmass * md.z + (1-this.params.dmass)*v.z
                }
                a = this.params.dmass * a + (1-this.params.dmass) * this.norm(this.sub(v, this.scale(md, Math.abs(this.dot(v,md))/this.dot(md,md))));
                
                mdir[iz*dim[1]*dim[0] + iy*dim[0] + ix] = md;
                this.params.anis[iz*dim[1]*dim[0] + iy*dim[0] + ix] = a;
            }

            // if part of the last 5k fibres, store it
            if(this.params.ns - count <5000) {
                let fib2 = [];
                for(i=0;i<fib.length;i+=Math.min(1, fib.length/10)) {
                    fib2.push([fib[i][0].x,fib[i][0].y,fib[i][0].z]);
                }
                fibres.push(fib2);
            }

            count += 1;
        }
        bar.update(this.params.ns);
        bar.stop();
        if( this.debug > 2 ) {
            console.log( '%cgenerate voxels finito', 'color: light-green; ' );
        }
        console.log("Min and Max fibre length:", min, max);
        writeTck(fibres, "test-sticky.tck");
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
      * @desc Muller's method adapted to an ellipsoid. Goes back to the case of a sphere if a=b=c
      */
    randomDirectionMullerEllipsoid: function randomDirectionMullerEllipsoid(a, b, c) {
        let v = {
            x: a * this.nrandom(),
            y: b * this.nrandom(),
            z: c * this.nrandom()
        };
        let lv = Math.sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
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

console.log('\nExtending HDDI from HDDILinAlg');
props=Object.keys(HDDILinAlg);
for(ind in props) {
    prop = props[ind];
    if(typeof prop !== 'undefined') {
        HDDI.prototype[prop] = HDDILinAlg[prop];
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

    console.log("Homogeneous sphere test");

    let hddi = new HDDI();

    const params = {
        w: 0.9, // stiffness parameter
        wh: 1-1e-5, // density homegeneity parameter
        step: 1, // step size
        minFibLength: 30,
        ns: 1e+6, // number of streamlines to throw
        dir: [      // array containing the diffusion directions 'measured' read in from bvec file; no more used
            {x: 0.817324, y: -0.49673, z: -0.29196},
            {x: 0.465087, y: -0.03533, z: 0.88456},
            {x: 0.820439, y: -0.31517, z: 0.477018},
            {x: -0.80334, y: 0.593293, z: -0.05141},
            {x: -0.15636, y: 0.788990, z: -0.59418},
            {x: -0.11253, y: -0.34483, z: -0.93189}
        ]
    };
    hddi.params = params;

    const wdir = 'experiments/04-homogeneous-sphere/results/';
    const exec = require('child_process').execSync;

    // generate sphere
    const dim = [80, 80, 80];
    let vol = new Int32Array(dim[0]*dim[1]*dim[2]);
    let i, j, k;
    let r = 35;
    for(i=0;i<dim[0];i++) {
        for(j=0;j<dim[1];j++) {
            for(k=0;k<dim[2];k++) {
                if( Math.pow(i-dim[0]/2,2) + Math.pow(j-dim[1]/2,2) + Math.pow(k-dim[2]/2,2) < r*r) {
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

    // initialise rho and dir
    hddi.initialise(dim);

    // generate a fibre density potential field
    let val;
    for(i=0;i<dim[0];i++) {
        for(j=0;j<dim[1];j++) {
            for(k=0;k<dim[2];k++) {
                val = Math.pow(i-dim[0]/2,2) + Math.pow(j-dim[1]/2,2) + Math.pow(k-dim[2]/2,2);
                if( val < r*r) {
                    hddi.setValue(hddi.rho, dim, i, j, k, 35*35-val);
                } else {
                    hddi.setValue(hddi.rho, dim, i, j, k, 0);
                }
            }
        }
    }

    // generate homogeneous streamlines
    hddi.genHomogeneousStreamlines(ident, dim);
    saveNiftiData(hddi.rho, dim, wdir + 'rho.nii.gz');

    // create a b0 image as max(dir)
    const b0 = hddi.computeB0(dim);

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