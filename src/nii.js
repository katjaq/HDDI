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

