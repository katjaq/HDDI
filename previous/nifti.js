/**
* A class to handle niftii files and provide slices for web visualisation
*
* @adapted from Ralph Schurade
* @license MIT License
*
*/
let Nifti = function() {
    "use strict";
    let me = {
        data: [],
        rawData: null,
        hdr: {},
        dimX: 0,
        dimY: 0,
        dimZ: 0,
        max: -100000,
        min: 100000,
        zero: 0,
        loaded: false,

        download: function( url, callback ) {
            var xhr = new XMLHttpRequest();
            xhr.open( 'GET', url, true );
            xhr.responseType = 'arraybuffer';

            xhr.onload = function() {
                me.rawData = me.response;
                me.data = new DataView( me.response ); // me.response == uInt8Array.buffer
                me.parseHeader();
                me.calcMinMax();

                me.loaded = true;
                if ( callback ) callback();
            };
            xhr.send();
        },

        loadFile: function( file, callback ) {
            var reader = new FileReader();
            reader.onload = function(e) {
                me.rawData = e.target.result;
                me.data = new DataView( e.target.result );
                me.parseHeader();
                me.calcMinMax();

                me.loaded = true;
                if ( callback ) callback();
            };
            reader.readAsArrayBuffer(file);
        },

        /**
         * @func makeNew
         * @desc Make a new nifti object
         * @param number nx X dimension
         * @param number dx Voxel size in the x dimension
         */
        makeNew: function( nx, ny, nz, dx, dy, dz, nt, datatype ) {
            var size = 352 + ( nx * ny * nz * nt * me.sizeOf( datatype ) );
            /*create byteArray*/
            me.rawData = new Uint8Array( size );
            me.data = new DataView( me.rawData.buffer );
            me.data.setInt32( 0, 348, true );
            me.data.setUint8( 39, 1 ); //toDo: realValue? to replace 1
            me.data.setInt16( 40, 5, true );
            me.data.setInt16( 42, nx, true );
            me.data.setInt16( 44, ny, true );
            me.data.setInt16( 46, nz, true );
            me.data.setInt16( 48, nt, true );
            me.data.setInt16( 50, nt, true );
            /* missing! fill in!! */
            me.data.setInt16( 68, 1007, true );
            me.data.setInt16( 70, me.datatypeCode( datatype ), true );
            /* bitpix... ( 72, */
            me.data.setFloat32( 80, dx, true );
            me.data.setFloat32( 84, dy, true );
            me.data.setFloat32( 88, dz, true );
            /* plus 5 more... float, should be 0 */
            me.parseHeader();
        },

        sizeOf: function sizeOf( datatype ) {
            switch( datatype ) {
            case 'UINT8':
                return 1;
            case 'INT16':
                return 2;
            case 'INT32':
                return 4;
            case 'FLOAT32':
                return 4;
            default:
                return -1;
            }
        },

        datatypeCode: function datatypeCode( datatype ) {
            switch( datatype ) {
            case 'UINT8':
                return 2;
            case 'INT16':
                return 4;
            case 'INT32':
                return 8;
            case 'FLOAT32':
                return 16;
            default:
                return -1;
            }
        },

        parseHeader: function parseHeader() {
            let i = 0;
            me.hdr.sizeof_hdr = me.data.getInt32( 0, true ); // 0
            me.hdr.data_type = []; // 4
            for ( i = 0; i < 10; ++i ) me.hdr.data_type.push( me.data.getUint8(  4 + i ) );
            me.hdr.db_name = []; // 14
            for ( i = 0; i < 18; ++i ) me.hdr.db_name.push( me.data.getUint8(  14 + i ) );
            me.hdr.extents = me.data.getInt32( 32, true ); // 32
            me.hdr.session_error = me.data.getInt16( 36, true ) // 36
            me.hdr.regular = me.data.getUint8( 38 ); // 38
            me.hdr.dim_info = me.data.getUint8( 39 ); // 39
            me.hdr.dim = []; // 40
            for ( i = 0; i < 8; ++i ) me.hdr.dim.push( me.data.getInt16(  40 + i * 2, true ) );
            me.hdr.intent_p1  = me.data.getFloat32( 56, true );
            me.hdr.intent_p2  = me.data.getFloat32( 60, true );
            me.hdr.intent_p3  = me.data.getFloat32( 64, true );
            me.hdr.intent_code  = me.data.getInt16( 68, true );
            me.hdr.datatype = me.data.getInt16( 70, true );
            me.hdr.bitpix = me.data.getInt16( 72, true );
            me.hdr.slice_start = me.data.getInt16( 74, true );
            me.hdr.pixdim = [];
            for ( i = 0; i < 8; ++i ) me.hdr.pixdim.push( me.data.getFloat32(  76 + i * 4, true ) );
            me.hdr.vox_offset = me.data.getFloat32( 108, true ); // 108
            me.hdr.scl_slope  = me.data.getFloat32( 112, true );
            me.hdr.scl_inter  = me.data.getFloat32( 116, true );
            me.hdr.slice_end = me.data.getInt16( 120, true );
            me.hdr.slice_code  = me.data.getUint8( 122 );
            me.hdr.xyzt_units  = me.data.getUint8( 123 );
            me.hdr.cal_max = me.data.getFloat32( 124, true );
            me.hdr.cal_min = me.data.getFloat32( 128, true );
            me.hdr.slice_duration = me.data.getFloat32( 132, true );
            me.hdr.toffset = me.data.getFloat32( 136, true );
            me.hdr.glmax = me.data.getInt32( 140, true );
            me.hdr.glmin = me.data.getInt32( 144, true );
            me.hdr.descrip = []; // 148
            for ( i = 0; i < 80; ++i ) me.hdr.descrip.push( me.data.getUint8(  148 + i ) );
            me.hdr.aux_file = []; // 228
            for ( i = 0; i < 24; ++i ) me.hdr.aux_file.push( me.data.getUint8(  228 + i ) );
            me.hdr.qform_code  = me.data.getInt16( 252, true );
            me.hdr.sform_code  = me.data.getInt16( 254, true );
            me.hdr.quatern_b  = me.data.getFloat32( 256, true );
            me.hdr.quatern_c  = me.data.getFloat32( 260, true );
            me.hdr.quatern_d  = me.data.getFloat32( 264, true );
            me.hdr.qoffset_x  = me.data.getFloat32( 268, true );
            me.hdr.qoffset_y  = me.data.getFloat32( 272, true );
            me.hdr.qoffset_z  = me.data.getFloat32( 276, true );
            me.hdr.srow_x = []; // 280
            for ( i = 0; i < 4; ++i ) me.hdr.srow_x.push( me.data.getFloat32(  280 + i * 4, true ) );
            me.hdr.srow_y = []; // 296
            for ( i = 0; i < 4; ++i ) me.hdr.srow_y.push( me.data.getFloat32(  296 + i * 4, true ) );
            me.hdr.srow_z = []; // 312
            for ( i = 0; i < 4; ++i ) me.hdr.srow_z.push( me.data.getFloat32(  312 + i * 4, true ) );
            me.hdr.intent_name = []; // 328
            for ( i = 0; i < 16; ++i ) me.hdr.intent_name.push( me.data.getUint8(  328 + i ) );
            me.hdr.magic = []; // 344
            for ( i = 0; i < 4; ++i ) me.hdr.magic.push( me.data.getUint8(  344 + i ) );

            me.dimX = me.hdr.dim[1];
            me.dimY = me.hdr.dim[2];
            me.dimZ = me.hdr.dim[3];
        },

        calcMinMax: function calcMinMax() {
            let i = null;
            switch( me.hdr.datatype ) {
            case 2: {
                for ( i = 348; i < me.data.byteLength; ++i ) {
                    if ( me.data.getUint8( i ) < me.min ) me.min = me.data.getUint8( i );
                    if ( me.data.getUint8( i ) > me.max ) me.max = me.data.getUint8( i );
                }
                console.log( 'min: ' + me.min + ' max: ' + me.max );
                //me.min = 0;
                //me.max = 255;
            }
            break;
            case 16: {
                for ( i = 348; i < me.data.byteLength; i+=4 ) {
                    if ( me.data.getFloat32( i ) < me.min ) me.min = me.data.getFloat32( i );
                    if ( me.data.getFloat32( i ) > me.max ) me.max = me.data.getFloat32( i );
                }
                //console.log( 'min: ' + me.min + ' max: ' + me.max );

                var div = me.max - me.min;
                me.zero = ( 0 - me.min ) / div;
                for ( var j = 348; j < me.data.length; j+=4 ) {
                    me.data.setFloat32(j, ( me.data.getFloat32(j) - me.min ) / div );
                }
            }
            break;
            default:
                console.log( 'Nifti calcMinMax(): datatype ' + me.hdr.datatype + ' not defined' );
            }
        },

        loadFinished: function loadFinished() {
            return me.loaded;
        },

        getImage: function getImage(orient, pos) {
            if ( !me.loaded ) console.log( 'DEBUG nifti file not finished loading');
            if ( orient === 'sagittal' && pos > me.hdr.dim[1] ) pos = 0;
            if ( orient === 'coronal' && pos > me.hdr.dim[2] ) pos = 0;
            if ( orient === 'axial' && pos > me.hdr.dim[3] ) pos = 0;

            if ( me.hdr.datatype === 2 ) {
                if (me.hdr.dim[4] === 1 ) {
                    return me.getImageGrayByte(orient,pos);
                }
                if (me.hdr.dim[4] === 3) {
                    return me.getImageRGBByte(orient,pos);
                }
            }
            else if ( me.hdr.datatype === 16 ) {
                if (me.hdr.dim[4] === 1 ) {
                    return me.getImageGrayFloat(orient,pos);
                }
            }
        },

        getImageGrayByte: function getImageGrayByte(orient, pos) {
            const c2d = document.createElement('canvas');

            ( orient === 'sagittal' ) ? c2d.width = me.hdr.dim[2] : c2d.width = me.hdr.dim[1];
            ( orient === 'axial' ) ? c2d.height = me.hdr.dim[2] : c2d.height = me.hdr.dim[3];
            var ctx = c2d.getContext('2d');
            var imageData = ctx.getImageData(0, 0, c2d.width, c2d.height);
            let x = null;
            let y = null;
            let z = null;
            let col = null;
            let index = null;

            if ( orient === 'axial' ) {
                for( x = 0; x < me.dimX; ++x )
                    for( y = 0; y < me.dimY; ++y )
                    {
                        col = me.data.getUint8( me.getId(x,(me.dimY-1)-y,pos) );
                        index = 4 * (y * imageData.width + x);
                        me.setImgData( index, col );
                    }
            }

            if ( orient === 'coronal' ) {
                for( x = 0; x < me.dimX; ++x )
                    for( z = 0; z < me.dimZ; ++z )
                    {
                        col = me.data.getUint8( me.getId( x, pos,(me.dimZ-1)-z) );
                        index = 4 * (z * imageData.width + x);
                        me.setImgData( index, col );
                    }
            }

            if ( orient === 'sagittal' ) {
                for( y = 0; y < me.dimY; ++y )
                    for( z = 0; z < me.dimZ; ++z )
                    {
                        col = me.data.getUint8( me.getId(parseInt(pos),y,(me.dimZ-1)-z) );
                        index = 4 * (z * imageData.width + y);
                        me.setImgData( index, col );
                    }
            }
            ctx.putImageData( imageData, 0, 0 );

            return imageData;

        },

/*
        setImgData: function setImgData( id, col ) {
            imageData.data[id] = col;
            imageData.data[id+1] = col;
            imageData.data[id+2] = col;
            imageData.data[id+3] = ( col > 0 ) ? 255 : 0;
        }
*/

        getId: function getId( x,y,z ) {
            return 352 + x + (y * me.hdr.dim[1]) + (z * me.hdr.dim[1] * me.hdr.dim[2]);
        },

        getIdFloat: function getIdFloat( x,y,z ) {
            return 352 + ( x + (y * me.hdr.dim[1]) + (z * me.hdr.dim[1] * me.hdr.dim[2]) ) * 4;
        },

        getValueId: function getValueId( id ) {
            switch( me.hdr.datatype ) {
                //UINT8
            case 2:
                return me.data.getUint8( id );
                //FLOAT32
            case 16:
                if( me.hdr.dim[5] == 1 ) {
                    return me.data.getFloat32( id, true );
                }
                if( me.hdr.dim[5] == 3 ) {
                    var out = [];
                    var blocksize = me.hdr.dim[1] * me.hdr.dim[2] * me.hdr.dim[3] * 4;
                    out[0] = me.data.getFloat32( id, true );
                    out[1] = me.data.getFloat32( id + blocksize, true );
                    out[2] = me.data.getFloat32( id + blocksize*2, true );
                    return out;
                }
                break;
            default:
                    console.log( "Nifti getValue(): datatype not defined" );
            }
        },

        getValue: function getValue( x, y, z ) {
            switch( me.hdr.datatype ) {
            //UINT8
            case 2:
                return me.getValueId( me.getId( x,y,z ) );
            //FLOAT32
            case 16:
                return me.getValueId( me.getIdFloat( x, y, z ) );
            }
        },

        setValue: function setValue( x, y, z, value ) {
            switch( me.hdr.datatype ) {
            case 2:
                me.data.setUint8( me.getId( x, y , z ), value );
                break;
            case 16:
                if( me.hdr.dim[5] == 1 ) {
                    me.data.setFloat32( me.getIdFloat( x, y, z ), value, true );
                }
                if( me.hdr.dim[5] == 3 ) {
                    var blocksize = me.hdr.dim[1] * me.hdr.dim[2] * me.hdr.dim[3] * 4;
                    me.data.setFloat32( me.getIdFloat( x, y, z ), value[0], true );
                    me.data.setFloat32( me.getIdFloat( x, y, z ) + blocksize, value[1], true );
                    me.data.setFloat32( me.getIdFloat( x, y, z ) + blocksize*2, value[2], true );
                }
                break;
            default:
                console.log( 'Nifti setValue(): datatype not defined' );
            }

        },

        getImageRGBByte: function getImageRGBByte(orient, pos) {
            var c2d = document.createElement('canvas');
            ( orient === 'sagittal' ) ? c2d.width = me.hdr.dim[2] : c2d.width = me.hdr.dim[1];
            ( orient === 'axial' ) ? c2d.height = me.hdr.dim[2] : c2d.height = me.hdr.dim[3];
            var ctx = c2d.getContext('2d');
            var imageData = ctx.getImageData(0, 0, c2d.width, c2d.height);
            let x = null;
            let y = null;
            let z = null;
            let r = null;
            let g = null;
            let b = null;
            let index = null;

            var gOff = me.hdr.dim[1] * me.hdr.dim[2] * me.hdr.dim[3];
            var bOff = 2 * gOff;

            if ( orient === 'axial' ) {
                for( x = 0; x < me.dimX; ++x )
                    for( y = 0; y < me.dimY; ++y )
                    {
                        r = me.data.getUint8( me.getId(x,y,pos) );
                        g = me.data.getUint8(parseInt(me.getId(x,y,pos))+parseInt(gOff) );
                        b = me.data.getUint8(parseInt(me.getId(x,y,pos))+parseInt(bOff) );
                        index = 4 * (y * imageData.width + x);
                        me.setImgData( index, r, g, b );
                    }
            }

            if ( orient === 'coronal' ) {
                for( x = 0; x < me.dimX; ++x )
                    for( z = 0; z < me.dimZ; ++z )
                    {
                        r = me.data.getUint8( me.getId(x,pos,z) );
                        g = me.data.getUint8( me.getId(x,pos,z)+gOff );
                        b = me.data.getUint8( me.getId(x,pos,z)+bOff );
                        index = 4 * (z * imageData.width + x);
                        me.setImgData( index, r, g, b );
                    }
            }

            if ( orient === 'sagittal' ) {
                for( y = 0; y < me.dimY; ++y )
                    for( z = 0; z < me.dimZ; ++z )
                    {
                        r = me.data.getUint8( me.getId(pos-1+1,y,z) );
                        g = me.data.getUint8( me.getId(pos-1+1,y,z)+gOff );
                        b = me.data.getUint8( me.getId(pos-1+1,y,z)+bOff );
                        index = 4 * (z * imageData.width + y);
                        me.setImgData( index, r, g, b );
                    }
            }

            return imageData;
        },

/*
        setImgData: function setImgData( id, r, g, b ) {
            imageData.data[id] = r;
            imageData.data[id+1] = g;
            imageData.data[id+2] = b;
            imageData.data[id+3] = 255;
        }
*/

        getImageGrayFloat: function getImageGrayFloat(orient, pos) {
            var c2d = document.createElement('canvas');

            ( orient === 'sagittal' ) ? c2d.width = me.hdr.dim[2] : c2d.width = me.hdr.dim[1];
            ( orient === 'axial' ) ? c2d.height = me.hdr.dim[2] : c2d.height = me.hdr.dim[3];
            var ctx = c2d.getContext('2d');
            var imageData = ctx.getImageData(0, 0, c2d.width, c2d.height);
            let x = 0;
            let y = 0;
            let z = 0;
            let col = 0;
            let index = 0;

            if ( orient === 'axial' ) {
                for( x = 0; x < me.dimX; ++x )
                    for( y = 0; y < me.dimY; ++y )
                    {
                        col = me.data.getFloat32( me.getIdFloat(x,(me.dimY-1)-y,pos) );
                        index = 4 * (y * imageData.width + x);
                        me.setImgData( index, col );
                    }
            }

            if ( orient === 'coronal' ) {
                for( x = 0; x < me.dimX; ++x )
                    for( z = 0; z < me.dimZ; ++z )
                    {
                        col = me.data.getFloat32( me.getIdFloat(x,pos,(me.dimZ-1)-z) );
                        index = 4 * (z * imageData.width + x);
                        me.setImgData( index, col );
                    }
            }

            if ( orient === 'sagittal' ) {
                for( y = 0; y < me.dimY; ++y )
                    for( z = 0; z < me.dimZ; ++z )
                    {
                        col = me.data.getFloat32( me.getIdFloat(pos,y,(me.dimZ-1)-z) );
                        index = 4 * (z * imageData.width + y);
                        me.setImgData( index, col );
                    }
            }

            return imageData;
        },

/*
        setImgData: function setImgData( id, col ) {
            imageData.data[id] = col*255;
            imageData.data[id+1] = col*255;
            imageData.data[id+2] = col*255;
            imageData.data[id+3] = 255;
        },
*/

        getRawData: function getRawData() {
            return me.rawData;
        },

        getMin: function getMin() {
            return me.min;
        },

        getMax: function getMax() {
            return me.max;
        },

        getDims: function getDims() {
            return { 'nx' : me.hdr.dim[1], 'ny' : me.hdr.dim[2], 'nz' : me.hdr.dim[3], 'dx' : me.hdr.pixdim[1], 'dy' : me.hdr.pixdim[2], 'dz' : me.hdr.pixdim[3] };
        }
    };

    return me;
};
