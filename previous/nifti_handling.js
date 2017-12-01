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

var HDDI = function HDDI() {
"use strict";
    var me = {
        loadMask: function loadMask( e ) {
            console.log( '%cload mask', 'color: green; ' );
            console.log( 'loaded element (mask): ' + e );
            var file = e.target.files[0];
            if ( !file ) {

                return;
            }

            nifti1 = new Nifti();
            nifti1.loadFile( file, maskLoadCallback );
        },

        maskLoadCallback: function maskLoadCallback() {
            console.log( '%cmask load callback', 'color: green; ' );

            var canvas = document.getElementById( 'c1' );
            var dims = nifti1.getDims();
            console.log( dims );
            canvas.width = dims.nx;
            canvas.height = dims.nz;
            var ctx = canvas.getContext('2d');

            var imageData = ctx.getImageData( 0, 0, dims.nx, dims.nz );
            var slice = nifti1.getImage( 'coronal', 60 );
            for( var i = 0; i < imageData.data.length; ++i ) {
                imageData.data[i] = slice.data[i];
            }
            ctx.putImageData( imageData, 0, 0 );
        },

        /**
         * @todo Rename invert as invertMask
         */
        invert: function invert() {
            var dims = nifti1.getDims();

            for( var z = 0; z < dims.nz; ++z ) {
                for( var y = 0; y < dims.ny; ++y ) {
                    for( var x = 0; x < dims.nx; ++x ) {
                        var val = nifti1.getValue( x, y, z ); //if mask --> val is 0 or 1
                        if( val > 0 ) {
                            nifti1.setValue( x, y, z, 255 - val );
                        }
                        else {
                            nifti1.setValue( x, y, z, [0, 0, 0] );
                        }
                        //nifti1.setValue( x, y, z, [val/255.0, 0, 0] );
                    }
                }
            }
        },

        /**
          * @desc save nifti3 to only set every visited voxel to 1 to 'reconstruct surface of fibers' afterwards.
          */
        saveFile: function saveFile() {
            if ( nifti3 ) {
                console.log( '%csave nifti3 (voxels set to 1)', 'color: green; ' );
                var filename = $('#input-fileName').val()
                var arrayBuffer = nifti3.getRawData();
                var blob = new Blob([arrayBuffer], { type: 'application/octet-binary' } );
                saveAs( blob, filename );
            }
        }
    };

    return me;
}
