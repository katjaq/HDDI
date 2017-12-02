var HDDIGUI = (function HDDIGUI() {
    "use strict";
    var me = {
        configure: function configure() {
            $('#fiMask').on('change', HDDI.loadMask, false);
            $('#btSurfaceVoxels').on('click', HDDI.identifyVoxels, false);
            $('#btVectors').on('click', HDDI.genVectors, false);
            $('#buttonSave').on('click', HDDI.saveFile, false);
            $('#buttonInvert').on('click', HDDI.invert, false);

            $('#slider1').on("input change", function() {
                if( nifti1 || nifti2 || nifti3 ) {
                    var dims = nifti1.getDims();
                    var maxRange = Math.floor(( dims.nx - 1 ) * 2);
                    console.log( 'maxRange: ', maxRange );
                    var minRange = 1;
                    document.getElementById('slider1').max = maxRange;
                    document.getElementById('slider1').min = minRange;
                    //document.getElementById('slider1').step = 1;
                }
                if( nifti3 ) {
                    sliderCor($(this).val());
                    sliderCor2($(this).val());
                    document.getElementById('slider2').value = $(this).val();
                    sliderCor3($(this).val());
                    document.getElementById('slider3').value = $(this).val();
                }
                else if( nifti2 ) {
                    sliderCor($(this).val());
                    sliderCor2($(this).val());
                    document.getElementById('slider2').value = $(this).val();
                }
                else if( nifti1 ) {
                    sliderCor($(this).val());
                }
                //console.log($(this).next().html($(this).val()));
            }); 

            $('#slider2').on("input change", function() {
                if( nifti1 || nifti2 || nifti3 ) {
                    var dims = nifti1.getDims();
                    var maxRange = Math.floor(( dims.nx - 1 ) * 2);
                    console.log( 'maxRange: ', maxRange );
                    var minRange = 1;
                    document.getElementById('slider2').max = maxRange;
                    document.getElementById('slider2').min = minRange;
                }
                if( nifti3 ) {
                    sliderCor($(this).val());
                    document.getElementById('slider1').value = $(this).val();
                    sliderCor2($(this).val());
                    document.getElementById('slider2').value = $(this).val();
                    sliderCor3($(this).val());
                }
                else if( nifti2 ) {
                    sliderCor($(this).val());
                    document.getElementById('slider1').value = $(this).val();
                    sliderCor2($(this).val());
                }
                //sliderSag($(this).val());
                //console.log($(this).next().html($(this).val()));
            });

            $('#slider3').on("input change", function() {
                if( nifti1 || nifti2 || nifti3 ) {
                    var dims = nifti1.getDims();
                    var maxRange = ( dims.nx - 1 ) * 2;
                    var minRange = 1;
                    document.getElementById('slider3').max = maxRange;
                    document.getElementById('slider3').min = minRange;
                }
                if( nifti3 ) {
                    sliderCor($(this).val());
                    document.getElementById('slider1').value = $(this).val();
                    sliderCor2($(this).val());
                    document.getElementById('slider2').value = $(this).val();
                    sliderCor3($(this).val());
                }
                //sliderAxi($(this).val());
                //console.log($(this).next().html($(this).val()));
            });
        },

        sliderCor: function sliderCor( value ) {
            console.log( '$%csliderCor', 'color: green; ' );
            var canvas = document.getElementById( 'c1' );
            var dims = nifti1.getDims();
            canvas.width = dims.nx;
            canvas.height = dims.nz;
            var ctx = canvas.getContext( '2d' );
            var imageData = ctx.getImageData( 0, 0, dims.nx, dims.nz );
            var slice = nifti1.getImage( 'coronal', value );

            for( var i = 0; i < imageData.data.length; ++i ) {
                imageData.data[i] = slice.data[i];
            }
            ctx.putImageData( imageData, 0, 0 );
        },

        sliderCor2: function sliderCor2( value ) {
            console.log( '%csliderCor2', 'color: green; ' );
            var canvas = document.getElementById( 'c2' );
            var dims = nifti2.getDims();
            canvas.width = dims.nx;
            canvas.height = dims.nz;
            var ctx = canvas.getContext( '2d' );
            var imageData = ctx.getImageData( 0, 0, dims.nx, dims.nz );
            var slice = nifti2.getImage( 'coronal', value );

            for( var i = 0; i < imageData.data.length; ++i ) {
                imageData.data[i] = slice.data[i];
            }
            ctx.putImageData( imageData, 0, 0 );
        },

        sliderCor3: function sliderCor3( value ) {
            console.log( '%csliderCor3', 'color: green; ' );
            var canvas = document.getElementById( 'c3' );
            var dims = nifti3.getDims();
            canvas.width = dims.nx;
            canvas.height = dims.nz;
            var ctx = canvas.getContext( '2d' );
            var imageData = ctx.getImageData( 0, 0, dims.nx, dims.nz );
            var slice = nifti3.getImage( 'coronal', value );

            for( var i = 0; i < imageData.data.length; ++i ) {
                imageData.data[i] = slice.data[i];
            }
            ctx.putImageData( imageData, 0, 0 );
        }
    };

    return me;
}());
