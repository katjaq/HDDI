"use strict";

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


var nifti1; 			//mask or T1 in format .nii
var nifti2;				//identify surface voxels (set to bv) and inside volume (set 1 if > 0) and outside (0), 
var nifti3; 			//contains value 1 for every visited voxel 
var nifti4; 			//stores direction values from one visit

var bv = 200; 			//border value
var bvox = [];			//bordervoxels

var v2;					//new random vector
var v;					//final direction vector
var w = 0.5;			//stiffness parameter
var step = 0.1;			//step size

var np = 1000;			//number of thrown particles

var nd = 3;				//number of directions

var vol = [];			//new volum array (storing direction value sum of all passed fibers) 

var vvC = [];			//voxelVisitCount: array (blocksize) containing integer: how many times has voxel been passed by a fiber

var dir = [];			//array containing the diffusion directions 'measured' read in from bvec file; no more used

var debug = 0;


/* usage */
/*
load html file
click load data 			--> choose either T1 or mask file in .nii format
click surface voxels 		--> identify surface voxels
click generate vectors 		--> computes random fibers
save the nifti file 		--> nifti3 will contain value of 1 in every visited voxel;
								nifti4 will contain color direction value
open saved nifti in fslview --> if nifti4: you can view it as *i* --> diffusion tensor, RGB lines 
								if nifti3: save it from fsl (nii.gz)
open in mango 				--> build surface to check 'fibers'
*/


/* next steps */
/*
keep vistitation map 
repell new fibers based on density
*/






function loadMask( e ) {
	console.log( '%cload mask', 'color: green;' );
	console.log( 'loaded element (mask): ' + e );
	var file = e.target.files[0];
	if ( !file ) {
		return;
	}
	
	nifti1 = new Nifti();
	nifti1.loadFile( file, maskLoadCallback );
	
}

function maskLoadCallback() {
	console.log( '%cmask load callback', 'color: green;' );
	
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
	ctx.putImageData( imageData, 0,0 );
}

function invert() {
	var dims = nifti1.getDims();

	for( var z = 0; z < dims.nz; ++z ) {
		for( var y = 0; y < dims.ny; ++y ) {
			for( var x = 0; x < dims.nx; ++x ) {
				var val = nifti1.getValue( x, y, z ); //if mask --> val is 0 or 1
				if( val > 0 ) {
					nifti1.setValue( x, y, z, 255 - val );	
				}
				else {
					nifti1.setValue( x, y, z, [0,0,0] );
				}
				//nifti1.setValue( x, y, z, [val/255.0,0,0] );	
			}
		}
	}
}

function sliderCor( value ) {
	console.log( '$%cliderCor', 'color: green;' );
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
	ctx.putImageData( imageData, 0,0 );
}

function sliderCor2( value ) {
	console.log( '%csliderCor2', 'color: green;' );
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
	ctx.putImageData( imageData, 0,0 );
}

function sliderCor3( value ) {
	console.log( '%csliderCor3', 'color: green;' );
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
	ctx.putImageData( imageData, 0,0 );
}

/* save nifti 4 to save color values */
/*
function saveFile() {
	if ( nifti4 ) {
		console.log( '%csave nifti4 (voxels set to direction)', 'color: green;' );
		var filename = $('#input-fileName').val()
		var arrayBuffer = nifti4.getRawData();
		var blob = new Blob([arrayBuffer], {type: 'application/octet-binary'} );
		saveAs( blob, filename );
	}
}
*/

/* save nifti3 to only set every visited voxel to 1 to 'reconstruct surface of fibers' afterwards */
function saveFile() {
	if ( nifti3 ) {
		console.log( '%csave nifti3 (voxels set to 1)', 'color: green;' );
		var filename = $('#input-fileName').val()
		var arrayBuffer = nifti3.getRawData();
		var blob = new Blob([arrayBuffer], {type: 'application/octet-binary'} );
		saveAs( blob, filename );
	}
}


function identifyVoxels() {
	console.log( '%cidentify mask voxels', 'color: green;' );
	var dims = nifti1.getDims();
	nifti2 = new Nifti();
	nifti2.makeNew( dims.nx, dims.ny, dims.nz, dims.dx, dims.dy, dims.dz, 1, 'UINT8' );
	
	var bvoxelId;	//index of border voxel
	var bvoxel; 	//coordinates
	
	for( var x = 1; x < (dims.nx - 1); ++x ) {
		for( var y = 1; y < (dims.ny -1); ++y ) {
			for( var z = 1; z < (dims.nz - 1); ++z) {
				var val = nifti1.getValue( x , y, z );
				if( val > 0 ) {
					nifti2.setValue( x, y, z, 1 );
				}
				else {
					nifti2.setValue( x, y, z, 0 );
				}
			}
		}
	}
	for( var x = 1; x < (dims.nx - 1); ++x ) {
		for( var y = 1; y < (dims.ny -1); ++y ) {
			top:
			for( var z = 1; z < (dims.nz - 1); ++z) {
				var val = nifti2.getValue( x , y, z );	
				
				if( val == 1 ) {
					for( var i = -1; i <= 1; ++i ) {
						for( var j = -1; j<= 1; ++j ) {
							for( var k= -1; k<= 1; ++k ) {
								if( nifti2.getValue( x+i, y+j, z+k ) == 0) {
									nifti2.setValue( x, y, z, bv );
									bvoxelId = nifti1.getId( x, y , z ); 
									bvox.push( { x:x, y:y, z:z, index:bvoxelId });
									continue top;
								}
							}
						}
					}
				}
			}
		} 
	}
	//if( debug > 2 ) { console.log( '%cborder voxels (index,x,y,z): '+ JSON.stringify(bvox), 'color: green;' ); }
	if( debug > 2 ) { console.log( '%cborder voxels (index x,y,z): ', 'color: orange;' ); console.log( bvox ); }
	if( debug > 2 ) { console.log( 'nifti2: ', nifti2 ); }
}



function genVectors() {
	console.log( '%cgenerate vectors', 'color: green;' );

	var dims = nifti2.getDims();
	var blocksize = dims.nx * dims.ny * dims.nz;
	
	//nifti3 to store 1 for every visited voxel
	nifti3 = new Nifti();
	nifti3.makeNew( dims.nx, dims.ny, dims.nz, dims.dx, dims.dy, dims.dz, 1, 'UINT8' );
	
	//nifti4 to store direction vector for every voxel
	nifti4 = new Nifti();
	nifti4.makeNew( dims.nx, dims.ny, dims.nz, dims.dx, dims.dy, dims.dz, 3, 'FLOAT32' );
	
	
	//empty box of size blocksize * nd
	//-->to get indices: k * dims.nx * dims.ny * dims.nz + z * dims.nx * dims.ny + y * dims.nx + x
	vol = new Float32Array( blocksize * nd );

	//'real' position inside one voxel
	var xPos;
	var yPos;
	var zPos;
	
	for( var i = 0; i < np; ++i ) {
		//take random surface voxel to start fiber
		var rindex = Math.floor( Math.random() * ( bvox.length - 1 ) )
		if( debug > 2 ) { console.log( '%crandom index: ', 'color: orange' ); console.log( rindex ); }
		
		if( debug > 2 ) { console.log( '%cborder voxels: index and coordinates: ', 'color: orange' ); console.log( bvox[ rindex ] ); }


		var voxel = { x: bvox[ rindex ].x, y: bvox[ rindex ].y, z: bvox[ rindex ].z };
		//var voxV; // voxel vector: to keep direction of border voxel at its index bvox[rindex] at position [0]
		//v1 is the border voxel direction, so
		var voxelId = bvox[ rindex ].index;
		
		xPos = voxel.x;
		yPos = voxel.y;
		zPos = voxel.z;
		
		//generate random direction for border voxel (random vector v)
		var v = { x:random(), y:random(), z:random() };
		if( debug > 2 ) { console.log( 'v = ', v ); }

var count = 0;
		while( nifti2.getValue( voxel.x, voxel.y, voxel.z ) > 0 ) {
	
		
			//var voxV (voxelVector) to keep the direction of the border voxel
			//var voxV = 
	
			//generate new random direction (new random vector v2)
			v2 = { x:random(), y:random(), z:random() };
			if( debug > 2 ) { console.log( 'v2 = ', v2 ); }
	
			//(new) direction 
			v = { x: w * v.x + (1-w)* v2.x , y: w * v.y + (1-w)* v2.y, z: w * v.z + (1-w)* v2.z };
			if( debug > 2 ) { console.log( 'new v = ', v ); }
	
			
			//set every visited voxel to 1
			//nifti3.setValue( voxel.x, voxel.y, voxel.z, [ 1, 1, 1 ] );
			nifti3.setValue( voxel.x, voxel.y, voxel.z, 1 );
			
			//write direction into border voxel
			nifti4.setValue( voxel.x, voxel.y, voxel.z, [ v.x, v.y, v.z ] );
			
			//[ bvox[ rindex ].index ] = { x:v1.x, y:v1.y, z:v1.z };
			//new voxel in next round nifti3.setValue( 
			//nifti3[ newVoxelInd ] = { x:v.x, y:v.y, z:v.z };



			//normalise our direction vector? ((to not stay like 'infinitely' inside one and the same voxel))
			//length of vector: lv
			var lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z ); 
			v = { x: v.x / lv, y: v.y / lv, z: v.z / lv };
			if( debug > 2 ) { console.log( '%cnormalised vector: ', 'color: orange' ); console.log( lv ); }			
			
			//go to the new voxel (which is index in 3*1dim box)
			xPos += v.x;
			yPos += v.y;
			zPos += v.z;

			//make int to get the index
			var newVoxelInt = { x: Math.floor( xPos ), y: Math.floor( yPos ), z: Math.floor( zPos ) };
			if( debug > 2 ) { console.log( 'new voxel int is at xyz ', newVoxelInt ); }

			//get index of new voxel in databox
			var newVoxelInd = nifti2.getId( newVoxelInt.x, newVoxelInt.y, newVoxelInt.z );
			if( debug > 2 ) { console.log( '%cnew voxel index ', 'color: orange;' ); console.log( newVoxelInd ); }
		
		
			
			voxel = newVoxelInt;
			voxelId = newVoxelInd;
			
			console.log( 'loop ' + count++ );

		}
	
		
		/*
		outcommented to first get the fiber running over the block, 
		we might end up not needing that at all 
		//get value for this voxel at newVoxelInd in k directions
		
		for( i = 0; i < nd; ++i ) {
		
			//get dirV accessible
			var dirV = { x: dir[i].x, y: dir[i].y, z: dir[i].z };
			//calculate dot product vol[currentVoxel] (v.x v.y and v.z) with dir[0++]
			//var val = dot( [v.x, v.y, v.z], [dir[i].x, dir[i].y, dir[i].z] );
			//var val = dot( [v.x, v.y, v.z], [dirV.x, dirV.y, dirV.z] );
			var val = dot( v, dirV );
			vol[newVoxelInd + i * blocksize] = val;
			
			if( debug > 2 ) { console.log( '%ccurrentVoxelIndex ', 'color: orange' ); console.log( newVoxelInd + i * blocksize ); }
			if( debug > 2 ) { console.log( '%cnew volume arry ', 'color: orange;' ); console.log( vol ); }
		}

		vvC[newVoxelInd] = vvC[newVoxelInd]?(vvC[newVoxelInd]+1):1;
		if( debug > 2 ) { console.log( '%cvoxel visit count ', 'color: orange;' ); console.log( vvC );
		*/
	}
	if( debug > 2 ) { console.log( '%cgenerate voxels finito', 'color: light-green;' ); }
	
	for( var x = 1; x < (dims.nx - 1); ++x ) {
		for( var y = 1; y < (dims.ny -1); ++y ) {
			for( var z = 1; z < (dims.nz - 1); ++z) {
				var val = nifti3.getValue( x , y, z );
				if( val != 1 ) {
					nifti3.setValue( x, y, z, 0 );
				}
				else {
					nifti3.setValue( x, y, z, 200 );
				}
			}
		}
	}
	if( debug > 2 ) { console.log( 'nifti3 after loop: ', nifti3 ); }
}
	
	





function random() {
	return ( Math.random() - 0.5 );
}



function dot( a, b ) {
	return( a.x * b.x + a.y * b.y + a.z * b.z )
}


/*
function dot( a, b ) { 
	var result = 0; 
	for( var i = 0; i < a.length; ++i) {
		result += a[i] * b[i];
		return result;
	}
}
*/



function getBvecsIn() {
	console.log( '%cget bvecs in', 'color: green;' );
	
	dir.push( { x:-0.26901271050312, y:0.426217878501591, z:-0.864094805089212 } );
	dir.push( { x:0.389308175473584, y:-0.464435056102534, z:-0.795882261217864 } );
	dir.push( { x:0.445081230516361, y:-0.1976420477346, z:0.873801848108661 } );
	//ToDO: read file from disc genre { bvec[i].x, bvec[i].y, bvec[i].z };
	
	if( debug > 2 ) { console.log( '%cdirections dir (x,y,z): ', 'color: orange;' ); console.log( 'directions: ', dir ); }
}





/* visitation map /*
/*
	as volume vector, initialized with aaaall zeroes öö
	--> count up voxel entry of 21536 by one
	--> visitation map of 21536++
	
	var dims = nifti2.getDims();
	var blocksize = dims.nx * dims.ny * dims.nz;
	
	vvC[newVoxelInd] = vvC[newVoxelInd]?(vvC[newVoxelInd]+1):1;
		if( debug > 2 ) { console.log( '%cvoxel visit count ', 'color: orange;' ); console.log( vvC );
*/












/*
some code to keep for more than one line at once?

function line(gx0, gy0, gz0, gx1, gy1, gz1) {

    var gx0idx = Math.floor(gx0);
    var gy0idx = Math.floor(gy0);
    var gz0idx = Math.floor(gz0);

    var gx1idx = Math.floor(gx1);
    var gy1idx = Math.floor(gy1);
    var gz1idx = Math.floor(gz1);
	//zweiter wert groesser als das erste
	//dann ist sx 1, ansonsten guckt er, ob der 
	//richtungsvektor, der in richtung wie mein vektor zeigt, unit vector erzeugt hier
    var sx = gx1idx > gx0idx ? 1 : gx1idx < gx0idx ? -1 : 0;
    var sy = gy1idx > gy0idx ? 1 : gy1idx < gy0idx ? -1 : 0;
    var sz = gz1idx > gz0idx ? 1 : gz1idx < gz0idx ? -1 : 0;

    var gx = gx0idx;
    var gy = gy0idx;
    var gz = gz0idx;

    //Planes for each axis that we will next cross
    var gxp = gx0idx + (gx1idx > gx0idx ? 1 : 0);
    var gyp = gy0idx + (gy1idx > gy0idx ? 1 : 0);
    var gzp = gz0idx + (gz1idx > gz0idx ? 1 : 0);

    //Only used for multiplying up the error margins
    var vx = gx1 === gx0 ? 1 : gx1 - gx0;
    var vy = gy1 === gy0 ? 1 : gy1 - gy0;
    var vz = gz1 === gz0 ? 1 : gz1 - gz0;

    //Error is normalized to vx * vy * vz so we only have to multiply up
    var vxvy = vx * vy;
    var vxvz = vx * vz;
    var vyvz = vy * vz;

    //Error from the next plane accumulators, scaled up by vx*vy*vz
    // gx0 + vx * rx === gxp
    // vx * rx === gxp - gx0
    // rx === (gxp - gx0) / vx
    var errx = (gxp - gx0) * vyvz;
    var erry = (gyp - gy0) * vxvz;
    var errz = (gzp - gz0) * vxvy;

    var derrx = sx * vyvz;
    var derry = sy * vxvz;
    var derrz = sz * vxvy;

    do {
        vol[gz*hdr.dim[1]*hdr.dim[2]+gy*hdr.dim[1]+gx];
        if (gx === gx1idx && gy === gy1idx && gz === gz1idx)
        	break;

        //Which plane do we cross first?
        var xr = Math.abs(errx);
        var yr = Math.abs(erry);
        var zr = Math.abs(errz);

        if (sx !== 0 && (sy === 0 || xr < yr) && (sz === 0 || xr < zr)) {
            gx += sx;
            errx += derrx;
        }
        else if (sy !== 0 && (sz === 0 || yr < zr)) {
            gy += sy;
            erry += derry;
        }
        else if (sz !== 0) {
            gz += sz;
            errz += derrz;
        }

    } while (true);
}
*/






$(document).ready(function() {
	document.getElementById('fiMask').addEventListener('change', loadMask, false);
	document.getElementById('btSurfaceVoxels').addEventListener('click', identifyVoxels, false);
	//document.getElementById('btBvecs').addEventListener('click', getBvecsIn, false);
	document.getElementById('btVectors').addEventListener('click', genVectors, false);
	//document.getElementById('buttonDownload').addEventListener('click', downloadFile, false);
	document.getElementById('buttonSave').addEventListener('click', saveFile, false);
	document.getElementById('buttonInvert').addEventListener('click', invert, false);
	$('#slider1').change(function() {
		sliderCor($(this).val());
    	//console.log($(this).next().html($(this).val()));
	});
	
	$('#slider2').change(function() {
		sliderCor2($(this).val());
		//sliderSag($(this).val());
    	//console.log($(this).next().html($(this).val()));
	});
	
	$('#slider3').change(function() {
		sliderCor3($(this).val());
		//sliderAxi($(this).val());
    	//console.log($(this).next().html($(this).val()));
	});
});



	/*##################*/
	
	/* older code to generate the vector and then compare how much the random directions would match one of the real measured directions */	
	/*
	//for( k=0; k <= 3; ++k ) {
	
	for( var x = 0; x < dims.nx; ++x ) {
		for( var y = 0; y < dims.ny; ++y ) {
			for( var z = 0; z < dims.nz; ++z ) {
				var val = nifti2.getValue( x , y, z );
				if( val == bv ) {
					makeVector( x, y, z, w );
				}
			}
		}
	}
}
*/


/*
function makeVector( x, y, z ) {
	//generate direction for border voxel
		v1 = { x:random(), y:random(), z:random() };
		if( debug > 2 ) { console.log( 'v1 = ', v1 ); }
	
		//generate random direction
		v2 = { x:random(), y:random(), z:random() };
		if( debug > 2 ) { console.log( 'v2 = ', v2 ); }
	
		//new direction 
		v = { x: w * v1.x + (1-w)* v2.x , y: w * v1.y + (1-w)* v2.y, z: w * v1.z + (1-w)* v2.z };
		if( debug > 2 ) { console.log( 'v = ', v ); }
	
		//go to the new voxel (which is index in 3*1dim box)
		var newVoxel = { x: bvox[ rindex ].x + v.x, y: bvox[ rindex ].y + v.y, z: bvox[ rindex ].z + v.z };
		if( debug > 2 ) { console.log( '%cnew voxel at xyz ', 'color: orange;' ); console.log( newVoxel ); }


		//make int to get the index
		var newVoxelInt = { x: Math.floor( bvox[ rindex ].x + v.x ), y: Math.floor( bvox[ rindex ].y + v.y ), z: Math.floor( bvox[ rindex ].z + v.z ) };
		if( debug > 2 ) { console.log( 'new voxel int is at xyz ', newVoxelInt ); }

		//get index of new voxel in databox
		var newVoxelInd = nifti2.getId( newVoxelInt.x, newVoxelInt.y, newVoxelInt.z );
		if( debug > 2 ) { console.log( '%cnew voxel index ', 'color: orange;' ); console.log( newVoxelInd ); }
		
		var nifti3 = new Nifti();
		
		

	
		
		//get value for this voxel at newVoxelInd in k directions
		for( i = 0; i < nd; ++i ) {
		
			//get dirV accessible
			var dirV = { x: dir[i].x, y: dir[i].y, z: dir[i].z };
			//calculate dot product vol[currentVoxel] with dir[0++]
			//var val = dot( [v.x, v.y, v.z], [dir[i].x, dir[i].y, dir[i].z] );
			//var val = dot( [v.x, v.y, v.z], [dirV.x, dirV.y, dirV.z] );
			var val = dot( v, dirV );
			vol[newVoxelInd + i * blocksize] = val;
			
			if( debug > 2 ) { console.log( '%cnew volume arry ', 'color: orange;' ); console.log( vol ); }
		}
}

*/


