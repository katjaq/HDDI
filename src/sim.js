let HDDISim = (function HDDISim() {
"use strict";
    var me = {
        getValue: function getValue(vol, dim, x, y, z) {
            return vol[z*dim[1]*dim[0] + y*dim[0] + x];
        },

        getValueN: function getValue(vol, dim, x, y, z) {
            let i, val = [];
            for(i=0; i<HDDI.params.dir.length; i++) {
                val.push(vol[i][z*dim[1]*dim[0] + y*dim[0] + x]);
            }

            return val;
        },

        setValue: function setValue(vol,dim, x, y, z, val) {
            vol[z*dim[1]*dim[0] + y*dim[0] + x] = val;
        },

        setValueN: function setValue(vol,dim, x, y, z, val) {
            let i;
            for(i=0; i<HDDI.params.dir.length; i++) {
                vol[i][z*dim[1]*dim[0] + y*dim[0] + x] = Math.abs(me.dot(HDDI.params.dir[i], val));
            }
        },

        addValue: function setValue(vol,dim, x, y, z, val) {
            vol[z*dim[1]*dim[0] + y*dim[0] + x] += val;
        },

        addValueN: function setValue(vol,dim, x, y, z, val) {
            let i;
            for(i=0; i<HDDI.params.dir.length; i++) {
                vol[i][z*dim[1]*dim[0] + y*dim[0] + x] += Math.abs(me.dot(HDDI.params.dir[i], val));
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

        /**
          * @func identifyVoxels
          * @desc Identify voxels as background (0), surface (bv) or core (1)
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

            // initialise 26-connected neighbourhoud
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
                        val = me.getValue( vol, dim, x , y, z );
                        if( val > 0 ) {
                            me.setValue( idvol, dim, x, y, z, 1 );
                        }
                        else {
                            me.setValue( idvol, dim, x, y, z, 0 );
                        }
                    }
                }
            }

            // mark "surface" values
            for( x = 1; x < dim[0]-1; ++x ) {
                for( y = 1; y < dim[1]-1; ++y ) {
                    for( z = 1; z < dim[2]-1; ++z) {
                        val = me.getValue( vol, dim, x , y, z );
                        if( val === 0 ) {
                            continue;
                        }
                        for( l=0; l<con26.length; l++ ) {
                            [i, j, k] = con26[l];
                            if( me.getValue( vol, dim, x+i, y+j, z+k ) === 0 ) {
                                me.setValue( idvol, dim, x, y, z, HDDI.bv );
                                break;
                            }
                        }
                    }
                }
            }

            // store coordinates of border voxels
            for( x = 0; x < dim[0]; ++x ) {
                for( y = 0; y < dim[1]; ++y ) {
                    for( z = 0; z < dim[2]; ++z) {
                        if(me.getValue(idvol, dim, x, y, z) === HDDI.bv ) {
                            HDDI.bvox.push({ x: x, y: y, z: z });
                        }
                    }
                }
            }

            if( HDDI.debug > 2 ) {
                console.log( '%cborder voxels (index x, y, z): ', 'color: orange; ' );
                console.log( HDDI.bvox );
            }

            return idvol;
        },

        randomDirection: function randomDirection() {
            let v = {
                x: me.random(), // green
                y: me.random(), // blue
                z: me.random()  // red
/*
                x: 0, //me.random(), // green
                y: 0, //me.random(), // blue
                z: 1, //me.random()  // red
*/
            };

            return v;
        },

        /**
          * @func genVectors
          * @desc Generates random fibres
          */
        /**
          * @todo Rename to genStreamlines
          */
        genStreamlines: function genStreamlines(vol, dim) {
            console.log( '%cgenerate vectors', 'color: green; ' );

            const sz = dim[0]*dim[1]*dim[2];
            let rho = new Float32Array(sz);
            let dir = [];
            let i, j;
            let count;
            let pos, rindex;
            let ix, iy, iz;
            let length, max, min;
            let v, v2;
            let bar = new progress.Bar({}, progress.Presets.shades_classic);

            bar.start(HDDI.params.ns, 0);

            // initialise rho and dir volumes to 0
            for(i=0;i<sz;i++) {
                rho[i] = 0;
            }
            for(i=0;i<HDDI.params.dir.length;i++) {
                dir[i] = new Float32Array(sz);
                for(j=0;j<sz;j++) {
                    dir[i][j] = 0;
                }
            }

            // generate fibres
            count = 0;
            while( count < HDDI.params.ns) {
                let fib = [];

                if( count%1000 === 0 ) {
                    bar.update(count);
                }

                //take random surface voxel to start fiber
                rindex = Math.round( Math.random() * ( HDDI.bvox.length - 1 ) );
                pos = {
                    x: HDDI.bvox[rindex].x,
                    y: HDDI.bvox[rindex].y,
                    z: HDDI.bvox[rindex].z
                };
                v = me.normalise(me.randomDirection(), HDDI.params.step);
                fib.push([pos, v]);
                length = 0;
                while( me.getValue( vol, dim, parseInt(pos.x), parseInt(pos.y), parseInt(pos.z) ) > 0 ) {
                    // random direction
                    v2 = me.normalise(me.randomDirection(), HDDI.params.step);
                    //new direction = mix of old plus random
                    v = me.normalise({
                        x: HDDI.params.w * v.x + (1-HDDI.params.w)* v2.x ,
                        y: HDDI.params.w * v.y + (1-HDDI.params.w)* v2.y,
                        z: HDDI.params.w * v.z + (1-HDDI.params.w)* v2.z
                    }, HDDI.params.step);
                    
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
                if(length < HDDI.params.minFibLength ) {
                    continue;
                }

                fs.writeFileSync('fib.json', JSON.stringify(fib.map((o) => [o[0].x, o[0].y, o[0].z].join(' '))));

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
                    me.addValue( rho, dim, ix, iy, iz, 1 );
                    me.addValueN( dir, dim, ix, iy, iz, v );
                }

                count += 1;
            }
            bar.update(HDDI.params.ns);
            bar.stop();
            if( HDDI.debug > 2 ) {
                console.log( '%cgenerate voxels finito', 'color: light-green; ' );
            }
            console.log("Min and Max fibre length:", min, max);

            return {
                rho: rho,
                dir: dir
            };
        },

        computeB0: function computeB0(vol, dim) {
            let i, j, k;
            let b0 = new Float32Array(dim[0] * dim[1] *dim[2]);
            for(i=0;i<dim[0];i++) {
                for(j=0;j<dim[1];j++) {
                    for(k=0;k<dim[2];k++) {
                        let sig = HDDI.getValueN(vol, dim, i, j, k);
                        HDDI.setValue(b0, dim, i, j, k, Math.max(...sig));
                    }
                }
            }

            return b0;
        },

        firstDirection: function firstDirection(vol, dim) {
            let i, j, sz = dim[0] * dim[1] * dim[2];
            let max;
            let dif, evals, evecs;
            let first = [new Float32Array(sz), new Float32Array(sz), new Float32Array(sz)];
            const iH = invHfromG(HDDI.params.dir);
            for(i=0; i<sz; i++) {
                let S = [0];
                for(j=0; j<HDDI.params.dir.length; j++) {
                    S[1+j] = vol[j][i];
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

        random: function random() {
            return ( Math.random() - 0.5 );
        },

        dot: function dot( a, b ) {
            return( a.x * b.x + a.y * b.y + a.z * b.z )
        },

        getBvecsIn: function getBvecsIn() {
            console.log( '%cget bvecs in', 'color: green; ' );

            HDDI.dir.push( { x:-0.26901271050312, y:0.426217878501591, z:-0.864094805089212 } );
            HDDI.dir.push( { x:0.389308175473584, y:-0.464435056102534, z:-0.795882261217864 } );
            HDDI.dir.push( { x:0.445081230516361, y:-0.1976420477346, z:0.873801848108661 } );
            // @todo read file from disc, genre { bvec[i].x, bvec[i].y, bvec[i].z };

            if( HDDI.debug > 2 ) {
                console.log( '%cdirections dir (x, y, z): ', 'color: orange; ' );
                console.log( 'directions: ', HDDI.dir );
            }
        }
    };

    return me;
}());
