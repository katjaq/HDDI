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

