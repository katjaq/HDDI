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
