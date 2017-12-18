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

