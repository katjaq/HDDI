var HDDIPCA = {
    firstPC: function firstPC(M) {
        let i, j, r = M[0];
        const niter = 20;
        for(i=0;i<niter;i++) {
            let s = {x: 0, y: 0, z: 0};
            for(j=0;j<M.length;j++) {
                s = this.add(s, this.scale(M[j], this.dot(M[j], r)));
            }
            r = this.norm(s);
        }

        return r;
    },
    nextPC: function nextPC(M, pc1) {
        let M2 = [];
        let i;
        for(i=0;i<M.length;i++) {
            M2[i] = this.sub(M[i], this.scale(pc1, this.dot(M[i],pc1)))
        }

        return this.firstPC(M2);
    }
}
