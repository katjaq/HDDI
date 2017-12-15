var HDDIRandom = {
    /**
     * @desc https://en.wikipedia.org/wiki/Marsaglia_polar_method
     */
    nrandomPrecomputed: 0,
    nrandomIsPrecomputedReady: false,
    nrandom: function nrandom() {
        if(this.nrandomIsPrecomputedReady) {
            this.nrandomIsPrecomputedReady = false;

            return this.nrandomPrecomputed;
        } else {
            let x1, x2, s, y1;
            do {
                x1 = 2*Math.random() - 1;
                x2 = 2*Math.random() - 1;
                s = x1 * x1 + x2 * x2;
            } while ( s >= 1 || s === 0 );
            s = Math.sqrt( (-2*Math.log(s) ) / s );
            y1 = x1 * s;
            this.nrandomPrecomputed = x2 * s;
            this.nrandomIsPrecomputedReady = true;

            return y1;
        }
    },

    randomDirectionNaive: function randomDirection() {
        let v = {
            x: Math.random() - 0.5,
            y: Math.random() - 0.5,
            z: Math.random() - 0.5
        };
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: v.x/lv,
            y: v.y/lv,
            z: v.z/lv
        };
    },

    /**
      * @desc Muller's method from http://mathworld.wolfram.com/SpherePointPicking.html
      */
    randomDirectionMuller: function randomDirection() {
        let v = {
            x: this.nrandom(),
            y: this.nrandom(),
            z: this.nrandom()
        };
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: v.x/lv,
            y: v.y/lv,
            z: v.z/lv
        };
    },

    /**
      * @desc Muller's method adapted to an ellipsoid. Goes back to the case of a sphere if a=b=c
      */
    randomDirectionMullerEllipsoid: function randomDirectionMullerEllipsoid(a, b, c) {
        let v = {
            x: a * this.nrandom(),
            y: b * this.nrandom(),
            z: c * this.nrandom()
        };
        let lv = Math.sqrt( v.x*v.x + v.y*v.y + v.z*v.z );
        return {
            x: v.x/lv,
            y: v.y/lv,
            z: v.z/lv
        };
    },

    /**
      * @desc Marsaglia's method from http://mathworld.wolfram.com/SpherePointPicking.html
      */
    randomDirection: function randomDirection() {
        let x1, x2, s;
        do {
            x1 = 2*Math.random() - 1;
            x2 = 2*Math.random() - 1;
            s = x1 * x1 + x2 * x2;
        } while ( s >= 1 || s === 0);
        return v = {
            x: 2 * x1 * Math.sqrt(1 - s),
            y: 2 * x2 * Math.sqrt(1 - s),
            z: 1 - 2 * s
        };
    }
}
