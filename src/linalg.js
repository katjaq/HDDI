var HDDILinAlg = {
    /**
      * @desc 3D vector addition
      */
    add: function add( a, b ) {
        return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z };
    },

    /**
      * @desc 3D vector subtraction
      */
    sub: function sub( a, b ) {
        return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
    },

    /**
      * @desc 3D vector dot product
      */
    dot: function dot( a, b ) {
        return( a.x * b.x + a.y * b.y + a.z * b.z )
    },

    normalise: function normalise(v, l) {
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: l * v.x / lv,
            y: l * v.y / lv,
            z: l * v.z / lv
        };
    },

    normVec: function normVec(v) {
        let lv = Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
        return {
            x: v.x / lv,
            y: v.y / lv,
            z: v.z / lv
        };
    },

    norm: function norm(v) {
        return Math.sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
    },

    scale: function scale(v, l) {
        return {
            x: l * v.x,
            y: l * v.y,
            z: l * v.z
        };
    },

    /**
      * @desc multiplication of matrices a and b
      */
    mult: function mult(a, b) {
        let r,c,i;
        let res = new Array(a.length).fill(0);
        res=res.map(o=>new Array(b[0].length).fill(0));
        for(r=0;r<a.length;r++) {
            for(c=0;c<b[0].length;c++) {
                for(i=0;i<a[0].length;i++) {
                    res[r][c] += a[r][i]*b[i][c];
                }
            }
        }

        return res;
    },

    /**
      * @desc transpose of matrix a
      */
    transpose: function transpose(a) {
        let i,j,r = new Array(a[0].length).fill(0);
        r=r.map(o=>new Array(a.length).fill(0));
        for(i=0;i<a.length;i++)
        for(j=0;j<a[0].length;j++) {
            r[j][i] = a[i][j];
        }
        return r;
    },

    /**
      * @desc make a diagonal matrix from vector v
      */
    diag: function diag(v) {
        let m = new Array(v.length).fill(0);
        m = m.map((o) => new Array(v.length).fill(0));
        v.map((o,i)=>m[i][i]=o);

        return m;
    },

    /**
      * @desc This function comes from the package numeric.js 
      *       Shanti Rao sent me this routine by private email. I had to modify it
      *       slightly to work on Arrays instead of using a Matrix object.
      *       It is apparently translated from http://stitchpanorama.sourceforge.net/Python/svd.py
      */

    svd: function svd(A) {
        //Compute the thin SVD from G. H. Golub and C. Reinsch, Numer. Math. 14, 403-420 (1970)
        let temp;
        let prec= 2.220446049250313e-16; //numeric.epsilon; //Math.pow(2,-52) // assumes double prec
        let tolerance= 1.e-64/prec;
        let itmax= 50;
        let c=0;
        let i=0;
        let j=0;
        let k=0;
        let l=0;

        let u= JSON.parse(JSON.stringify(A)); // numeric.clone(A);
        let m= u.length;
        let n= u[0].length;

        if (m < n) {
            throw "Need more rows than columns";
        }

        let e = new Array(n);
        let q = new Array(n);
        for (i=0; i<n; i++) {
            e[i] = 0.0;
            q[i] = 0.0;
        }
        let v = new Array(n).fill(0);
        v = v.map((o) => new Array(n).fill(0));

        function pythag(a,b) {
            a = Math.abs(a)
            b = Math.abs(b)
            if (a > b)
                return a*Math.sqrt(1.0+(b*b/a/a))
            else if (b == 0.0)
                return a
            return b*Math.sqrt(1.0+(a*a/b/b))
        }

        //Householder's reduction to bidiagonal form

        let f= 0.0;
        let g= 0.0;
        let h= 0.0;
        let x= 0.0;
        let y= 0.0;
        let z= 0.0;
        let s= 0.0;

        for (i=0; i < n; i++) {
            e[i]= g;
            s= 0.0;
            l= i+1;
            for (j=i; j < m; j++)
                s += (u[j][i]*u[j][i]);
            if (s <= tolerance) {
                g= 0.0;
            } else {
                f= u[i][i];
                g= Math.sqrt(s);
                if (f >= 0.0) g= -g;
                h= f*g-s
                u[i][i]=f-g;
                for (j=l; j < n; j++)
                {
                    s= 0.0
                    for (k=i; k < m; k++)
                        s += u[k][i]*u[k][j]
                    f= s/h
                    for (k=i; k < m; k++)
                        u[k][j]+=f*u[k][i]
                }
            }
            q[i]= g
            s= 0.0
            for (j=l; j < n; j++) {
                s= s + u[i][j]*u[i][j]
            }
            if (s <= tolerance) {
                g= 0.0;
            } else {
                f= u[i][i+1]
                g= Math.sqrt(s)
                if (f >= 0.0) g= -g
                h= f*g - s
                u[i][i+1] = f-g;
                for (j=l; j < n; j++) {
                    e[j]= u[i][j]/h;
                }
                for (j=l; j < m; j++) {
                    s=0.0
                    for (k=l; k < n; k++) {
                        s += (u[j][k]*u[i][k])
                    }
                    for (k=l; k < n; k++) {
                        u[j][k]+=s*e[k]
                    }
                }
            }
            y= Math.abs(q[i])+Math.abs(e[i])
            if (y>x)
                x=y
        }

        // accumulation of right hand transformations
        for (i=n-1; i != -1; i+= -1) {
            if (g != 0.0) {
                h= g*u[i][i+1];
                for (j=l; j < n; j++) {
                    v[j][i]=u[i][j]/h;
                }
                for (j=l; j < n; j++) {
                    s=0.0
                    for (k=l; k < n; k++) {
                        s += u[i][k]*v[k][j];
                    }
                    for (k=l; k < n; k++) {
                        v[k][j]+=(s*v[k][i]);
                    }
                }
            }
            for (j=l; j < n; j++) {
                v[i][j] = 0;
                v[j][i] = 0;
            }
            v[i][i] = 1;
            g= e[i]
            l= i
        }
        // accumulation of left hand transformations
        for (i=n-1; i != -1; i+= -1) {
            l= i+1
            g= q[i]
            for (j=l; j < n; j++)
                u[i][j] = 0;
            if (g != 0.0) {
                h= u[i][i]*g
                for (j=l; j < n; j++) {
                    s=0.0
                    for (k=l; k < m; k++) s += u[k][i]*u[k][j];
                    f= s/h
                    for (k=i; k < m; k++) u[k][j]+=f*u[k][i];
                }
                for (j=i; j < m; j++) u[j][i] = u[j][i]/g;
            }
            else
                for (j=i; j < m; j++) u[j][i] = 0;
            u[i][i] += 1;
        }

        // diagonalization of the bidiagonal form
        prec= prec*x
        for (k=n-1; k != -1; k+= -1) {
            for (let iteration=0; iteration < itmax; iteration++)
            { // test f splitting
                let test_convergence = false
                for (l=k; l != -1; l+= -1)
                {
                    if (Math.abs(e[l]) <= prec) {
                        test_convergence= true;
                        break;
                    }
                    if (Math.abs(q[l-1]) <= prec) {
                        break
                    }
                }
                if (!test_convergence)
                { // cancellation of e[l] if l>0
                    c= 0.0
                    s= 1.0
                    let l1= l-1
                    for (i =l; i<k+1; i++)
                    {
                        f= s*e[i]
                        e[i]= c*e[i]
                        if (Math.abs(f) <= prec)
                            break
                        g= q[i]
                        h= pythag(f,g)
                        q[i]= h
                        c= g/h
                        s= -f/h
                        for (j=0; j < m; j++)
                        {
                            y= u[j][l1]
                            z= u[j][i]
                            u[j][l1] =  y*c+(z*s)
                            u[j][i] = -y*s+(z*c)
                        }
                    }
                }
                // test f convergence
                z= q[k]
                if (l == k)
                { //convergence
                    if (z<0.0)
                    { //q[k] is made non-negative
                        q[k]= -z
                        for (j=0; j < n; j++)
                            v[j][k] = -v[j][k]
                    }
                    break  //break out of iteration loop and move on to next k value
                }
                if (iteration >= itmax-1) {
                    throw 'Error: no convergence.'
                }
                // shift from bottom 2x2 minor
                x= q[l]
                y= q[k-1]
                g= e[k-1]
                h= e[k]
                f= ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
                g= pythag(f,1.0)
                if (f < 0.0)
                    f= ((x-z)*(x+z)+h*(y/(f-g)-h))/x
                else
                    f= ((x-z)*(x+z)+h*(y/(f+g)-h))/x
                // next QR transformation
                c= 1.0
                s= 1.0
                for (i=l+1; i< k+1; i++) {
                    g= e[i]
                    y= q[i]
                    h= s*g
                    g= c*g
                    z= pythag(f,h)
                    e[i-1]= z
                    c= f/z
                    s= h/z
                    f= x*c+g*s
                    g= -x*s+g*c
                    h= y*s
                    y= y*c
                    for (j=0; j < n; j++)
                    {
                        x= v[j][i-1]
                        z= v[j][i]
                        v[j][i-1] = x*c+z*s
                        v[j][i] = -x*s+z*c
                    }
                    z= pythag(f,h)
                    q[i-1]= z
                    c= f/z
                    s= h/z
                    f= c*g+s*y
                    x= -s*g+c*y
                    for (j=0; j < m; j++) {
                        y= u[j][i-1]
                        z= u[j][i]
                        u[j][i-1] = y*c+z*s
                        u[j][i] = -y*s+z*c
                    }
                }
                e[l]= 0.0
                e[k]= f
                q[k]= x
            }
        }
        //vt= transpose(v)
        //return (u,q,vt)
        for (i=0;i<q.length; i++)
          if (q[i] < prec) q[i] = 0

        //sort eigenvalues
        for (i=0; i< n; i++) {
        //writeln(q)
         for (j=i-1; j >= 0; j--) {
          if (q[j] < q[i]) {
        //  writeln(i,'-',j)
           c = q[j]
           q[j] = q[i]
           q[i] = c
           for(k=0;k<u.length;k++) {
             temp = u[k][i];
             u[k][i] = u[k][j];
             u[k][j] = temp;
           }
           for(k=0;k<v.length;k++) {
             temp = v[k][i];
             v[k][i] = v[k][j];
             v[k][j] = temp;
           }
        //   u.swapCols(i,j)
        //    v.swapCols(i,j)
           i = j
          }
         }
        }

        return {U:u,S:q,V:v}
    },

    /**
      * @desc The inverse of matrix a, where a=USV', is a^(-1) = VSU': https://www.cse.unr.edu/~bebis/CS791E/Notes/SVD.pdf
      */
    svdinv: function svdinv(a) {
        const {U:u,S:s,V:v} = this.svd(a);
        let i, j;
        let d = [];
        for(i=0;i<s.length;i++) {
            d[i] = v[i].map((o,j)=>o/s[j]);
        }
        return this.mult(d,this.transpose(u));
    }
}
