const dim = [5, 5, 5]
const vol = new Float32Array(dim[0]*dim[1]*dim[2]);

function generateXGradient() {
    let i, j, k;

    for(i=0; i<dim[0]; i++) {
        for(j=0; j<dim[1]; j++) {
            for(k=0; k<dim[2]; k++) {
                vol[k*dim[1]*dim[0] + j*dim[0] + i] = i;
            }
        }
    }
}

function generateYGradient() {
    let i, j, k;

    for(i=0; i<dim[0]; i++) {
        for(j=0; j<dim[1]; j++) {
            for(k=0; k<dim[2]; k++) {
                vol[k*dim[1]*dim[0] + j*dim[0] + i] = j;
            }
        }
    }
}

function generateZGradient() {
    let i, j, k;

    for(i=0; i<dim[0]; i++) {
        for(j=0; j<dim[1]; j++) {
            for(k=0; k<dim[2]; k++) {
                vol[k*dim[1]*dim[0] + j*dim[0] + i] = k;
            }
        }
    }
}

function generateXYGradient() {
    let i, j, k;

    for(i=0; i<dim[0]; i++) {
        for(j=0; j<dim[1]; j++) {
            for(k=0; k<dim[2]; k++) {
                vol[k*dim[1]*dim[0] + j*dim[0] + i] = i + j;
            }
        }
    }
}

function generateXYZGradient() {
    let i, j, k;

    for(i=0; i<dim[0]; i++) {
        for(j=0; j<dim[1]; j++) {
            for(k=0; k<dim[2]; k++) {
                vol[k*dim[1]*dim[0] + j*dim[0] + i] = i + j + k;
            }
        }
    }
}


function print(res) {
    let x, y, z;
    for(x=0;x<dim[0];x++) {
        for(y=0;y<dim[1];y++) {
            let str = ""
            for(z=0;z<dim[2];z++) {
                str += res[z*dim[1]*dim[0] + y*dim[0] + x] + " ";
            }
            console.log(str);
        }
        console.log("\n");
    }
}

console.log("x derivative of x gradient");
generateXGradient();
print(xder());

console.log("x derivative of y gradient");
generateYGradient();
print(xder());

console.log("x derivative of xy gradient");
generateXYGradient();
print(xder());

console.log("y derivative of xy gradient");
generateXYGradient();
print(yder());
