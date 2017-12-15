const fs = require('fs');
const code = fs.readFileSync('nii.js').toString();
eval(code);

loadNiftiData('../data/baboon.nii.gz')
.then((res) => {
    console.log(res);
    saveNiftiData(res.vol, res.dim, 'baboon.nii.gz')
    .then( () => console.log("done"))
    .catch( (err) => console.log(err) );
})
.catch((err) => console.log(err));
