var R = require('../rlab');

var RBM = R.NN.RBM;

var data = [
[1,1,1,0,0,0],
[1,0,1,0,0,0],
[1,1,1,0,0,0],
[0,0,1,1,1,0],
[0,0,1,1,0,0],
[0,0,1,1,1,0]];

var rbm = new RBM(R, {
    input : data,
    nVisible : 6,
    nHidden : 2
});

rbm.train({
    lr : 0.6,
    k : 1,
    epochs : 500
});

var v = [[1, 1, 0, 0, 0, 0],
    [0, 0, 0, 1, 1, 0]];

console.log(rbm.reconstruct(v).strM());
console.log(rbm.sampleHgivenV(v)[0].strM());