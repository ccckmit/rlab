// http://cs.stanford.edu/people/karpathy/convnetjs/
// https://github.com/junku901/dnn
var NN = R.NN = module.exports = {}

NN.sigmoid = function(x) {
	var sigmoid = (1. / (1 + Math.exp(-x)))
	if(sigmoid ==1) {
	 //   console.warn("Something Wrong!! Sigmoid Function returns 1. Probably javascript float precision problem?\nSlightly Controlled value to 1 - 1e-14")
			sigmoid = 0.99999999999999; // Javascript Float Precision Problem.. This is a limit of javascript.
	} else if(sigmoid ==0) {
		//  console.warn("Something Wrong!! Sigmoid Function returns 0. Probably javascript float precision problem?\nSlightly Controlled value to 1e-14")
			sigmoid = 1e-14;
	}
	return sigmoid; // sigmoid cannot be 0 or 1;
}

NN.dSigmoid = function(x){
	a = NN.sigmoid(x);
	return a * (1. - a);
}

NN.softmaxVec = function(vec) {
    var max = vec.max();
    var preSoftmaxVec = vec.map((x)=>Math.exp(x - max));
    var preSoftmaxSum = preSoftmaxVec.sum();
		return preSoftmaxVec.map((x)=>(x/preSoftmaxSum));
}

NN.softmaxMat = function(m) {
    var len=m.length, r=new Array();
    for(var i=0; i<len; i++)
        r[i]=NN.softmaxVec(mat[i]);
    return r;
}

NN.binarySample = function(m) {
	return m.map1((x)=>(Math.random()<x)?1:0);
}

NN.RBM = require("./neural/rbm");
NN.MLP = require("./neural/mlp");
NN.HiddenLayer = require("./neural/hiddenLayer");
NN.DBN = require("./neural/dbn");
NN.CRBM = require("./neural/crbm");
NN.CDBN = require("./neural/cdbn");

