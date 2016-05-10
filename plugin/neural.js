// http://cs.stanford.edu/people/karpathy/convnetjs/
// https://github.com/junku901/dnn
var NN = {};

NN.sigmoid = function(x) {
	var sigmoid = (1. / (1 + Math.exp(-x)))
	if(sigmoid ==1) {
	 //   console.warn("Something Wrong!! Sigmoid Function returns 1. Probably javascript float precision problem?\nSlightly Controlled value to 1 - 1e-14")
			sigmoid = 0.99999999999999; // Javascript Float Precision Problem.. This is a limit of javascript.
	} else if(sigmoid ==0) {
		//  console.warn("Something Wrong!! Sigmoid Function returns 0. Probably javascript float precision problem?\nSlightly Controlled value to 1e-14")
			sigmoid = 1e-14;
	}
	return sigmoid; // sigmoid cannot be 0 or 1;;
}

NN.dSigmoid = function(x){
	a = B.sigmoid(x);
	return a * (1. - a);
}

module.exports = NN;