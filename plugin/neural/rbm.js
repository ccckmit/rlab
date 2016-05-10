var NN = require("../neural");
var R, M;

var RBM = function (rlab, settings) {
  R = rlab;
  M = R.M;
  Object.assign(this, settings);
	this.settings = {
			'log level' : 1 // 0 : nothing, 1 : info, 2: warn
	};
	var a = 1. / this.nVisible;
	this.W = this.W || M.randomM(this.nVisible,this.nHidden,-a,a);
	this.hbias = this.hbias || M.newV(this.nHidden);
  this.vbias = this.vbias || M.newV(this.nVisible);
}

RBM.prototype.train = function(settings) {
	var lr = settings.lr||0.8;
	var k  = settings.k||1;
	var epochs = settings.epochs||1500;
	this.input = settings.input||this.input;

	var i,j;
	var currentProgress = 1;
	for(i=0;i<epochs;i++) {
			/* CD - k . Contrastive Divergence */
		var ph = this.sampleHgivenV(this.input);
		var phMean = ph[0], phSample = ph[1];
		var chainStart = phSample;
		var nvMeans, nvSamples, nhMeans, nhSamples;

		for(j=0 ; j<k ; j++) {
			if (j==0) {
				var gibbsVH = this.gibbsHVH(chainStart);
				nvMeans = gibbsVH[0], nvSamples = gibbsVH[1], nhMeans = gibbsVH[2], nhSamples = gibbsVH[3];
			} else {
				var gibbsVH = this.gibbsHVH(nhSamples);
				nvMeans = gibbsVH[0], nvSamples = gibbsVH[1], nhMeans = gibbsVH[2], nhSamples = gibbsVH[3];
			}
		}
    // ((input^t*phMean)-(nvSample^t*nhMeans))*1/input.length
	  var deltaW = this.input.tr().dot(phMean).sub(nvSamples.tr().dot(nhMeans)).mul(1./this.input.length);
    // deltaW = (input*phMean)-(nvSample^t * nhMeans)*1/input.length
		var deltaVbias = this.input.sub(nvSamples).colMean();
    // deltaHbias = (phSample - nhMeans).mean(row)
    var deltaHbias = phSample.sub(nhMeans).colMean();
    // W += deltaW*lr
		this.W = this.W.add(deltaW.mul(lr));
    // vbias += deltaVbias*lr
		this.vbias = this.vbias.add(deltaVbias.dot(lr));
    // hbias += deltaHbias*lr
		this.hbias = this.hbias.add(deltaHbias.dot(lr));
		if(this.settings['log level'] > 0) {
			var progress = (1.*i/epochs)*100;
			if(progress > currentProgress) {
				console.log("RBM",progress.toFixed(0),"% Completed.");
				currentProgress+=8;
			}
		}
	}
	if(this.settings['log level'] > 0)
		console.log("RBM Final Cross Entropy : ",this.getReconstructionCrossEntropy())
};

RBM.prototype.propup = function(v) {
  // sigmoid(v*W+hbias)
	return v.dot(this.W).addMV(this.hbias).mapM(NN.sigmoid);
};

RBM.prototype.probToBinaryMat = function(m) {
	return M.mapM(m, (x)=>(Math.random() < m[i][j])?1:0);
};

RBM.prototype.propdown = function(h) {
	return h.dot(this.W.tr()).addMV(this.vbias).mapM(NN.sigmoid);
};

RBM.prototype.sampleHgivenV = function(v0_sample) {
	var h1_mean = this.propup(v0_sample);
	var h1_sample = this.probToBinaryMat(h1_mean);
	return [h1_mean,h1_sample];
};

RBM.prototype.sampleVgivenH = function(h0_sample) {
	var v1_mean = this.propdown(h0_sample);
	var v1_sample = this.probToBinaryMat(v1_mean);
	return [v1_mean,v1_sample];
};

RBM.prototype.gibbsHVH = function(h0_sample) {
	var v1 = this.sampleVgivenH(h0_sample);
	var h1 = this.sampleHgivenV(v1[1]);
	return [v1[0],v1[1],h1[0],h1[1]];
};

RBM.prototype.reconstruct = function(v) {
	var h = v.dot(this.W).addMV(this.hbias).mapM(NN.sigmoid);
	return h.dot(this.W.tr()).addMV(this.vbias).mapM(NN.sigmoid);
};

RBM.prototype.getReconstructionCrossEntropy = function() {
	var reconstructedV = this.reconstruct(this.input);
	var a = M.mapMM(this.input, reconstructedV, function(x,y){
		return x*Math.log(y);
	});

	var b = M.mapMM(this.input,reconstructedV,function(x,y){
		return (1-x)*Math.log(1-y);
	});
  var crossEntropy = -a.add(b).rowSum().mean();
	return crossEntropy

};

RBM.prototype.set = function(property,value) {
	this.settings[property] = value;
}

module.exports = RBM;