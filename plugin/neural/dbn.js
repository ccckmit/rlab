/**
 * Created by joonkukang on 2014. 1. 13..
 */
// var math = require('./utils').math;
// HiddenLayer = require('./HiddenLayer');
// RBM = require('./RBM');
// MLP = require('./MLP');
var NN=R.NN, RBM = NN.RBM, MLP = NN.MLP, HiddenLayer=NN.HiddenLayer, sigmoid=NN.sigmoid;

DBN = module.exports = function (settings) {
    this.x = settings['input'];
    this.y = settings['label'];
    this.sigmoidLayers = [];
    this.rbmLayers = [];
    this.hiddenLayerSizes = settings['hiddenLayerSizes'];
    this.nLayers = settings['hiddenLayerSizes'].length;
    this.nIns = settings['nIns'];
    this.nOuts = settings['nOuts'];

    // Constructing Deep Neural Network
    var i;
    for(i=0 ; i<this.nLayers ; i++) {
        var inputSize, layerInput;
        inputSize = (i===0)?settings['nIns']:settings['hiddenLayerSizes'][i-1];

        if(i == 0)
            layerInput = this.x;
        else
            layerInput = this.sigmoidLayers[this.sigmoidLayers.length-1].sampleHgivenV();

        var sigmoidLayer = new HiddenLayer({
            'input' : layerInput,
            'nIn' : inputSize,
            'nOut' : settings['hiddenLayerSizes'][i],
            'activation' : sigmoid
        });
        this.sigmoidLayers.push(sigmoidLayer);

        var rbmLayer = new RBM({
            'input' : layerInput,
            'nVisible' : inputSize,
            'nHidden' : settings['hiddenLayerSizes'][i]
        });
        this.rbmLayers.push(rbmLayer);
    }
    this.outputLayer = new HiddenLayer({
        'input' : this.sigmoidLayers[this.sigmoidLayers.length-1].sampleHgivenV(),
        'nIn' : settings['hiddenLayerSizes'][settings['hiddenLayerSizes'].length - 1],
        'nOut' : settings['nOuts'],
        'activation' : sigmoid
    });
};

DBN.prototype.pretrain = function (settings) {
    var lr = 0.6, k = 1, epochs = 2000;
    if(typeof settings['lr'] !== 'undefined')
        lr = settings['lr'];
    if(typeof settings['k'] !== 'undefined')
        k = settings['k'];
    if(typeof settings['epochs'] !== 'undefined')
        epochs = settings['epochs'];

    var i,j;
    for(i=0; i<this.nLayers ; i++) {
        var layerInput ,rbm;
        if (i==0)
            layerInput = this.x;
        else
            layerInput = this.sigmoidLayers[i-1].sampleHgivenV(layerInput);
        rbm = this.rbmLayers[i];
        rbm.logLevel = 0;
        rbm.train({
            'lr' : lr,
            'k' : k,
            'input' : layerInput,
            'epochs' : epochs
        });

        console.log("DBN RBM",i,"th Layer Final Cross Entropy: ",rbm.getReconstructionCrossEntropy());
        console.log("DBN RBM",i,"th Layer Pre-Training Completed.");

        // Synchronization between RBM and sigmoid Layer
        this.sigmoidLayers[i].W = rbm.W;
        this.sigmoidLayers[i].b = rbm.hbias;
    }
    console.log("DBN Pre-Training Completed.")
};

DBN.prototype.finetune = function (settings) {
    var lr = 0.2, epochs = 1000;
    if(typeof settings['lr'] !== 'undefined')
        lr = settings['lr'];
    if(typeof settings['epochs'] !== 'undefined')
        epochs = settings['epochs'];

    //Fine-Tuning Using MLP (Back Propagation)
    var i;
    var pretrainedWArray = [], pretrainedBArray = []; // HiddenLayer W,b values already pretrained by RBM.
    for(i=0; i<this.nLayers ; i++) {
        pretrainedWArray.push(this.sigmoidLayers[i].W);
        pretrainedBArray.push(this.sigmoidLayers[i].b);
    }
    // W,b of Final Output Layer are not involved in pretrainedWArray, pretrainedBArray so they will be treated as undefined at MLP Constructor.
    var mlp = new MLP({
        'input' : this.x,
        'label' : this.y,
        'n_ins' : this.nIns,
        'nOuts' : this.nOuts,
        'hiddenLayerSizes' : this.hiddenLayerSizes,
        'w_array' : pretrainedWArray,
        'b_array' : pretrainedBArray
    });
    mlp.train({
        'lr' : lr,
        'epochs' : epochs
    });
    for(i=0; i<this.nLayers ; i++) {
        this.sigmoidLayers[i].W = mlp.sigmoidLayers[i].W;
        this.sigmoidLayers[i].b = mlp.sigmoidLayers[i].b;
    }
    this.outputLayer.W = mlp.sigmoidLayers[this.nLayers].W;
    this.outputLayer.b = mlp.sigmoidLayers[this.nLayers].b;

};

DBN.prototype.getReconstructionCrossEntropy = function() {
    var reconstructedOutput = this.predict(this.x);
    var a = R.map2(this.y,reconstructedOutput,function(x,y){
        return x*Math.log(y);
    });

    var b = R.map2(this.y,reconstructedOutput,function(x,y){
        return (1-x)*Math.log(1-y);
    });

    var crossEntropy = -a.madd(b).colSum().mean();
    return crossEntropy
};

DBN.prototype.predict = function (x) {
    var layerInput = x, i;
    for(i=0; i<this.nLayers ; i++) {
        layerInput = this.sigmoidLayers[i].output(layerInput);
    }
    var output = this.outputLayer.output(layerInput);
    return output;
};

