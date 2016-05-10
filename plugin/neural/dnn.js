var RBM = require("./rbm");

module.exports = dnn = {};

dnn.RBM = function(setting) {
	setting.nVisible = setting.input[0].length;
  return new RBM(setting); 
}
