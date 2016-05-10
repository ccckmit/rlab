var B = require("./base");
var N = require("numeric");
var M = {};

M.tr = N.transpose;
M.str = N.prettyPrint;
M.rows=function(m) { return m.length; }
M.cols=function(m) { return m[0].length; }
M.row =function(m,i) { return m[i]; }
M.col =function(m,j) { 
  var cols = m.cols();
  var c = M.newV(cols);
	for (var i=0;i<cols;i++) {
		c[i] = m[i][j];
	}
	return c;
}

M.newV = function(n, value) {
	return M.rep([n], value||0);
}

M.newM = function(rows, cols, value) {
	return M.rep([rows, cols], value||0);
}

M.randomV = function(n, a, b) { 
  return M.random([n]).mul(b-a).add(a); 
}

M.randomM = function(rows, cols, a, b) { 
  return M.random([rows, cols]).mul(b-a).add(a); 
}

M.rowSum = function(m) {
	var rows = M.rows(m);
	var s=M.newV(rows, 0);
	for (var i=0; i<rows; i++) {
		s[i] = m[i].sumM();
	}
	return s;
}

M.colSum = function(m) {
	var mt = M.tr(m);
	return M.rowSum(mt);
}

M.rowMean = function(m) { 
  return M.rowSum(m).div(m.cols());
}

M.colMean = function(m) { 
  return M.colSum(m).div(m.rows());
}

M.addMV = function(m,v) {
  var result = [];
  for(var i=0;i<m.length;i++) {
    result.push(m[i].add(v));
	}
  return result;
}

M.mapM = function(m, f) {
	var fm = M.clone(m);
	var rows = M.rows(m), cols=M.cols(m);
  for(i=0;i<rows;i++) {
    for(j=0;j<cols;j++)
      fm[i][j]=f(m[i][j]);
  }
  return fm;
}

M.mapMM = function(m1,m2,f) {
	var fm = M.clone(m1);
	var rows = m1.rows(), cols=m1.cols();
  for(i=0;i<rows;i++) {
    for(j=0;j<cols;j++)
      fm[i][j]=f(m1[i][j],m2[i][j]);
  }
  return fm;
}

M.flatM=function(m) {
	var a=[];
	var ai = 0;
	for (var i=0; i<m.length;i++)
		for (var j=0; j<m[i].length; j++)
			a[ai++] = m[i][j];
	return a;
}

M.eigR=function(m) {
	var E = M.eig(m);
	return {lambda:E.lambda.x, E:E.E.x};
}

B.mix(M, N);

module.exports = M;