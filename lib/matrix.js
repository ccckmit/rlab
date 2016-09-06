var N = require("numeric");
var M = N;
// Advance mathematics
M.ODE=N.dopri; // dopri(x0,x1,y0,f,tol,maxit,event) #Ordinary Diff Eq
M.minimize=N.uncmin; // uncmin(f,x0,tol,gradient,maxit,callback,options) # Unconstrained optimization
M.sparse=N.ccsSparse; // Matrix => Sparse
M.sparse2full=N.ccsFull; // Sparse => Matrix
M.complex=N.t;
// matrix
M.tr = N.transpose;
M.str = M.strM = N.prettyPrint;
M.sumM = N.sum;
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
		s[i] = m[i].sum();
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

M.fillVM=function(v,rows,cols) {
	var m = M.newM(rows,cols);
	for (var r=0; r<rows; r++) {
		for (var c=0; c<cols; c++) {
			m[r][c] = v[r*cols+c];
		}
	}
	return m;
}

M.fillMM=function(m,rows,cols) {
	var v = M.flatM(m);
	return M.fillVM(m,rows,cols);
}

M.eigR=function(m) {
	var E = M.eig(m);
	return {lambda:E.lambda.x, E:E.E.x};
}

module.exports = M;