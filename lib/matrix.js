var N, R, M, V;
N = require("numeric");
R = module.exports = require("./vector");
M = R.M = R.Matrix = {};
V = R.V;
// var M = N;
// Advance mathematics
R.ode = N.dopri; // dopri(x0,x1,y0,f,tol,maxit,event) #Ordinary Diff Eq
R.minimize=N.uncmin; // uncmin(f,x0,tol,gradient,maxit,callback,options) # Unconstrained optimization
R.solveLP = N.solveLP; // dopri(x0,x1,y0,f,tol,maxit,event) #Ordinary Diff Eq
M.sparse=N.ccsSparse; // Matrix => Sparse
M.sparse2full=N.ccsFull; // Sparse => Matrix
R.spline = N.spline;
R.linspace = N.linspace;
// M.complex=N.t;
// matrix
M.svd = N.svd;
// M.det = N.det;
// M.inv = N.inv;
M.lu = N.cLU;
M.luSolve = N.cLUsolve;
// M.dot = N.dot;
// M.rep = N.rep;
// M.tr = N.transpose;
// M.diag = N.diag;
// M.sumM = N.sum;

M.str = M.mstr = M.strM = N.prettyPrint;
M.rows=function(m) { return m.length; }
M.cols=function(m) { return m[0].length; }
M.row =function(m,i) { return m[i]; }
M.col =function(m,j) { 
  var rows = m.length, c = new Array(rows);
  for (var i=0;i<rows;i++) {
    c[i]=m[i][j];
  }
  return c;
}

R.newM = M.newM = function(rows, cols, value=0) {
  return V.repeat([rows, cols], value);
}

R.randomM = M.random = function(rows, cols, a, b) {
  return V.random([rows, cols], a, b);
}

M.rowSum = function(m) {
	var rows=m.length, s=new Array(rows);
  for (var i=0; i<rows; i++) {
    s[i]=V.sum(m[i]);
  }
  return s;
}

M.colSum = function(m) {
	var rows = m.length;
	if (rows === 0) return [];
	var s=m[0];
  for (var i=1; i<rows; i++) {
    s = V.vadd(s, m[i]);
  }
  return s;
}

M.rowMean = function(m) { 
  return M.rowSum(m).div(m.cols());
}

M.colMean = function(m) { 
  return M.colSum(m).div(m.rows());
}

R.addMV = M.addMV = function(m,v) {
  var rows = m.length, r = new Array(rows);
  for(var i=0;i<rows;i++) {
//    r.push(m[i].add(v)); // 使用多型速度會變慢
    r[i]=V.vadd(m[i], v);  // 這行比較快
  }
  return r;
}

R.madd = M.add = function(m1, m2) {
	var rows = m1.length, r=new Array(rows); // r = [];
	for (var i = 0; i < rows; i++) {
		r[i] = R.vadd(m1[i], m2[i]);
	}
	return r;
}

R.msub = M.sub = function(m1, m2) {
	var rows = m1.length, r=new Array(rows); // r = [];
	for (var i = 0; i < rows; i++) {
		r[i] = R.vsub(m1[i], m2[i]);
	}
	return r;
}

R.mmul = M.mul = function(m1, m2) {
	var rows = m1.length, r=new Array(rows); // r = [];
	for (var i = 0; i < rows; i++) {
		r[i] = R.vmul(m1[i], m2[i]);
	}
	return r;
}

R.mdiv = M.div = function(m1, m2) {
	var rows = m1.length, r=new Array(rows); // r = [];
	for (var i = 0; i < rows; i++) {
		r[i] = R.vdiv(m1[i], m2[i]);
	}
	return r;
}

M.fillMV = function(v,rows,cols) {
  var m = new Array(rows); // []; 
  for (var r=0; r<rows; r++) {
		var mr = m[r] = new Array(cols);;
    for (var c=0; c<cols; c++) {
      mr[c] = v[r*cols+c];
    }
  }
  return m;
}

M.eig=function(m) {
  var E = N.eig(m);
  return {lambda:E.lambda.x, E:E.E.x};
}

M.tr = M.mtranspose = function(m){
  var r = [], rows = m.length, cols = m[0].length;
  for(var j = 0; j < cols; j++){
    var rj = r[j] = [];
    for(var i = 0; i < rows; i++){
        rj[i] = m[i][j];
    }
  }
  return r;
}

R.mdot = M.dot = function(a, b, isComplex=false){
	var arows = a.length, acols=a[0].length, bcols=b[0].length;
  var r = [], bt = M.tr(b);
	for (var i = 0; i < arows; i++) {
		var ri = r[i] = [];
		for (var j = 0; j < bcols; j++) {
			ri.push(V.vdot(a[i],bt[j], isComplex));
		}
	}
  return r;
}

R.isMatrix = M.isMatrix = function(m) {
	return (m instanceof Array && m[0] instanceof Array);
}

R.dot = function(a,b,isComplex=false) {
	return (R.isMatrix(a))?R.mdot(a,b,isComplex):R.vdot(a,b,isComplex);
}

M.diag = function(v) {
	var rows = v.length;
  var r = M.newM(rows, rows);
  for(var i = 0; i < rows; i++){
    r[i][i] = v[i];
  }
  return r;
}

M.id = M.identity = function(n) {
  return M.diag(V.repeat([n], ()=>1));
}

M.inv = function(x) {
  var s = V.dim(x), abs = Math.abs, m = s[0], n = s[1];
  var A = R.clone(x), Ai, Aj;
  var I = M.identity(m), Ii, Ij;
  var i,j,k,x;
  for(j=0;j<n;++j) {
    var i0 = -1;
    var v0 = -1;
    for(i=j;i!==m;++i) { k = abs(A[i][j]); if(k>v0) { i0 = i; v0 = k; } }
    Aj = A[i0]; A[i0] = A[j]; A[j] = Aj;
    Ij = I[i0]; I[i0] = I[j]; I[j] = Ij;
    x = Aj[j];
    for(k=j;k!==n;++k)    Aj[k] /= x; 
    for(k=n-1;k!==-1;--k) Ij[k] /= x;
    for(i=m-1;i!==-1;--i) {
      if(i!==j) {
        Ai = A[i];
        Ii = I[i];
        x = Ai[j];
        for(k=j+1;k!==n;++k)  Ai[k] -= Aj[k]*x;
        for(k=n-1;k>0;--k) { Ii[k] -= Ij[k]*x; --k; Ii[k] -= Ij[k]*x; }
        if(k===0) Ii[0] -= Ij[0]*x;
      }
    }
  }
  return I;
}

M.det = function(x) {
  var s = V.dim(x);
  if(s.length !== 2 || s[0] !== s[1]) { throw new Error('numeric: det() only works on square matrices'); }
  var n = s[0], ret = 1,i,j,k,A = R.clone(x),Aj,Ai,alpha,temp,k1,k2,k3;
  for(j=0;j<n-1;j++) {
    k=j;
    for(i=j+1;i<n;i++) { if(Math.abs(A[i][j]) > Math.abs(A[k][j])) { k = i; } }
    if(k !== j) {
      temp = A[k]; A[k] = A[j]; A[j] = temp;
      ret *= -1;
    }
    Aj = A[j];
    for(i=j+1;i<n;i++) {
      Ai = A[i];
      alpha = Ai[j]/Aj[j];
      for(k=j+1;k<n-1;k+=2) {
          k1 = k+1;
          Ai[k] -= Aj[k]*alpha;
          Ai[k1] -= Aj[k1]*alpha;
      }
      if(k!==n) { Ai[k] -= Aj[k]*alpha; }
  }
    if(Aj[j] === 0) { return 0; }
    ret *= Aj[j];
  }
  return ret*A[j][j];
}

R.mixThisMap(Array.prototype, M, {
lu:"lu",
luSolve:"luSolve",
svd:"svd",
// "cdelsq",
// "clone",
rows:"rows",
cols:"cols",
row:"row",
col:"col",
tr:"tr",
inv:"inv",
// "all",
// "any",
// "same",
// "isFinite",
// "isNaN",
// "mapreduce",
// "complex",
det:"det",
// "norm2",
// "norm2Squared",
// "norm2inf",
madd:"add",
msub:"sub",
mmul:"mul",
mdiv:"div",
mdot:"dot",
// "dim",
eig:"eig",
// "sum",
rowSum:"rowSum",
colSum:"colSum",
rowMean:"rowMean",
colMean:"colMean",
// mapM:"mmap1",
// mapMM:"mmap2",
flat:"flat",
fillVM:"fillVM",
addMV:"addMV",
// fillMM:"fillMM",
getBlock:"getBlock",
setBlock:"setBlock",
getDiag:"getDiag",
diag:"diag",
// "parseFloat",
// "parseDate",
// "parseCSV",
// "toCSV",
strM:"strM", 
mstr:"mstr", 
// "sumM",
isMatrix:"isMatrix",
});

R.mixThisMap(Array.prototype, R, {
dot:"dot",
});

// console.log('msub=', [[1,2],[3,4]].msub);