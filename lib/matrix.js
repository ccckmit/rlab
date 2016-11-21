var N = require("numeric");
var M = require("./algebra");
// var M = N;
// Advance mathematics
M.ode=N.dopri; // dopri(x0,x1,y0,f,tol,maxit,event) #Ordinary Diff Eq
M.minimize=N.uncmin; // uncmin(f,x0,tol,gradient,maxit,callback,options) # Unconstrained optimization
M.sparse=N.ccsSparse; // Matrix => Sparse
M.sparse2full=N.ccsFull; // Sparse => Matrix
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

M.mstr = M.strM = N.prettyPrint;
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
  return M.repeat([n], 0);
}

M.newM = function(rows, cols, value) {
  return M.repeat([rows, cols], 0);
}

M.randomV = function(n, a, b) {
  return M.random([n]).mul(b-a).add(a); 
}

M.randomM = function(rows, cols, a, b) {
  return N.random([rows, cols]).mul(b-a).add(a); 
}

M.rowSum = function(m) {
	var s=M.repeat(M.dim(m[0]),0);
  for (var i=0; i<m[0].length; i++) {
    s[i] = m[i].sum();
  }
  return s;
}

M.colSum = function(m) {
	var s=M.repeat(M.dim(m[0]),0);
  for (var i=0; i<m.length; i++) {
    s = s.add(m[i]);
  }
  return s;
}

/*
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
*/
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
  return M.fillVM(v,rows,cols);
}

M.eig=function(m) {
  var E = N.eig(m);
  return {lambda:E.lambda.x, E:E.E.x};
}

M.tr = M.transpose = function(m){
  var r = [];
  for(var j = 0; j < m.cols(); j++){
    r[j] = [];
    for(var i = 0; i < m.rows(); i++){
        r[j][i] = m[i][j];
    }
  }
  return r;
}

M.dot = function(a, b){
  var r = [];
  for(var i = 0; i < a.rows(); i++){
    r[i] = [];
    for(var j = 0; j < b.cols(); j++){
      r[i][j] = 0;
      for(var k = 0; k < a.cols(); k++){
        r[i][j] = r[i][j].add(a[i][k].mul(b[k][j]));
      }
    }
  }
  return r;
} 

M.diag = function(v) {
  var r = M.newM(v.length, v.length);
  for(var i = 0; i < v.length; i++){
    r[i][i] = v[i];
  }
  return r;
}

M.I = M.identity = function(n) {
  return M.diag(M.repeat([n], ()=>1));
}

M.inv = function(x) {
  var s = M.dim(x), abs = Math.abs, m = s[0], n = s[1];
  var A = M.clone(x), Ai, Aj;
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
  var s = M.dim(x);
  if(s.length !== 2 || s[0] !== s[1]) { throw new Error('numeric: det() only works on square matrices'); }
  var n = s[0], ret = 1,i,j,k,A = M.clone(x),Aj,Ai,alpha,temp,k1,k2,k3;
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

module.exports = M;

/*
M.mstr = function(m) {
	var s = ""; 
	for (var i=0; i<M.rows(m); i++) {
		var line = "";
		for (var j=0; j<M.cols(m); j++) {
			line += m[i][j].str()+",";
		}
		s+=line+"\n";
	}
	return s;
}
*/
/*
// 6. Coordinate matrices
numeric.cLU = function LU(A) {
    var I = A[0], J = A[1], V = A[2];
    var p = I.length, m=0, i,j,k,a,b,c;
    for(i=0;i<p;i++) if(I[i]>m) m=I[i];
    m++;
    var L = Array(m), U = Array(m), left = numeric.rep([m],Infinity), right = numeric.rep([m],-Infinity);
    var Ui, Uj,alpha;
    for(k=0;k<p;k++) {
        i = I[k];
        j = J[k];
        if(j<left[i]) left[i] = j;
        if(j>right[i]) right[i] = j;
    }
    for(i=0;i<m-1;i++) { if(right[i] > right[i+1]) right[i+1] = right[i]; }
    for(i=m-1;i>=1;i--) { if(left[i]<left[i-1]) left[i-1] = left[i]; }
    var countL = 0, countU = 0;
    for(i=0;i<m;i++) {
        U[i] = numeric.rep([right[i]-left[i]+1],0);
        L[i] = numeric.rep([i-left[i]],0);
        countL += i-left[i]+1;
        countU += right[i]-i+1;
    }
    for(k=0;k<p;k++) { i = I[k]; U[i][J[k]-left[i]] = V[k]; }
    for(i=0;i<m-1;i++) {
        a = i-left[i];
        Ui = U[i];
        for(j=i+1;left[j]<=i && j<m;j++) {
            b = i-left[j];
            c = right[i]-i;
            Uj = U[j];
            alpha = Uj[b]/Ui[a];
            if(alpha) {
                for(k=1;k<=c;k++) { Uj[k+b] -= alpha*Ui[k+a]; }
                L[j][i-left[j]] = alpha;
            }
        }
    }
    var Ui = [], Uj = [], Uv = [], Li = [], Lj = [], Lv = [];
    var p,q,foo;
    p=0; q=0;
    for(i=0;i<m;i++) {
        a = left[i];
        b = right[i];
        foo = U[i];
        for(j=i;j<=b;j++) {
            if(foo[j-a]) {
                Ui[p] = i;
                Uj[p] = j;
                Uv[p] = foo[j-a];
                p++;
            }
        }
        foo = L[i];
        for(j=a;j<i;j++) {
            if(foo[j-a]) {
                Li[q] = i;
                Lj[q] = j;
                Lv[q] = foo[j-a];
                q++;
            }
        }
        Li[q] = i;
        Lj[q] = i;
        Lv[q] = 1;
        q++;
    }
    return {U:[Ui,Uj,Uv], L:[Li,Lj,Lv]};
};
*/