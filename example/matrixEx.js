var M = require("../rlab").M;
var v = [1,2,3];
log("v.sin()=", v.sin());
log("v.norm2()=", v.norm2());
log("v.norm2Squared()=", v.norm2Squared());

var A = [[1,2,3],[4,5,6],[7,3,9]];
var AiA = A.inv().dot(A);
log("AiA=\n", AiA.strM());
log("AiA.tr()=\n", AiA.tr().strM());
log("A=\n", A.str());
log("A.mul(0.1)=\n", A.mul(0.1).strM());
log("A.row(1)=", A.row(1));
log("A.col(1)=", A.col(1));
log("A.sumM()=", A.sumM());
log("A.rowSum()=", A.rowSum());
log("A.colSum()=", A.colSum());
log("A.mean(row)=", A.rowMean().str());
log("A.mean(col)=", A.colMean().str());

var D = M.diag(v);
log("D=", D);

log("===========eigen================");
var Eλ = M.eigR(A);
var E = Eλ.E, λ=Eλ.lambda;
log("E*[λ]*E-1=", E.dot(λ.diag()).dot(E.inv()).strM());

log("===========LU================");
var lu = M.cLU([[0,0,1,1,1,2,2],[0,1,0,1,2,1,2],[2,-1,-1,2,-1,-1,2]]);
log('lu:', lu.strM());
var luSolve = M.cLUsolve(lu,[5,-8,13]);
log('luSolve:', luSolve.strM());

log("===========Sparse================");
var S = [
[0,0,0,0,0,0],
[0,3,0,0,0,0],
[0,0,0,6,0,0],
[0,0,9,0,0,0],
[0,0,0,0,12,0],
[0,0,0,0,0,5],
[0,0,1,1,0,0],
[0,0,0,0,0,0],
];
log("sparse(S)=", M.sparse(S)); 
// The relation between A and its sparse representation SA is:
//  A[i][SA[1][k]] = SA[2][k] with SA[0][i] ≤ k < SA[0][i+1]
// 雖然有點不同，但可參考： http://openhome.cc/Gossip/AlgorithmGossip/SparseMatrix.htm
// sparse(S)=[[0,  0, 1,     3,    5, 6,  7 ], // (行)
//                 [ 1,  3, 6, 2, 6, 4,  5 ],  // (列)
//                 [ 3,  9, 1, 6, 1, 12, 5 ] ] // (值)
