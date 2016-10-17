var M = require("../rlab").M;
var v = [1,2,3];
print("v.sin()=", v.sin().str());
print("v.norm2()=", v.norm2().str());
print("v.norm2Squared()=", v.norm2Squared().str());

var A = [[1,2,3],[4,5,6],[7,3,9]];
var AiA = A.inv().dot(A);
print("AiA=\n", AiA.strM());
print("AiA.tr()=\n", AiA.tr().strM());
print("A=\n", A.str());
// print("A.mul(0.1)=\n", A.mul(0.1).strM()); 
// 目前只允許矩陣乘矩陣，要加上矩陣 +-*/ 數量 (R, C, ....)
print("A.row(1)=", A.row(1));
print("A.col(1)=", A.col(1));
print("A.sumM()=", A.sum());
print("A.rowSum()=", A.rowSum());
print("A.colSum()=", A.colSum());
print("A.mean(row)=", A.rowMean().str());
print("A.mean(col)=", A.colMean().str());

var D = M.diag(v);
print("D=", D);

print("===========eigen================");
var Eλ = M.eigR(A);
var E = Eλ.E, λ=Eλ.lambda;
print("E*[λ]*E-1=", E.dot(λ.diag()).dot(E.inv()).strM());

print("===========LU================");
var lu = M.cLU([[0,0,1,1,1,2,2],[0,1,0,1,2,1,2],[2,-1,-1,2,-1,-1,2]]);
print('lu:', M.str(lu));
var luSolve = M.cLUsolve(lu,[5,-8,13]);
print('luSolve:', luSolve.strM());

print("===========Sparse================");
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
print("sparse(S)=", M.sparse(S)); 
// The relation between A and its sparse representation SA is:
//  A[i][SA[1][k]] = SA[2][k] with SA[0][i] ≤ k < SA[0][i+1]
// 雖然有點不同，但可參考： http://openhome.cc/Gossip/AlgorithmGossip/SparseMatrix.htm
// sparse(S)=[[0,  0, 1,     3,    5, 6,  7 ], // (行)
//                 [ 1,  3, 6, 2, 6, 4,  5 ],  // (列)
//                 [ 3,  9, 1, 6, 1, 12, 5 ] ] // (值)
