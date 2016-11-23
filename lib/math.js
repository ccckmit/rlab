// var M = require("./calculus");
var R = require("./geometry");

// 線性最小平方迴歸 : 線性代數, Larson, 翁慶昌 5/e , 131 頁
// https://en.wikipedia.org/wiki/Least_squares
R.minSquare = function(x,y) {
	var Xt = [R.newV(x[0].length, 1)].concat(x), X = Xt.tr(); // X = [1,x]t
	var Yt = [ y ], Y = Yt.tr();               // Y = [y]t
	var A = Xt.dot(X).inv().dot(Xt).dot(Y);    // A = (XtX)-1 Xt Y
	return A;
}

// 多項式最小平方迴歸 : https://en.wikipedia.org/wiki/Polynomial_least_squares

module.exports=R;
