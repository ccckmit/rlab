module.exports = O = require("./matrix");

// =========== Calculus =================
O.dx = 0.01;

// 微分 differential calculus
O.fdiff = O.fdifferential = function(f, x, dx) {
	dx = dx || O.dx;
  var dy = f(x.add(dx)).sub(f(x.sub(dx)));
  return dy.div(dx.mul(2));
}

// 積分 integral calculus
O.fint = O.fintegral = function(f, a, b, dx) {
	dx = dx || O.dx;
  var area = 0.0;
  for (var x=a; x<b; x=x+dx) {
    area = area + f(x).mul(dx);
  }
  return area;
}

// 偏微分 partial differential calculus
// f=[f1,f2,....] , x=[x1,x2,...] , dx=[dx1,dx2,....]
O.pdiff = O.pdifferential = function(f, x, i) {
	f = O.fa(f);
	var dx = O.vfill(x.length, 0);
	dx[i] = O.dx;
//	var df = f.apply(null, x.add(dx)).sub(f.apply(null, x.sub(dx)));
	var df = f(x.add(dx)).sub(f(x.sub(dx)));
	return df.div(dx.norm().mul(2));
}

// multidimensional integral calculus
// f=[f1,f2,....] , a=[a1,a2,...] , b=[b1,b2,....]
O.pint = O.pintegral = function(f, a, b) {
	
}

// 梯度 gradient : grad(f,x)=[pdiff(f,x,0), .., pdiff(f,x,n)]
O.fgrad = O.fgradient = function(f, x) {
	var gf = [];
	for (var i=0; i<x.length; i++) {
		gf[i] = O.pdiff(f, x, i);
	}
	return gf;
}

// 散度 divergence : div(F,x) = sum(pdiff(F[i],x,i))
O.fdiv = O.fdivergence = function(F, x) {
	var Fx = F(x);
	var f=[], d=[];
	for (var i=0; i<x.length; i++) {
		f[i] = (xt)=>F(xt)[i];
		d[i] = O.pdiff(f[i],x,i);
	}
	return d.sum();
}

// 旋度 curl : curl(F) = div(F)xF 
O.fcurl = O.fcurlance = function(F, x) {
	
}

// 線積分： int F●dr = int F(r(t))●r'(t) dt
O.vint = O.vintegral = function(F, r, a, b, dt) {
	dt = dt||O.dx;
	var sum = 0;
	for (var t=a; t<b; t+=dt) {
		sum += F(r(t)).dot(r.diff(t));
	}
	return sum;
}

// 定理：int F●dr  = int(sum(Qxi dxi))

// 向量保守場： F=grad(f) ==> int F●dr = f(B)-f(A)

// 定理： 向量保守場 F, pdiff(Qxi,xj) == pdiff(Qxj,xi)  for any i, j
// ex: F=[x^2y, x^3] 中， grad(F)=[x^2, 3x^2] 兩者不相等，所以不是保守場

// 格林定理：保守場中 《線積分=微分後的區域積分》
// int P dx + Q dy = int int pdiff(Q,x) - pdiff(P, y) dx dy

// 散度定理：通量 = 散度的區域積分
// 2D : int F●n ds = int int div F dx dy
// 3D : int int F●n dS = int int int div F dV

// 史托克定理：  F 在曲面 S 上的旋度總和 = F 沿著邊界曲線 C 的線積分
// int int_D curl(F)●n dS = int_C F●dr

