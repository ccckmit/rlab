// module : Field & Group Theory
var Op = require("./op");
var extend = Object.assign;
var F=Op.Field, Field=F.Field, Group=F.Group, FieldObj=F.FieldObj;
var A = { Op:Op }
/*
A.gcd = function(a, b) {
  if (!b) return a;
  return I.gcd(b, a % b);
}

A.lcm = function(a, b) {
  return (a * b) / gcd(a, b);   
}

A.isPrime=function(n) {
	for (var i=2; i<n/2; i++)
		if (n%i===0) return false;
	return n%1===0;
}
*/
// =========== Polynomial Field ==============
A.PolynomialRing=extend({}, Field);

class Polynomial extends FieldObj {
  constructor(c) {
    super(A.PolynomialRing);
    this.c = c; // sum(ci x^i)
  }
	eval(x) {
		var c = this.c, i=c.length-1, sum=c[i];
		for (i--;i>=0; i--) {
			sum = c[i].add(x.mul(sum));
		}
		return sum;
	}
	size() { return this.c.length }
	toString() {
		var s = [], c=this.c;
		for (var i=c.length-1; i>=0; i--) {
			s.push(c[i]+'x^'+i);
		}
		return s.join('+');
	}
	root2(b,c) {
		var d = b.mul(b).sub(c.mul(4)).sqrt();
		return -b
	}
}

A.Polynomial = Polynomial;
/*
function root2(a,b,c) { // ax^2+bx+c=0
	var t = b*b-4*a*c
	if (t < 0) throw Error('沒有實根');
	var t2 = Math.sqrt(t);
	return [(-b+t2)/(2*a), (-b-t2)/(2*a)];
}

console.log("root2(1x^2+4x+0)=", root2(1,4,0));

function root3(a,b,c) { // x^3+ax^2+bx+c=0
  var q=((2*a*a*a/27)-(a*b/3)+c/2);
	var p=(b-a*a/3)/3;
	var D=p*p*p+q*q;
	var u_p=Math.sqrt(3, -q+Math.sqrt(D));
	var u_m=Math.sqrt(3, -q-Math.sqrt(D));
	var rho_p = 1/2*(-1+3i);
	var rho_m = 1/2*(-1-3i);
  var y1=u_p+u_m;
  var y2=rho_p+u_p+ rho_m*u_m;
  var y3=rho_p-u_p+ rho_p*u_m;
}
*/
A.PolynomialAddGroup={
  e:new Polynomial([0]), // 到底應該設幾個？
  op:function(x,y) {
		var size = Math.max(x.size(), y.size());
		var a=Op.vextend(x.c,size), b=Op.vextend(y.c,size);
	  var c= Op.add(a,b);
		return new Polynomial(c);
	},
  inv:function(x) {
		var c = Op.neg(x.c);
		return new Polynomial(c);
	},
}
  
extend(A.PolynomialAddGroup, Group);

A.PolynomialMulSemiGroup={
  e:new Polynomial([1]),
  op:function(x,y) { 
	  var c=[];
	  for (var xi=0; xi<x.c.length; xi++) {
			for (var yi=0; yi<y.c.length; yi++) {
				var cxy = (typeof c[xi+yi]==='undefined')?0:c[xi+yi];
				c[xi+yi] = cxy+x.c[xi]*y.c[yi];
			}
		}
		return new Polynomial(c);
	},
  inv:function(x) { throw Error('PolynomialMulSemiGroup.inv not defined') },
}

extend(A.PolynomialMulSemiGroup, Group);

A.PolynomialRing.init(A.PolynomialAddGroup, A.PolynomialMulSemiGroup);

module.exports = A;

/*
// =========== Vector Field ==============
var VectorField=extend({}, Field);

class Vector extends FieldObj {
  constructor(a) {
    super(VectorField);
    this.a = a;
  }
	size() { return a.length }
	get(i) { return a[i] }
}

var VectorAddGroup={
  e:{ get:function(i) { return 0} },
  op:function(x,y) { 
	  var c = [];
	  for (var i=0; i<x.size(); i++) {
			c[i] = x.get(i).add(y.get(i));
		}
		return new Vector(c);
	},
  inv:function(x) { 
	  var a = [];
		for (var i=0; i<x.size(); i++) {
			a[i] = x.get(i).neg();
		}
	},
}
  
extend(VectorAddGroup, Group);

var VectorMulGroup={
  e:{ get:function(i) { return 1} },
  op:function(x,y) { 
	  var c = [];
	  for (var i=0; i<x.size(); i++) {
			c[i] = x.get(i).mul(y.get(i));
		}
		return new Vector(c);
	},
  inv:function(x) { 
	  var a = [];
		for (var i=0; i<x.size(); i++) {
			a[i] = x.get(i).inv();
		}
		return new Vector(a);
	},
}

extend(VectorMulGroup, Group);

VectorField.init(VectorAddGroup, VectorMulGroup);

*/
/*
// =========== Function Field ==============
var FunctionRing=extend({}, Field);

class FunctionObj extends FieldObj {
  constructor(f) {
    super(FunctionRing);
    this.f = f;
  }
	
	eval(v) { return this.f(v) }
}

var FunctionAddGroup={
  e:new FunctionObj(function(x) {return 0}),
  op:function(x,y) { 
	  return new FunctionObj(function f(v) {
			return x.eval(v).add(y.eval(v));
		})
	},
  inv:function(x) { 
	  return new FunctionObj(function f(v) {
			return -1*x.eval(v);
		})
	},
}
  
extend(FunctionAddGroup, Group);

var FunctionMulSemiGroup={
  e:new FunctionObj(function f(v) { return v }),
  op:function(x,y) { 
	  return new FunctionObj(function f(v) {
			return x.eval(y.eval(v));
		})
	},
  inv:function(x) { throw Error('Function SemiGroup has no inverse!') },
}

extend(FunctionMulSemiGroup, Group);

FunctionRing.init(FunctionAddGroup, FunctionMulSemiGroup);

*/