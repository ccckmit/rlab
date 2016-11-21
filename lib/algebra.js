var O = require("./field");
var extend = O.extend;
// ============= Array Operation =========================
var acc = {
	sum:{init:0,op:'+'}, 
	product:{init:1,op:'*'},
	min:{init:Number.MAX_VALUE,op:'min'},
	max:{init:-Number.MAX_VALUE,op:'max'},
}

O.recursiveDim = function(a,d) { // dimension by recursion
	if (O.isArray(a)) {
		d.push(a.length)
		O.recursiveDim(a[0],d);
	}
}

O.dim = function(a) {
	var d=[];
	O.recursiveDim(a,d);
	return d;
}

O.vfill=function(size, value) {
	var v=[];
	for (var i=0; i<size; i++)
		v[i] = value;
	return v;
}

O.recursiveRepeat = function(dim, i, v) {
	if (i===dim.length) {
		return O.isFunction(v)?v():v;
	} else {
		var a = [];
		for (var k=0; k<dim[i]; k++) {
			a.push(O.recursiveRepeat(dim, i+1, v));
		}
		return a;
	}
}

O.repeat = function(dim, v) {
	v = v||0;
	return O.recursiveRepeat(dim, 0, v);
}

O.vrandom = function(dim) {
	return O.repeat(dim, ()=>O.random(0,1));
}

O.vextend=function(a, size) {
	var v = a.slice();
	for (var i=a.length; i<size; i++) {
		v.push(0);
	}
	return v;
}

O.vdot=function(x,y) {
	var sum = 0;
	for (var i=0; i<x.length; i++) {
		sum = sum.add(x[i].mul(y[i]));
	}
	return sum;
}

O.vop=function(op,x,y) {
	var c=[],result=(O.isUndefined(acc[op]))?undefined:acc[op].init;
  if (O.isArray(x)) {
		for (var i=0; i<x.length; i++) {
			if (!O.isUndefined(result)) {
				var xi = x[i];
				if (O.isArray(x[i])) // x[i] is subArray
					xi = O.op(op, x[i]);
				else {
					if (!O.isUndefined(y))
						xi = y(x[i]); // y should be a function
				}
				result = O.op(acc[op].op, xi, result);
//				if (op==='min') console.log("xi=%d result=%d", xi, result);
			}
			else {
				var yi = O.isArray(y)?y[i]:y;
				c[i] = O.op(op, x[i], yi);
			}
		}
	} else {
		throw Error('vop fail:op,x,y=', op, x, y);
	}
	if (x.length === c.length)
		return c;
	else
		return result;
}

// ============= Operation =========================
O.op = function(op,x,y) {
  if (O.isField(x)) {
		if (!O.isUndefined(acc[op])) {
			switch (op) {
				case 'max': return Math.max(x,y);
				case 'min': return Math.min(x,y);
				default:return x;
			}
		}	else if (x instanceof O.Complex || 
		    (O.isNumber(x) && y instanceof O.Complex)) {
			x = O.Complex.toComplex(x); 
			y = O.isNumber(y)?O.Complex.toComplex(y):y;
			switch (op) {
				case 'eval': var exp = y; return exp(x);
				case 'neg':return x.neg();
				case 'inv':return x.inv();
				case 'bnot':return x.bnot();				
				case '+':return x.add(y);
				case '-':return x.sub(y);
				case '*':return x.mul(y);
				case '/':return x.div(y);
				case 'power':return x.power(y);
				case 'sqrt':return x.sqrt();
				case 'eq':return x.eq(y);
				case 'neq':return x.neq(y);
				case 'geq':return x.geq(y);
				case 'leq':return x.leq(y);
				case 'gt':return x.gt(y);
				case 'lt':return x.lt(y);
			}
		} else if (O.isBool(x) || O.isNumber(x)) {
			switch (op) {
				case 'eval': var exp = y; return exp(x);
				case 'not':return !x;
				case 'neg':return -x;
				case 'inv':return 1/x;
				case 'bnot':return ~x;
				case 'and':return x&&y;
				case 'or':return x||y;
				case 'xor':return x!==y;
				case '+':return x+y;
				case '-':return x-y;
				case '*':return O.isArray(y)?O.op(op,y,x):x*y;
				case '/':return O.isArray(y)?O.op('*',1/x,y):x/y;
				case '%':return x%y;
				case 'eq':return x===y;
				case 'neq':return x!==y;
				case 'geq':return x>=y;
				case 'leq':return x<=y;
				case 'gt':return x>y;
				case 'lt':return x<y;
				case '&':return x&y;
				case '|':return x|y;
				case 'bxor':return x^y;
				case '<<':return x<<y;
				case '>>':return x>>y;
				case 'and':return x&&y;
				case 'or':return x||y;
				case 'xor':return x!==y;
				case 'sqrt':
				  if (x>=0) return Math.sqrt(x);
					else {
						var c = new O.Complex(x,0);
						return c.sqrt();
					}
				case 'power':
				  if (x>=0 || O.isInteger(y))
						return Math.pow(x,y);
					else
						return new O.Complex(x,0).power(y);
				case 'log':return Math.log(x);
				case 'exp':return Math.exp(x);
				case 'abs':return Math.abs(x);;
				case 'sin':return Math.sin(x);
				case 'cos':return Math.cos(x);
				case 'tan':return Math.tan(x);
				case 'asin':return Math.asin(x);
				case 'acos':return Math.acos(x);
				case 'atan':return Math.atan(x);
				case 'sinh':return Math.sinh(x);
				case 'cosh':return Math.cosh(x);
				case 'tanh':return Math.tanh(x);
				case 'ceil':return Math.ceil(x);
				case 'floor':return Math.floor(x);
				case 'round':return Math.round(x);
				case 'log1p':return Math.log1p(x);
				case 'log10':return Math.log10(x);
				case 'log2':return Math.log2(x);
				case 'random':return Math.random();
				case 'sign':return Math.sign(x);
				case 'abs':return Math.abs(x);
				case 'cbrt':return Math.cbrt(x); // cubic root
			}
		}
	} else if (O.isFunction(x)) {
		if (O.isFunction(y)) {
			switch (op) {
				case 'neg':return O.fneg(x);
				case 'inv':return O.finv(x);
				case '+':return O.fadd(x,y);
				case '-':return O.fsub(x,y);
				case '*':return O.fmul(x,y);
				case '/':return O.fdiv(x,y);
				case 'compose':return O.fcompose(x,y);
			}			
		} else {
			switch (op) {
				case 'eval':return x(y);
			}
		}
	}	else if (O.isArray(x)) {
		return O.vop(op,x,y);
	}
	throw Error('op fail:op,x,y=', op, x, y);
}

// ============= Array Function =======================
O.sum=function(x) { return O.vop('sum', x) }
O.product=function(x) { return O.vop('product', x) }
O.max=function(x) { return O.vop('max', x) }
O.min=function(x) { return O.vop('min', x) }
O.eval=O.map=function(x,y) { return O.op('eval', x,y) }
O.norm=function(x) {
	var norm2 = O.vop('sum', x, (x)=>x*x);
	return Math.sqrt(norm2);
}

// +-*/%^
O.add=function(x,y) {	return O.op('+',x,y) }
O.sub=function(x,y) {	return O.op('-',x,y) }
O.mul=function(x,y) {	return O.op('*',x,y) }
O.div=function(x,y) {	return O.op('/',x,y) }
O.mod=function(x,y) {	return O.op('%',x,y) }
O.power=function(x,y) { return O.op('power', x, y) }
O.neg=function(x) { return O.op('neg', x) }
O.inv=function(x) { return O.op('inv', x) }
// logical
O.not=function(x) { return O.op('not', x) }
O.and=function(x,y) { return O.op('&&', x, y) }
O.or=function(x,y) { return O.op('||', x, y) }
O.xor=function(x,y) { return O.op('xor', x, y) }
O.bnot=function(x) { return O.op('bnot', x) }
O.band=function(x,y) { return O.op('&', x, y) }
O.bor=function(x,y) { return O.op('|', x, y) }
O.bxor=function(x,y) { return O.op('bxor', x, y) }
O.lshift=function(x,y) { return O.op('<<', x, y) }
O.rshift=function(x,y) { return O.op('>>', x, y) }
// compare
O.eq=function(x,y) { return O.op('eq', x, y) }
O.neq=function(x,y) { return O.op('neq', x, y) }
O.geq=function(x,y) { return O.op('geq', x, y) }
O.leq=function(x,y) { return O.op('leq', x, y) }
O.gt=function(x,y) { return O.op('gt', x, y) }
O.lt=function(x,y) { return O.op('lt', x, y) }
// number function
O.sqrt=function(x) { return O.op('sqrt', x) }
O.log=function(x) { return O.op('log', x) }
O.exp=function(x) { return O.op('exp', x) }
O.abs=function(x) { return O.op('abs', x) }
O.sin=function(x) { return O.op('sin', x) }
O.cos=function(x) { return O.op('cos', x) }
O.tan=function(x) { return O.op('tan', x) }
O.asin=function(x) { return O.op('asin', x) }
O.acos=function(x) { return O.op('acos', x) }
O.atan=function(x) { return O.op('atan', x) }
O.atan2=function(x) { return O.op('atan2', x) }
O.ceil=function(x) { return O.op('ceil', x) }
O.floor=function(x) { return O.op('floor', x) }
O.round=function(x) { return O.op('round', x) }
 
O.parse = function(s) {
	if (s.indexOf(';')>=0) {
		var m = split(s, ";"), matrix;
		for (var i=0; i<m.length; i++) {
			matrix[i] = O.parse(m[i]);
		}
		return matrix;
	} if (s.indexOf(',')>=0) {
		var a = split(s, ","), array;
		for (var i=0; i<a.length; i++) {
			array[i] = O.parse(a[i]);
		}
		return array;
	}
	else if (s.indexOf('/')>=0)
		return O.Ratio.parse(s);
	else if (s.indexOf('i')>=0)
		return O.Complex.parse(s);
	else {
		return parseFloat(s);
	}
}

// =========== Polynomial Ring ==============
O.PolynomialRing=extend({}, O.Field);

class Polynomial extends O.FieldObj {
  constructor(c) {
    super(O.PolynomialRing);
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
	
	root() {
		var p = this.normalize(); // 正規化
		switch (this.size()) {
		  case 2:return p.c[0].neg();
		  case 3:return p.root2(p.c[1], p.c[0]);
			case 4:return p.root3(p.c[2], p.c[1], p.c[0]);
			default:throw Error('root() fail');
		}
			
	}

	root2(b,c) { // x^2 + bx + c	=0
		var d = b.mul(b).sub(c.mul(4)).sqrt();
		return [b.neg().add(d), b.neg().sub(d)];
	}
	
	root3(a,b,c) { // x^3+ax^2+bx+c=0
		var q=((2*a*a*a/27)-(a*b/3)+c)/2;
		var p=(b-a*a/3)/3;
		var D=p*p*p+q*q;
		var Dsqrt = D.sqrt(), _q=q.neg();
		var u_p=_q.add(Dsqrt).power(1/3); // (-q+sqrt(D))^1/3
		var u_m=_q.sub(Dsqrt).power(1/3); // (-q-sqrt(D))^1/3
		console.log("q=%s p=%s D=%s u+=%s u-=%s", q, p, D, u_p, u_m);
		var rho_p = (1/2).mul(O.parse('-1+3i'));
		var rho_m = (1/2).mul(O.parse('-1-3i'));
		console.log("rho_p=%s rho_m=%s", rho_p, rho_m);
		var y1=u_p.add(u_m);
		var y2=rho_p.add(u_p).add(rho_m.mul(u_m));
		var y3=rho_p.sub(u_p).add(rho_p.mul(u_m));
		return [y1, y2, y3];
	}	
	
	normalize() {
		var a = this.c[this.size()-1];
		var nc = this.c.div(a);
		return new Polynomial(nc);
	}
}

O.Polynomial = Polynomial;

O.PolynomialAddGroup={
  e:new Polynomial([0]), // 到底應該設幾個？
  op:function(x,y) {
		var size = Math.max(x.size(), y.size());
		var a=O.vextend(x.c,size), b=O.vextend(y.c,size);
	  var c=O.add(a,b);
		return new Polynomial(c);
	},
  inv:function(x) {
		var c = O.neg(x.c);
		return new Polynomial(c);
	},
}
  
extend(O.PolynomialAddGroup, O.Group);

O.PolynomialMulSemiGroup={
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

extend(O.PolynomialMulSemiGroup, O.Group);

O.PolynomialRing.init(O.PolynomialAddGroup, O.PolynomialMulSemiGroup);

module.exports = O;
