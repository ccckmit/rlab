// module : Field & Group Theory
var I=require("./integer");
var R=require("./rule");
var S=require("./set");
var extend = Object.assign;
var eq = R.eq;
var F = {}

// ========== Group =================
// 注意： 箭頭函數會自動將 this 變數綁定到其定義時所在的物件，因此以下很多地方不能用箭頭函數。
// 參考： https://developer.mozilla.org/zh-TW/docs/Web/JavaScript/Reference/Functions/Arrow_functions
var Group=F.Group={ 
  invOp:function(x,y) { 
    return this.op(x,this.inv(y)); 
  },
  power:function(x,n) {
    var p=this.e;
    for (var i=0;i<n;i++) {
      p=this.op(p,x);
    }
    return p;
  },
  // ref:https://en.wikipedia.org/wiki/Group_(mathematics)
  // 封閉性：For all a, b in G, a • b, is also in G
  closability:function(a,b) {
    return this.has(this.op(a,b));
  },
  // 結合性：For all a, b and c in G, (a • b) • c = a • (b • c).
  associativity:function(a,b,c) {
    var op = this.op.bind(this);
    return eq(op(op(a,b),c), op(a,op(b,c)))
  },
  // 單位元素：Identity element
  identity:function(a) {
    return eq(this.op(this.e,a),a)
  },
  // 反元素：Inverse element
  inversability:function(a) {
    return eq(this.op(a,this.inv(a)),this.e);
  },
}

// ========== Field =================
var Field=F.Field={
  sub:function(x,y) { return this.addGroup.invOp(x,y) },
  div:function(x,y) { return this.mulGroup.invOp(x,y) }, // 這個用箭頭函數會錯， Why ?
  power:function(x,n) { return this.mulGroup.power(x,n) },
  init:function(addGroup, mulGroup) {
    this.addGroup = addGroup;
    this.mulGroup = mulGroup;
    this.zero = addGroup.e;
    this.add  = function(x,y) { return this.addGroup.op(x,y) }
    this.neg  = function(x) { return this.addGroup.inv(x) }
    this.one  = mulGroup.e;
    this.mul  = function(x,y) { return this.mulGroup.op(x,y) }
    this.inv  = function(x) { return this.mulGroup.inv(x) }
    this.power= function(x,n) { return this.mulGroup.power(x,n) }
  }
}

// ref : https://en.wikipedia.org/wiki/Group_homomorphism
//  https://en.wikipedia.org/wiki/Fundamental_theorem_on_homomorphisms
// 同態：h(a • b) = h(a) x h(b) 
F.homomorphism=function(h, g1, g2) {
  var a=g1.random(), b=g2.random();
  return eq(h(group1.op(a,b)), group2.op(h(a), h(b)))
}

// ref : https://en.wikipedia.org/wiki/Isomorphism
//  https://en.wikipedia.org/wiki/Isomorphism_theorem
// 同構：h(a • b) = h(a) • h(b)
R.isomorphism=function(h1, h2, g1, g2) {
  var a1=g1.random(), b1=g2.random();
  var a2=g1.random(), b2=g2.random();
  return homorphism(h1,g1,g2)&&homorphism(h2,g2,g1);
}

// ========== Float Field =================
var FloatAddGroup={
  e:0,
  op:function(x,y) { return x+y },
  inv:function(x) { return -x},
}

extend(FloatAddGroup, Group, S.Float);

var FloatMulGroup={
  e:1,
  op:function(x,y) { return x*y },
  inv:function(x) { return 1/x},
}

extend(FloatMulGroup, Group, S.Float);

var FloatField=extend({}, Field, S.Float);

FloatField.init(FloatAddGroup, FloatMulGroup);

// ========== Finite Field =================
var FiniteAddGroup={
  e:0,
  op:function(x,y) { return (x+y)%this.n },
  inv:function(x) { return (this.n-x) }
}

extend(FiniteAddGroup, Group, S.Float);

var FiniteMulGroup={
  e:1,
  op:function(x,y) { return (x*y)%this.n }, 
  inv:function(x) { return this.invMap[x] },
  setOrder:function(n) {
    this.n = n;
    let invMap = new Map();
    for (var x=1; x<n; x++) {
      var y = this.op(x,x);
      invMap.set(x,y);
    }
    this.invMap = invMap;
  }
}

extend(FiniteMulGroup, Group, S.Float);

var FiniteField=extend({}, Field, S.Float);

FiniteField.create=function(n) {
  var F = extend({}, FiniteField);
  var addGroup = extend({n:n}, FiniteAddGroup);
  var mulGroup = extend({n:n}, FiniteMulGroup);
  F.init(addGroup, mulGroup);
  mulGroup.setOrder(n);
  return F;
}

class MathObj {
  constructor() {}
  str() { return this.toString() }
  eq(y) { return eq(this, y) }
}

// =========== Field Object ==============
class FieldObj extends MathObj {
  constructor(field) { 
    super();
    this.field = field;
  }
  
  add(y) { return this.field.add(this,y) }
  mul(y) { return this.field.mul(this,y) }
  neg() { return this.field.neg(this) }
  inv() { return this.field.inv(this) }
  div(y) { return this.field.div(this,y) }
  sub(y) { return this.field.sub(this,y) }
  power(n) { return this.field.power(this,n) }
}

// =========== Complex Field ==============
var ComplexField=extend({}, Field);

class Complex extends FieldObj {
  constructor(a,b) {
    super(ComplexField);
    this.a = a; this.b = b; 
  }
  conj() { return new Complex(this.a, -1*this.b); }
  
  toString() { 
    var op = (this.b<0)?'':'+';
    return this.a+op+this.b+'i'; 
  }
  
  parse(s) {
    var m = s.match(/^([^\+]*)(\+(.*))?$/);
    var a = parseFloat(m[1]);
    var b = typeof m[3]==='undefined'?1:parseFloat(m[3]);
    return new Complex(a, b)
  }
  
  ln() {
    var a=this.a, b=this.b, r=a*a+b*b;
    var w = 1/2*Math.log(r);
    var x = Math.acos(a/Math.sqrt(r));
    return new Complex(w, x);
  }
  
  exp() {
    var a=this.a, b=this.b;
    var r=Math.exp(a);
    return new Complex(r*Math.cos(b), r*Math.sin(b));
  }
}

var ComplexSet=S.create(function(a) { return a.constructor===Complex });

var ComplexAddGroup={
  e:new Complex(0,0),
  op:function(x,y) { return new Complex(x.a+y.a, x.b+y.b) },
  inv:function(x) { return new Complex(-x.a, -x.b) }
}

extend(ComplexAddGroup, Group, ComplexSet);

var ComplexMulGroup={
  e:new Complex(1,0),
  op:function(x,y) {
    return new Complex(x.a*y.a-x.b*y.b, x.a*y.b+x.b*y.a);
  },
  inv:function(x) {
    var a=x.a,b=x.b, r=(a*a+b*b);
    return new Complex(a/r, -b/r);
  } 
}

extend(ComplexMulGroup, Group, ComplexSet);

extend(ComplexField, ComplexSet);

ComplexField.init(ComplexAddGroup, ComplexMulGroup);

// =========== Ratio Field ==============
var RatioField=extend({}, Field);

class Ratio extends FieldObj {
  constructor(a,b) {
    super(RatioField);
    this.a = a; this.b = b; 
  }

  reduce() {
    var a = this.a, b=this.b;
    var c = I.gcd(a, b);
    return new Ratio(a/c, b/c);
  }
  
  toString() { return this.a+'/'+this.b; }

  parse(s) {
    var m = s.match(/^(\d+)(\/(\d+))?$/);
    var a = parseInt(m[1]);
    var b = typeof m[3]==='undefined'?1:parseInt(m[3]);
    return new Ratio(a, b)
  } 
}

var RatioAddGroup={
  e:new Ratio(0,1),
  op:function(x,y) { return new Ratio(x.a*y.b+x.b*y.a, x.b*y.b) },
  inv:function(x) { return new Ratio(-x.a, x.b); },
}
  
extend(RatioAddGroup, Group);

var RatioMulGroup={
  e:new Ratio(1,1),
  op:function(x,y) { return new Ratio(x.a*y.a, x.b*y.b) },
  inv:function(x) { return new Ratio(x.b, x.a) },
}

extend(RatioMulGroup, Group);

RatioField.init(RatioAddGroup, RatioMulGroup);

// export 
module.exports = {
  Group:Group,
  Field:Field,
  FloatField:FloatField,
  FiniteField:FiniteField,
  ComplexField:ComplexField,
  Complex:Complex,
  RatioField:RatioField,
  Ratio:Ratio,
};


