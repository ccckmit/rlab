// ========== Group =================
class Group {
  constructor(op) { 
    for (var key in op) {
      this[key] = op[key];
    }
  }
  invOp(x,y) { return this.op(x,this.inv(y)); }
  power(x,n) {
    var p=this.identity;
    for (var i=0;i<n;i++) {
      p=this.op(p,x);
    }
    return p;
  }
}

// ========== Field =================
class Field {
  constructor(addGroup, mulGroup) {
    this.addGroup = addGroup;
    this.mulGroup = mulGroup;
    this.zero=addGroup.identity;
    this.one=mulGroup.identity;
    this.add=addGroup.op;
    this.mul=mulGroup.op;
    this.neg=addGroup.inv;
    this.inv=mulGroup.inv;
  }
  sub(x,y) { return this.add(x, this.neg(y)) }
  div(x,y) { return this.mul(x, this.inv(y)) }
  power(x,n) { return this.mulGroup.power(x,n); }
}

// ========== Float Field =================
var FloatAddGroup = new Group({
  identity:0,
  op:(x,y)=>x+y,
  inv:(x)=>-x,
});

var FloatMulGroup = new Group({
  identity:1,
  op:(x,y)=>x*y,
  inv:(x)=>1/x,
});

var FloatField = new Field(FloatAddGroup, FloatMulGroup);

// ========== Finite Field =================
class FiniteAddGroup extends Group{
  constructor(n) {
    super({
      identity:0,
      op:(x,y)=>(x+y)%this.n,
      inv:(x)=>(x)=>(this.n-x),
    });
    this.n = n;
  }
}

class FiniteMulGroup extends Group {
  constructor(n) {
    super({
      identity:1,
      op:(x,y)=>(x*y)%n,
      inv:(x)=>this.invMap[x],
    });
    this.setOrder(n);
  }
  
  setOrder(n) {
    this.n = n;
    let invMap = new Map();
    for (var x=1; x<n; x++) {
      var y = this.op(x,x);
      invMap.set(x,y);
    }
    this.invMap = invMap;
  }
}

// ========== Finite Field =================
class FiniteField extends Field {
  constructor(n) {
    super(new FiniteAddGroup(n), new FiniteMulGroup(n));
  } 
}

// =========== Object => Group ==============
class AddGroup extends Group {
	op(x,y) { return x.add(y) }
	inv(x) { return x.neg() }
}

class MulGroup extends Group {
	op(x,y) { return x.mul(y) }
	inv(x) { return x.inv() }
}

class FieldObj {
	constructor() {}
	toString() { return this.str() }
}

// ========== Complex Number =================
class Complex extends FieldObj {
  constructor(a,b) { super(); this.a = a; this.b = b; }
	
  conj() { return new Complex(this.a, -1*this.b); }
  
  neg() { return new Complex(-this.a, -this.b); }
  
  add(c2) { return new Complex(this.a+c2.a, this.b+c2.b); }
  
  sub(c2) { return new Complex(this.a-c2.a, this.b-c2.b); }
  
  mul(c2) {
    var a=this.a, b=this.b, c=c2.a, d=c2.b;
    return new Complex(a*c-b*d, a*d+b*c);
  }
  
  inv() {
    var c=this.a, d=this.b;
    var r=(c*c+d*d);
    return new Complex(c/r, -d/r); 
  }
  
  div(c2) {
    var a=this.a, b=this.b, c=c2.a, d=c2.b;
    var r=(c*c+d*d);
    return new Complex((a*c+b*d)/r, (b*c-a*d)/r);
  }
  
  str() { 
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

class Ratio extends FieldObj {
  constructor(a,b) { super(); this.a = a; this.b = b; }
  
  div(r2) { return new Ratio(this.a*r2.b, this.b*r2.a); }
	
  mul(r2) { return new Ratio(this.a*r2.a, this.b*r2.b); }
 
  inv() { return new Ratio(this.b, this.a); }
  
  neg() { return new Ratio(-this.a, this.b); }
	
  add(r2) {
	  if (this.b === r2.b)
			return new Ratio(this.a+r2.a, this.b);
	  return new Ratio(this.a*r2.b+this.b*r2.a, this.b*r2.b); 
	}
	
  sub(r2) { 
	  return this.add(r2.neg());
	}
  
  toString() { return this.a+'/'+this.b; }

  parse(s) {
    var m = s.match(/^(\d+)(\/(\d+))?$/);
    var a = parseInt(m[1]);
    var b = typeof m[3]==='undefined'?1:parseInt(m[3]);
    return new Ratio(a, b)
  } 
}

// ========== Complex Field =================
var ComplexField = new Field(
	new AddGroup({identity:new Complex(0,0)}), 
	new MulGroup({identity:new Complex(1,0)})
);

	// ========== Rational Field =================
var RatioField = new Field(
	new AddGroup({identity:new Ratio(0,1)}), 
	new MulGroup({identity:new Ratio(1,1)})
);

module.exports = {
  Group:Group,
  Field:Field,
  FloatField:FloatField,
  FiniteField:FiniteField,
  ComplexField:ComplexField,
	RatioField:RatioField,
  Complex:Complex,
	Ratio:Ratio,
};

/*
// =========== Object => Group ==============
class AddGroup extends Group {
	op(x,y) { return x.add(y) }
	inv(x) { return x.neg() }
}

class MulGroup extends Group {
	op(x,y) { return x.mul(y) }
	inv(x) { return x.inv() }
}

// ========== Complex Number =================
class Complex {
  constructor(a,b) { this.a = a; this.b = b; }
  
  conj() { return new Complex(this.a, -1*this.b); }
  
  neg() { return new Complex(-this.a, -this.b); }
  
  add(c2) { return new Complex(this.a+c2.a, this.b+c2.b); }
  
  sub(c2) { return new Complex(this.a-c2.a, this.b-c2.b); }
  
  mul(c2) {
    var a=this.a, b=this.b, c=c2.a, d=c2.b;
    return new Complex(a*c-b*d, a*d+b*c);
  }
  
  inv() {
    var c=this.a, d=this.b;
    var r=(c*c+d*d);
    return new Complex(c/r, -d/r); 
  }
  
  div(c2) {
    var a=this.a, b=this.b, c=c2.a, d=c2.b;
    var r=(c*c+d*d);
    return new Complex((a*c+b*d)/r, (b*c-a*d)/r);
  }
  
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

class Ratio {
  constructor(a,b) { this.a = a; this.b = b; }
  
  
  div(r2) { return new Ratio(this.a*r2.b, this.b*r2.a); }
   mul(r2) { return new Ratio(this.a*r2.a, this.b*r2.b); }
 
  inv() { return new Ratio(this.b, this.a); }
  
  neg() { return new Ratio(-this.a, this.b); }
	
  add(r2) {
	  if (this.b === r2.b)
			return new Ratio(this.a+r2.a, this.b);
	  return new Ratio(this.a*r2.b+this.b*r2.a, this.b*r2.b); 
	}
	
  sub(r2) { 
	  return this.add(r2.neg());
//	  return new Ratio(this.a*r2.b-this.b*r2.a, this.b*r2.b); 
	}
  
  toString() { return this.a+'/'+this.b; }

  parse(s) {
    var m = s.match(/^(\d+)(\/(\d+))?$/);
    var a = parseInt(m[1]);
    var b = typeof m[3]==='undefined'?1:parseInt(m[3]);
    return new Ratio(a, b)
  } 
}
class RatioField extends Field {
  constructor() {
    super(
		  new Group({
        identity:{a:0,b:0},
				op :function(x,y) { 
				  return (x.b===y.b)?{a:x.a+y.a, b:x.b}
					                  :{a:x.a*y.b+x.b*y.a, b:x.b*y.b};
				},
				inv:function(x) { return {a:-x.a,b:x.b} },
		  }), 
		  new Group({
				identity:{a:1,b:1},
				op :function(x,y) { 
				  return (x.b===y.a)?{a:x.a,b:y.b}
					                  :{a:x.a*y.a, b:x.b*y.b};
				},
				inv:function(x) { return {a:x.b, b:x.a} },
			})
		)
  }
}


*/
