var R = require("../rlab");
var A = R.Algebra;
var Complex=R.Field.Complex;
var F = A.FunctionRing;
var be = R.Rule.be;
var times3=new A.FunctionObj(function(x) { 
  return x.mul(3);
})
var times2=new A.FunctionObj(function(x) { 
  return x*2;
})
var timesC3=new A.FunctionObj(function(x) { 
  return x.mul(new Complex(3,0));
})

var power2=new A.FunctionObj(function(x) { 
  return x.power(2);
})
var p2=power2.eval(2);
print('p2=', p2.toString());
var p32=F.add(times3,power2);
print('p32(2)=', p32.eval(2));

print('t3(p2)=', times3.eval(2));
print('t2(p2)=', times2.eval(2));
print('t3(3+2i)=', timesC3.eval(new Complex(3,2)).toString());
// var t32=FF.add(times3,times2);
// print('t32(2)=', t32.eval(2));
