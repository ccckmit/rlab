var R  = require("../rlab");
var M  = R.Math;
// var A  = R.Algebra;
var op = R.Math.op;
var p  = R.Math.parse;
var c1 = p('5+2i');

print('isField(%s)=%s', c1, M.isField(c1));
print('%s instanceof Complex=%s', c1, c1 instanceof Complex);
print('op(+, 3, 5)=%s', op('+', 3, 5));
print('op(+, 3, 5+2i)=%s', op('+', 3, c1));
var v1 = [1,2,3], v2=[2,3,c1], m=[[1,2], [3,p('1+2i')]];
print('op(+, %s, %s)=%s', v1, v2, op('+', v1, v2));
// print('op(dot, %s, %s)=%s', v1, v2, op('dot', v1, v2));
print('%s.sqrt()=%s', v1, op('sqrt', v1).str());
print('%s.sqrt()=%s', v2, op('sqrt', v2).str());
print('v2.sum()=%s', op('sum', v2).str())
print('v2.product()=%s', op('product', v2).str())
print('m.sum()=%s', op('sum', m).str())
print('m.product()=%s', op('product', m).str())
print('v1.max()=%s', op('max', v1).str())
var a = 8, b=3;
print('op(%, %s, %s)=%s', a, b, op('%',a,b));
print('(%s).mod(%s)=%s', a, b, a.mod(b));


/*
var F2 = (x,y,z)=>[3*x*x*z, 4*x*y*z, y*z*z]; // F2=[3x^2z, 4xyz, yz^2]
var F2a = M.fa(F2); // curl(F2) = [z^2-4xy, 3x^2, 4yz]
print('curl(%s, [1,0,0])=%s', F1, M.fcurl(F2a, [0,1,1])); // = [1,0,4]
*/
/*
var R = require("../rlab");
var A = R.Algebra;
var P = A.Polynomial;
var p1 = new P([1,2]), p2=new P([1,1,1]);
console.log('p1=%s p2=%s', p1, p2);
console.log('p1+p2=%s', p1.add(p2));
console.log('p1-p2=%s', p1.sub(p2));
console.log('p1*p2=%s', p1.mul(p2));
console.log('p1(3)=', p1.eval(3));
console.log('p2(2)=', p2.eval(2));
*/