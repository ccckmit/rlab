var R = require("../rlab");
var op = R.op;
var p  = R.parse;
var c1 = p('5+2i');

print('=========== Operator ===============');
print('isField(%s)=%s', c1, R.isField(c1));
print('%s instanceof Complex=%s', c1, c1 instanceof R.Complex);
print('op(+, 3, 5)=%s', op('+', 3, 5));
print('op(+, 3, 5+2i)=%s', op('+', 3, c1));
var v1 = [1,3,2], v2=[2,3,c1], m=[[1,2], [3,p('1+2i')]];
print('op(+, %s, %s)=%s', v1, v2, op('+', v1, v2));
// print('op(dot, %s, %s)=%s', v1, v2, op('dot', v1, v2));
print('%s.sqrt()=%s', v1, op('sqrt', v1).str());
print('%s.sqrt()=%s', v2, op('sqrt', v2).str());
print('v2.sum()=%s', op('sum', v2).str())
print('v2.product()=%s', op('product', v2).str())
print('m.sum()=%s', op('sum', m).str())
print('m.product()=%s', op('product', m).str())
print('v1.max()=%s', op('max', v1).str())
print('v1.min()=%s', op('min', v1).str())
// print('v1.min=', v1.min);

var exp = (x)=>3*x*x+5;
print('%s.eval(%d)=%s', exp, 1, op('eval', 1, exp).str());
print('eval(1,%s)=', exp, op('eval', 1, exp).str());
print('eval(%s,1)=', exp, op('eval', exp, 1).str());
print('eval(%s,%s)=', v1, exp, op('eval', v1, exp).str());
print('%s.eval(%s)=', v1.eval(exp));

print('norm(%s)=%s', v1, R.norm(v1));

var a = 8, b=3;
print('op(%, %s, %s)=%s', a, b, op('%',a,b));
print('(%s).mod(%s)=%s', a, b, a.mod(b));
var P = R.Polynomial;
var p1 = new P([1,2]), p2=new P([2,2,4]), p3=new P([1,-6,21,-52]);
print('=========== Polynomial ===============');
print('p1=%s p2=%s', p1, p2);
print('p1+p2=%s', p1.add(p2));
print('p1-p2=%s', p1.sub(p2));
print('p1*p2=%s', p1.mul(p2));
print('p1(3)=', p1.eval(3));
print('p2(2)=', p2.eval(2));
print('p2.root()=%s', p2.root());
print('p3.root()=%s', p3.root());