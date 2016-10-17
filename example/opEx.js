var R = require("../rlab");
var O = R.Op;
var op = R.Op.op;
var p  = R.Op.parse;
var c1 = p('5+2i');

print('isField(%s)=%s', c1, O.isField(c1));
print('%s instanceof Complex=%s', c1, c1 instanceof O.Field.Complex);
print('op(+, 3, 5)=%s', op('+', 3, 5));
print('op(+, 3, 5+2i)=%s', op('+', 3, c1));
var v1 = [1,2,3], v2=[2,3,c1];
print('op(+, %s, %s)=%s', v1, v2, op('+', v1, v2));
// print('op(dot, %s, %s)=%s', v1, v2, op('dot', v1, v2));
print('%s.sqrt()=%s', v1, op('sqrt', v1).str());
print('%s.sqrt()=%s', v2, op('sqrt', v2).str());
var exp = (x)=>3*x*x+5;
print('%s.eval(%d)=%s', exp, 1, op('eval', 1, exp).str());
print('eval(1,%s)=', exp, op('eval', 1, exp).str());
print('eval(%s,1)=', exp, op('eval', exp, 1).str());
print('eval(%s,%s)=', v1, exp, op('eval', v1, exp).str());
print('%s.eval(%s)=', v1.eval(exp));

var a = 8, b=3;
print('op(%, %s, %s)=%s', a, b, op('%',a,b));
print('(%s).mod(%s)=%s', a, b, a.mod(b));

// var exp=(x)=>x;
print('(%s).diff(1)=', exp, exp.diff(1));
print('(%s).integral(0,1)=', exp, exp.integral(0,1));

