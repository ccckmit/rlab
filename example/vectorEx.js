var R = require("../lib/vector");
var V = R.V, p = R.parse;
var c1 = p('5+2i');
var v1 = [1,3,2], v2=[2,3,c1], m=[[1,2], [3,p('1+2i')]];

print('=========== Operator ===============');
print('add(3,5)=%s', V.add(3, 5));
print('oadd(3,5+2i)=%s', V.oadd(c1, 3));
print('v1=%s v2=%s m=%s', v1, v2, m);
print('add(%s, %s)=%s', v1, v2, V.add(v1, v2));
print('%s.sqrt()=%s', v1, V.sqrt(v1));
print('%s.osqrt()=%s', v2, V.osqrt(v2));
print('v1.sum()=%s', V.sum(v1))
print('v2.osum()=%s', V.osum(v2))
print('v1.product()=%s', V.product(v1))
print('m.osum()=%s', V.osum(m))
print('m.oproduct()=%s', V.oproduct(m))
print('v1.max()=%s', V.max(v1))
print('v1.min()=%s', V.min(v1))
print('=========== Map Reduce===============');
print('map1(v1,x^2)=%s', V.map1(v1, (x)=>x*x));
print('map2(v1,v1,x+y)=%s', V.map2(v1, v1, (x,y)=>x+y));
print('reduce(v1,add,0)=%s', V.reduce(v1, (x,y)=>x+y, 0));
print('max(v1)=%s', V.max(v1));
print('min(v1)=%s', V.min(v1));
print('add(v1,v2)=%s', V.add(v1,v1));
print('sin(v1)=%s', V.sin(v1));
print('=========== Eval ===============');
var exp = (x)=>3*x*x+5;
print('%s.map1(%s)=%s', v1, exp, V.map1(v1, exp));
print('norm(%s)=%d', v1, V.norm(v1));
print('=========== repeat ===============');
print('repeat([3,2,2],0)=%s\n', V.repeat([3,2,2]));
print('random([2,2],[0,1])=%s\n', V.repeat([2,2],()=>R.random(0,1)));

