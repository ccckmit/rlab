var R = require("../rlab").Math;

var v1 = [1,2,3];
var exp = (x)=>3*x*x+5;
print('%s.eval(%d)=%s', exp, 1, exp.eval(1).str());
print('%s.eval(%s)=%s', v1, exp, v1.eval(exp).str());

print('norm(%s)=%s', v1, R.norm(v1));

print('(%s).diff(1)=', exp, exp.diff(1));
print('(%s).integral(0,1)=', exp, exp.integral(0,1));

var f1 = (x,y)=>2*x+y;

print('pdiff(%s, [1,0], 0)=', f1, R.pdiff(f1, [1,0], 0));
print('pdiff(%s, [0,1], 1)=', f1, R.pdiff(f1, [0,1], 1));
print('pdiff(%s, [1,1], 1)=', f1, R.pdiff(f1, [1,1], 1));

var f2 = (x,y)=>x*y-x*x;
print('grad(%s, [1,0], 0)=', f2, R.fgrad(f2, [1,0], 0));
/*
var F1 = (x,y,z)=>[3*x*x*x*z,4*x*y*z,y*z*z]; // F1=[3x^3z,4xyz,yz^2]
print('div(%s, [1,0,0])=%s', F1, R.fdiv(F1, [0,1,1])); // div(F1) = 9x^2+4xz+2yz = 2
*/