var R = require("../rlab");
var A = R.A;
// var print = console.log;
// var A = require("../lib/algebra");
var F = A.FloatField;
print('F:2+3=',F.add(2,3));
print('F:2^3=',F.power(2,3));

var F7= new A.FiniteField(7);
print('F7:3+6=', F7.add(3,6));
print('F7:3*6=', F7.mul(3,6));

var C = A.ComplexField;
var c1 = new A.Complex(2,3);
print('C:c1=%s', c1);
print('C:c1+c1=%s', C.add(c1,c1));
print('C:c1*c1=%s', C.mul(c1,c1));
print('C:(c1*c1)/c1=%s', C.div(C.mul(c1,c1),c1));

var co = new A.Complex(2,3);
print('C:co=%s', co);
print('C:co+co=%s', co.add(co));
print('C:co*co=%s', co.mul(co));
print('C:co/co=%s', co.div(co));
print('C:co*co/co=%s', co.mul(co).div(co));

var Q = A.RatioField;
var q1 = new A.Ratio(2,3);
print('Q:q1=%s', q1);
print('Q:q1+q1=%s', Q.add(q1,q1).toString());
print('Q:q1-q1=%s', Q.sub(q1,q1));
print('Q:q1*q1=%s', Q.mul(q1,q1));
print('Q:q1/q1*q1=%s', Q.mul(Q.div(q1,q1),q1));
