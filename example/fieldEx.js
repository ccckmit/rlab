var R = require("../rlab");
var F = R.Field, be=R.Rule.be, p=R.p;
/*
var FF = F.FloatField;
var be = R.Rule.be, eq = R.Rule.eq;
be('FF:2+3=5',FF.add(2,3)===5);
be('FF:2^3=8',FF.power(2,3)===8);
*/
var F7= F.FiniteField.create(7);
var a = 3, b=6;
be('F7:3+6=2', F7.add(a,b)===2);
be('F7:3*6=4', F7.mul(a,b)===4);
be('F7:closability(+)=> a+b in F7', F7.addGroup.closability(3,6));
be('F7:associativity(+)=>(a+b)+c=a+(b+c)', F7.addGroup.associativity(3,6,4));
be('F7:identity(+)=>a+0=a', F7.addGroup.identity(3));
be('F7:inversability(+)=a-a=0', F7.addGroup.inversability(3));

var C = F.ComplexField;
var c1 = new F.Complex(2,3);
be('C :c1==2+3i', c1.str()==='2+3i');
be('C :c1+c1=4+6i', C.add(c1,c1).str()==='4+6i');
be('C :c1*c1=-5+12i', C.mul(c1,c1).str()==='-5+12i');
be('C :(c1*c1)/c1=2+3i', C.div(C.mul(c1,c1),c1).str()==='2+3i');

var Q = F.RatioField;
var q1 = new F.Ratio(2,3);
be('Q:q1=2/3', q1.str()==='2/3');
be('Q:q1+q1=4/3', Q.add(q1,q1).reduce().str()==='4/3');
be('Q:q1-q1=0', Q.sub(q1,q1).a===0);
be('Q:q1*q1=4/9', Q.mul(q1,q1).str()==='4/9');
be('Q:q1/q1*q1=2/3', Q.mul(Q.div(q1,q1),q1).reduce().str()==='2/3');
be('Q:q1^3=8/27', q1.power(3).str()==='8/27');

var co = new F.Complex(2,3);
be('C:co=2+3i', co.str()==='2+3i');
be('C:co+co=4+6i', co.add(co).str()==='4+6i');
be('C:co*co=-5+12i', co.mul(co).str()==='-5+12i');
be('C:co/co=1+0i', co.div(co).str()==='1+0i');
be('C:co*co/co=2+3i', co.mul(co).div(co).str()==='2+3i');

var c1=p('1+2i'), c2=p('2+1i'), c3=p('10+0i');
print('%s * %s=%s', c1, c2, c1.mul(c2));
print('(%s)*3=%s', c1, c1.mul(3));
print('3 * 3=%s', p('3').mul(3));
// print('c1.exp().log()=', c1.exp().log().str());

var sqrt2 = Math.sqrt(2);
var c=new F.Complex(sqrt2, sqrt2);
print('c=%s', c);
print('c.toPolar=%j', c.toPolar().str());
print('c*c=%s', c.mul(c));
print('c^2=%s', c.power(2));
print('c^2.sqrt()=%s', c.power(2).sqrt());
// be('c1.exp().ln()=', c1.exp().ln().str()==='1+2i');