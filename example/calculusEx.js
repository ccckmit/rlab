var R = require("../rlab");
var x2=(x)=>x*x;

with (R) {
print('d(x^2,2)=', differential(x2, 2));
print('d(sin(x/4),pi/4)=', differential(sin, PI/4));
print('i(x^2,0,1)=', integral(x2,0,1));
print('i(sin(x),0,pi/2)=', integral(sin,0,PI/2));	
}
