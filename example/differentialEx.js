var R = require("../rlab");

var d = R.D.d, i=R.D.i, sin=R.sin, PI = R.PI, x2=(x)=>x*x;

log('d(x^2,2)=', d(x2, 2));
log('d(sin(x/4),pi/4)=', d(sin, PI/4));
log('i(x^2,0,1)=', i(x2,0,1));
log('i(sin(x),0,pi/2)=', i(sin,0,PI/2));
