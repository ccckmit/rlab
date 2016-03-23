var R = require("./rlab");
var c = console;

var x = R.samples(R.range(1,6), {size:10, replace:false});
c.log("x=%j max=%d min=%d mean=%d", x, R.max(x), R.min(x), R.mean(x));

var x = R.samples(R.range(0,100, 10), {size:10, replace:true});
c.log("x=", x);