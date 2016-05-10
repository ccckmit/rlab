var R = require("../rlab");
var v = [1,3,5];

var x = R.rnorm(10, 0, 0.1);
log("x=", x.str());
log("x.sort()=", x.sort().str());

var t1=R.ttest({x:x, mu:0});
R.report(t1);

