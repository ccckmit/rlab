var R = require("../rlab");
var v = [1,3,5];
log("v.max()=", v.max());
log("v.min()=", v.min());
log("v.sum()=", v.sum());
log("v.normalize()=", v.normalize());
log("v.normalize().sum()=", v.normalize().sum());
log("v.product()=", v.product());
log("v.mean()=", v.mean());
log("v.range()=", v.range());
log("v.median()=", v.median());
log("v.variance()=", v.variance());
log("v.sd()=", v.sd(), " sd^2=", v.sd()*v.sd());
log("v.cov(v)=", v.cov(v), "v.cor(v)=", v.cor(v));
log("factorial(5)=", R.factorial(5));
