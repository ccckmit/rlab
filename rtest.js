var R = require("./rlab");
var c = console;

var x = R.samples(R.steps(0,100, 10), 10, {replace:true}).value();
c.log("x=", x);

var x = R.samples(R.steps(1,6), 10, {replace:false}).value();
c.log("x=%j max=%d min=%d mean=%d", x, R.max(x), R.min(x), R.mean(x));

c.log("cov(x,x)=", R.cov(x,x).toFixed(2));
c.log("cor(x,x)=", R.cor(x,x).toFixed(2)); // 相關係數
c.log("factorial(10)=", R.factorial(10)); // 階層 n!
c.log("lfactorial(10)=", R.lfactorial(10).toFixed(2)); // log(n!)
c.log("choose(5,2)=", R.choose(5,2)); // 組合 C(n,m)
c.log("lchoose(5,2)=", R.lchoose(5,2).toFixed(2)); // log C(n,m)
c.log("permutation(5,2)=", R.permutation(5,2)); // P(n,m)

c.log("runif(10, -5, -1)=", R.runif(10, -5, -1).str()); 
c.log("dunif(-3, -5, -1)=", R.dunif(-3, -5, -1)); 
c.log("punif(-3, -5, -1)=", R.punif(-3, -5, -1)); 
c.log("qunif(0.5, -5, -1)=", R.qunif(0.5, -5, -1)); 

var x = R.rnorm(10, 0, 1);
c.log("x.str()=", x.str());
c.log("x.sortBy()=", x.sortBy().str());
c.log("str(x)=", R.str(x.value()));

c.log("rbinom(10, 5, 0.5)=", R.rbinom(10,5,0.5));
c.log("dbinom(4, 5, 0.5)=", R.dbinom(4,5,0.5));
c.log("dbinom(5, 5, 0.5)=", R.dbinom(5,5,0.5));
c.log("pbinom(4, 5, 0.5)=", R.pbinom(4,5,0.5));
c.log("qbinom(0.9, 5, 0.5)=", R.qbinom(0.9,5,0.5));

c.log("sd(x)=", x.sd().toFixed(2));

var t1=R.ttest({x:x.value(), mu:0} );
R.report(t1);

var A = [[1,2,3],[4,5,6],[7,3,9]];
var iA = R.M.inv(A);
c.log("A=", R.str(A));
c.log("iA=", R.str(iA));
var AiA = R.M.dot(A, iA);
c.log("AiA=", R.str(AiA));




