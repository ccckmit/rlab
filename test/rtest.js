var R = require("../rlab");
var M = R.M;
var c = console;

var dice = R.steps(1,6);
c.log("dice=", dice);

var x = R.samples(dice, 6, {replace:false});
c.log("x=", x);

var x = R(dice).samples(6, {replace:false}).value();
c.log("chain1:x=", x);

var x = R(dice).samples(10).str();
c.log("chain2:x=", x);

var x = R.samples(dice, 10);
c.log("x=", x, "max=", R.max(x), "min=", R.min(x), 
      "mean=", R.mean(x), "sd=", R.sd(x));
			
c.log("cov(x,x)=", R.cov(x,x));
c.log("cor(x,x)=", R.cor(x,x)); // 相關係數
c.log("factorial(10)=", R.factorial(10)); // 階層 n!
c.log("lfactorial(10)=", R.lfactorial(10).toFixed(4)); // log(n!)
c.log("choose(5,2)=", R.choose(5,2)); // 組合 C(n,m)
c.log("lchoose(5,2)=", R.lchoose(5,2)); // log C(n,m)
c.log("permutation(5,2)=", R.permutation(5,2)); // P(n,m)
c.log("runif(10, -5, -1)=", R.runif(10, -5, -1).str()); 
// c.log(".chain(10).runif(-5,-1)=", R.chain(10).runif(-5,-1)); 
c.log("dunif(-3, -5, -1)=", R.dunif(-3, -5, -1)); 
c.log("punif(-3, -5, -1)=", R.punif(-3, -5, -1)); 
c.log("qunif(0.5, -5, -1)=", R.qunif(0.5, -5, -1)); 

var x = R.rnorm(10, 0, 1);
c.log("x=", R.str(x));
c.log("x.sort()=", x.sort().str());

c.log("rbinom(10, 5, 0.5)=", R.rbinom(10,5,0.5));
c.log("dbinom(4, 5, 0.5)=", R.dbinom(4,5,0.5));
c.log("dbinom(5, 5, 0.5)=", R.dbinom(5,5,0.5));
c.log("pbinom(4, 5, 0.5)=", R.pbinom(4,5,0.5));
c.log("qbinom(0.9, 5, 0.5)=", R.qbinom(0.9,5,0.5));

var t1=R.ttest({x:x, mu:0} );
R.report(t1);

var A = [[1,2,3],[4,5,6],[7,3,9]];
var iA = M.inv(A);
c.log("A=", R.str(A));
c.log("iA=", R.str(iA));
var AiA = M.dot(A, iA);

c.log("AiA=", R.str(AiA));

c.log("====iA=====\n", M.str(iA))

