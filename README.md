rlab -- A JavaScript Scientific Library like R based on lodash and jStat


## install

```
npm install rlab
```

## use rlab

file : rtest.js

```javascript
var R = require("rlab");
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
```

## run

```
D:\Dropbox\github\rlab>node rtest
dice= [ 1, 2, 3, 4, 5, 6 ]
x= [ 2, 1, 3, 4, 6, 5 ]
chain1:x= [ 6, 2, 3, 5, 1, 4 ]
chain2:x= [1, 2, 5, 6, 6, 3, 6, 6, 2, 5]
x= [ 5, 6, 3, 4, 6, 3, 4, 2, 5, 2 ] max= 6 min= 2 mean= 4 sd= 1.4907119849998598

cov(x,x)= 1.4907119849998598
cor(x,x)= 1
factorial(10)= 3628800
lfactorial(10)= 15.1044
choose(5,2)= 10
lchoose(5,2)= 2.302585092994045
permutation(5,2)= 20
runif(10, -5, -1)= [-3.3, -2.68, -3.5, -2.96, -4.48, -1.9, -2.12, -2.02, -4.59,
-4.09]
dunif(-3, -5, -1)= 0.25
punif(-3, -5, -1)= 0.5
qunif(0.5, -5, -1)= -3
x= [0.79, 0.49, 1.01, -1.13, 0.19, 0.4, -0.14, 1.01, 0.1, -1]
x.sort()= [-0.14, -1, -1.13, 0.1, 0.19, 0.4, 0.49, 0.79, 1.01, 1.01]
rbinom(10, 5, 0.5)= [ 3, 2, 3, 2, 3, 2, 3, 1, 2, 1 ]
dbinom(4, 5, 0.5)= 0.15625
dbinom(5, 5, 0.5)= 0.03125
pbinom(4, 5, 0.5)= 0.96875
qbinom(0.9, 5, 0.5)= 4
=========== report ==========
name    : "ttest(X)"
h       : "H0:mu=0"
alpha   : 0.05
op      : "="
pvalue  : 0.49
ci      : [-0.37, 0.71]
df      : 9
mean    : 0.17
sd      : 0.75
A= [[1, 2, 3], [4, 5, 6], [7, 3, 9]]
iA= [[-0.9, 0.3, 0.1], [-0.2, 0.4, -0.2], [0.77, -0.37, 0.1]]
AiA= [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
====iA=====
 [[     -0.9,      0.3,      0.1],
 [     -0.2,      0.4,     -0.2],
 [     0.77,    -0.37,      0.1]]
```


