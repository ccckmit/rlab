rlab -- A JavaScript Scientific Library like R based on lodash and jStat


## install

```
npm install rlab
```

## use rlab

file : rtest.js

```javascript
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
```

## run

```
D:\Dropbox\github\rlab>node rtest
x= [ 10, 50, 20, 0, 80, 100, 90, 30, 40, 60 ]
x=[5,5,5,4,1,3,4,3,2,1] max=5 min=1 mean=3.3
cov(x,x)= 1.57
cor(x,x)= 1.00
factorial(10)= 3628800
lfactorial(10)= 15.10
choose(5,2)= 10
lchoose(5,2)= 2.30
permutation(5,2)= 20
runif(10, -5, -1)= [-3.16, -2.57, -2.12, -4.13, -2.3, -2.14, -2.55, -4.78, -1.77
, -3.66]
dunif(-3, -5, -1)= 0.25
punif(-3, -5, -1)= 0.5
qunif(0.5, -5, -1)= -3
x.str()= [-0.15, 0.2, 0.53, -0.45, -0.34, -0.98, 1.25, 0.44, 1.61, 0.11]
x.sortBy()= [-0.98, -0.45, -0.34, -0.15, 0.11, 0.2, 0.44, 0.53, 1.25, 1.61]
str(x)= [-0.15, 0.2, 0.53, -0.45, -0.34, -0.98, 1.25, 0.44, 1.61, 0.11]
rbinom(10, 5, 0.5)= [ 2, 3, 2, 3, 1, 5, 2, 2, 1, 2 ]
dbinom(4, 5, 0.5)= 0.15625
dbinom(5, 5, 0.5)= 0.03125
pbinom(4, 5, 0.5)= 0.96875
qbinom(0.9, 5, 0.5)= 4
sd(x)= 0.78
=========== report ==========
name    : "ttest(X)"
h       : "H0:mu=0"
alpha   : 0.05
op      : "="
pvalue  : 0.39
ci      : [-0.34, 0.78]
df      : 9
mean    : 0.22
sd      : 0.78
A= [[1, 2, 3], [4, 5, 6], [7, 3, 9]]
iA= [[-0.9, 0.3, 0.1], [-0.2, 0.4, -0.2], [0.77, -0.37, 0.1]]
AiA= [[1, 0, -0.], [0, 1, -0.], [0, 0., 1.]]

```


