rlab -- A JavaScript Scientific Library like R

## introduction

The rlab is a A JavaScript Scientific Library like R.

It's based on `lodash.js , jStat.js and numeric.js`

## install

```
npm install rlab
```

## use rlab

file : probabilityEx.js

```javascript
var R = require("rlab");
var c = console;
var dice = R.steps(1,6);
c.log("sample(1:6, 10)", R.samples(dice, 10));
c.log("runif(10,0,1)=", R.runif(10, 0, 1).str());
c.log("rnorm(10,5,1)=", R.rnorm(10, 5, 1).str());
c.log("dnorm(5,5,1)=", R.dnorm(5, 5, 1));
c.log("pnorm(5,5,1)=", R.pnorm(5, 5, 1));
c.log("qnorm(0.5,5,1)=", R.qnorm(0.5, 5, 1));
c.log("rbinom(10, 5, 0.5)=", R.rbinom(10,5,0.5));
c.log("dbinom(4, 5, 0.5)=", R.dbinom(4,5,0.5));
c.log("dbinom(5, 5, 0.5)=", R.dbinom(5,5,0.5));
c.log("pbinom(4, 5, 0.5)=", R.pbinom(4,5,0.5));
c.log("qbinom(0.9, 5, 0.5)=", R.qbinom(0.9,5,0.5));
```

run :

```
$ node probabilityEx.js
sample(1:6, 10) [ 3, 5, 3, 2, 3, 3, 1, 2, 4, 3 ]
runif(10,0,1)= [0.9119,0.5899,0.6839,0.1350,0.6894,0.9512,0.8186,0.5826,0.4279,0
.5125]
rnorm(10,5,1)= [5.8961,5.4312,6.0002,5.3623,5.5281,4.4413,6.2144,5.7173,5.3111,1
.3146]
dnorm(5,5,1)= 0.3989422804014327
pnorm(5,5,1)= 0.5
qnorm(0.5,5,1)= 5
rbinom(10, 5, 0.5)= [ 2, 1, 2, 2, 4, 4, 1, 4, 3, 2 ]
dbinom(4, 5, 0.5)= 0.15625
dbinom(5, 5, 0.5)= 0.03125
pbinom(4, 5, 0.5)= 0.96875
qbinom(0.9, 5, 0.5)= 4
```

file : statisticsEx.js

```javascript
var R = require("rlab");
var c = console;
var v = [1,3,5];
c.log("v.max()=", v.max());
c.log("v.min()=", v.min());
c.log("v.sum()=", v.sum());
c.log("v.normalize()=", v.normalize());
c.log("v.normalize().sum()=", v.normalize().sum());
c.log("v.product()=", v.product());
c.log("v.mean()=", v.mean());
c.log("v.range()=", v.range());
c.log("v.unique()=", v.unique());
c.log("v.median()=", v.median());
c.log("v.variance()=", v.variance());
c.log("v.deviation()=", v.deviation());
c.log("v.sd()=", v.sd(), " sd^2=", v.sd()*v.sd());
c.log("v.cov(v)=", v.cov(v), "v.cor(v)=", v.cor(v));
c.log("factorial(5)=", R.factorial(5));
```

run : 

```
$ node statisticsEx.js
v.max()= 5
v.min()= 1
v.sum()= 9
v.normalize()= [ 0.1111111111111111, 0.3333333333333333, 0.5555555555555556 ]
v.normalize().sum()= 1
v.product()= 15
v.mean()= 3
v.range()= 4
v.unique()= [ 1, 3, 5 ]
v.median()= 3
v.variance()= 2.6666666666666665
v.deviation()= [ -2, 0, 2 ]
v.sd()= 1.632993161855452  sd^2= 2.6666666666666665
v.cov(v)= 2 v.cor(v)= 1
factorial(5)= 120
```

file : matrixEx.js

```javascript
var M = require("rlab").M;

var c = console;
var v = [1,2,3];
c.log("v.sin()=", v.sin());
c.log("v.norm2()=", v.norm2());
c.log("v.norm2Squared()=", v.norm2Squared());

var A = [[1,2,3],[4,5,6],[7,3,9]];
var AiA = A.inv().dot(A);
c.log("AiA=\n", AiA.strM());
c.log("AiA.tr()=\n", AiA.tr().strM());
c.log("A=\n", A.str());
c.log("A.mul(0.1)=\n", A.mul(0.1).strM());
c.log("A.row(1)=", A.row(1));
c.log("A.col(1)=", A.col(1));
c.log("A.sumM()=", A.sumM());
c.log("A.rowSum(2)=", A.rowSum(2));
c.log("A.colSum(2)=", A.colSum(2));
c.log("A.mean(row)=", A.rowMean().str());
c.log("A.mean(col)=", A.colMean().str());

var D = M.diag(v);
c.log("D=", D);

var Eλ = M.eigR(A);
var E = Eλ.E, λ=Eλ.lambda;
c.log("E*[λ]*E-1=", E.dot(λ.diag()).dot(E.inv()).strM());
```

run : 

```
$ node matrixEx.js
v.sin()= [ 0.8414709848078965, 0.9092974268256817, 0.1411200080598672 ]
v.norm2()= 3.7416573867739413
v.norm2Squared()= 14
AiA=
 [[          1,   1.11e-16,  -1.11e-16],
 [          0,          1,  4.441e-16],
 [ -3.331e-16, -3.331e-16,          1]]
AiA.tr()=
 [[          1,          0, -3.331e-16],
 [   1.11e-16,          1, -3.331e-16],
 [  -1.11e-16,  4.441e-16,          1]]
A=
 [[1.0000,2.0000,3.0000],[4.0000,5.0000,6.0000],[7.0000,3.0000,9.0000]]
A.mul(0.1)=
 [[        0.1,        0.2,        0.3],
 [        0.4,        0.5,        0.6],
 [        0.7,        0.3,        0.9]]
A.row(1)= [ 4, 5, 6 ]
A.col(1)= [ 2, 5, 3 ]
A.sumM()= 40
A.rowSum(2)= [ 6, 15, 19 ]
A.colSum(2)= [ 12, 10, 18 ]
A.mean(row)= [2.0000,5.0000,6.3333]
A.mean(col)= [4.0000,3.3333,6.0000]
D= [ [ 1, 0, 0 ], [ 0, 2, 0 ], [ 0, 0, 3 ] ]
E*[λ]*E-1= [[          1,          2,          3],
 [          4,          5,          6],
 [          7,          3,          9]]
```

## IDE

There is a webIDE server based on rlab , you may start it by :

```
node rserver.js
```

## Author

Author: ccckmit

Email : ccckmit@gmail.com 


