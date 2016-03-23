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

var x = R.samples(R.range(0,100, 10), {size:10, replace:true});
c.log("x=", x);

var x = R.samples(R.range(1,6), {size:10, replace:false});
c.log("x=%j max=%d min=%d mean=%d", x, R.max(x), R.min(x), R.mean(x));

c.log("cov(x,x)=", R.cov(x,x));
c.log("cor(x,x)=", R.cor(x,x)); // 相關係數
c.log("factorial(10)=", R.factorial(10)); // 階層 n!
c.log("lfactorial(10)=", R.lfactorial(10)); // log(n!)
c.log("choose(5,2)=", R.choose(5,2)); // 組合 C(n,m)
c.log("lchoose(5,2)=", R.lchoose(5,2)); // log C(n,m)
c.log("permutation(5,2)=", R.permutation(5,2)); // P(n,m)

c.log("runif(10, -5, -1)=", R.runif(10, -5, -1)); 
c.log("rnorm(10, 5, 2)=", R.rnorm(10, 5, 2)); 
```

## run

```
$ node rtest
x=[1,3,2,6,3,2,5,4,4,3] max=6 min=1 mean=3.3
x= [ 20, 50, 10, 0, 80, 60, 70, 30, 100, 90 ]
```


