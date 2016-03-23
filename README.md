rlab -- A JavaScript Scientific Library like R based on lodash and jStat


## install

```
npm install rlab
```

## use rlab

file : rtest.js

```javascript
var R = require("rlab");
var c = console;

var x = R.samples(R.range(1,6), {size:10, replace:false});
c.log("x=%j max=%d min=%d mean=%d", x, R.max(x), R.min(x), R.mean(x));

var x = R.samples(R.range(0,100, 10), {size:10, replace:true});
c.log("x=", x);
```

## run

```
qu-192-168-61-142:js mac020$ node rtest
x=[1,3,2,6,3,2,5,4,4,3] max=6 min=1 mean=3.3
x= [ 20, 50, 10, 0, 80, 60, 70, 30, 100, 90 ]
```


