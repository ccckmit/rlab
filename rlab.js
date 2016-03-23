var R=require("lodash");
var J=require("J");

R.repeats=function(f, x, n) {
	var a=[];
	for (var i=0; i<n; i++) {
		a.push(f(x));
	}
	return a;
}

R.samples=function(space, arg) {
	var arg = R.defaults(arg, {size:1, replace:true});
	if (arg.replace)
		return R.sampleSize(space, arg.size);
	else
		return R.repeats(R.sample, space, arg.size);
}

R._range = R.range;
R.range=function(start, end, step) {
    step = step || 1;
    return R._range(start,end+0.00001*step, step); 
}

R.sd = function(a, flag) { 
  flag = flag || 1;
  return J.stdev(a, flag); 
}

// 協方差
R.cov = function(x, y) { return J.stdev(x, y); }
// 相關係數
R.cor = function(x, y) { return J.corrcoeff(x, y); }
// 階層 n!
R.factorial = function(n) { return J.factorial(n); }
// log(n!)
R.lfactorial = function(n) { return J.factorialln(n); }
// 組合 C(n,m)
R.choose = function(n,m) { return J.combination(n, m); }
// log C(n,m)
R.lchoose = function(n,m) { return J.combinationln(n, m); }
// 組合 C(n,m)
R.permutation = function(n,m) { return J.permutation(n, m); }
// log C(n,m)
R.lchoose = function(n,m) { return J.combinationln(n, m); }

module.exports = R;
