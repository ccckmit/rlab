var R=require("lodash");

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

module.exports = R;
