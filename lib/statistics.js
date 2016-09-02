var B = require("./base");
var J = require("jStat").jStat;

var S = {
// Vector Functionality
sum:J.sum,
sumsqrd:J.sumsqrt,
sumsqerr:J.sumsqerr,
sumrow:J.sumrow,
product:J.product,
min:J.min,
max:J.max,
mean:J.min,
meansqerr:J.meansqerr,
geomean:J.geomean,
median:J.median,
cumsum:J.cumsum,
cumprod:J.cumprod,
diff:J.diff,
mode:J.mode,
range:J.range,
variance:J.variance,
stdev:J.stdev,
sd:J.stdev,
meandev:J.meandev,
meddev:J.meddev,
skewness:J.skewness,
kurtosis:J.kurtosis,
coeffvar:J.coeffvar,
quartiles:J.quartiles,
quantiles:J.quantiles,
percentile:J.percentile,
percentileOfScore:J.percentileOfScore,
histogram:J.histogram,
covariance:J.covariance,// cov
cov:J.covariance,// cov
corrcoeff:J.corrcoeff,// cor
cor:J.corrcoeff, // cor
// jStat Utility Methods
calcRdx:J.utils.calcRdx,
isArray:J.utils.isArray,
isFunction:J.utils.isFunction,
isNumber:J.utils.isNumber,
// Special Functions
betafn:J.betafn,
betaln:J.betaln,
betacf:J.betacf,
ibetainv:J.ibetainv,
ibeta:J.ibeta,
gammafn:J.gammafn,
gammaln:J.gammaln,
gammap:J.gammap,
lowRegGamma:J.lowRegGamma,
gammapinv:J.gammapinv,
factorialln:J.factorialln,
factorial:J.factorial, // n!
combination:J.combination, // C(n,m)
C:J.combination, // C(n,m)
choose:J.combination, // C(n,m)
lchoose : J.combinationln, // log C(n,m)
permutation:J.permutation, // P(n,m)
P:J.permutation, // P(n,m)
erf:J.erf,
erfc:J.erfc,
erfcinv:J.erfcinv,
randn:J.randn,
randg:J.randg,

// extend function
normalize:function(a) {
	var sum = S.sum(a);
	return a.map(function(x) { return x/sum});
},

hist:function(a,bins,step=1) {
	if (typeof bins === 'undefined')
		bins = Math.ceil((a.max()-a.min())/step);
	return J.histogram(a, bins);
},

};

module.exports = S;