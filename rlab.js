var R=require("lodash");
var J=require("jStat").jStat;

R.repeats=function(f, x, n) {
	var a=[];
	for (var i=0; i<n; i++) {
		a.push(f(x));
	}
	return a;
}

R.calls=function() {
  var args = Array.prototype.slice.call(arguments);
  var n = args[0];
  var f = args[1];
  var params = args.slice(2, args.length);
  var a=[];
  for (var i=0; i<n; i++)
    a.push(f.apply(null, params));
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


// 均等分布
R.runif=function(n, a, b) { return R.calls(n, J.uniform.sample, a, b); }
R.dunif=function(x, a, b) { return J.uniform.pdf(x, a, b); }
R.punif=function(q, a, b) { return J.uniform.cdf(q, a, b); }
R.qunif=function(p, a, b) { return J.uniform.inv(p, a, b); }
// 常態分布
R.rnorm=function(n, mean, sd) { return R.calls(n, J.normal.sample, mean, sd); }
R.dnorm=function(x, mean, sd) { return J.normal.pdf(x, mean, sd); }
R.pnorm=function(q, mean, sd) { return J.normal.cdf(q, mean, sd); }
// R.qnorm=function(p, mean, sd) { return R.q2x(p, function (q) { return R.pnorm(q, mean, sd);});}
R.qnorm=function(p, mean, sd) { return J.normal.inv(p, mean, sd); }
// 布瓦松分布
R.rpois=function(n, l) { return R.calls(n, J.poisson.sample, l); }
R.dpois=function(x, l) { return J.poisson.pdf(x, l); }
R.ppois=function(q, l) { return J.poisson.cdf(q, l); }
R.qpois=function(p, l) { return J.poisson.inv(p, l); }

// F 分布
R.rf=function(n, df1, df2) { return R.calls(n, J.centralF.sample, df1, df2); }
R.df=function(x, df1, df2) { return J.centralF.pdf(x, df1, df2); }
R.pf=function(q, df1, df2) { return J.centralF.cdf(q, df1, df2); }
R.qf=function(p, df1, df2) { return J.centralF.inv(p, df1, df2); }
// T 分布
R.rt=function(n, dof) { return R.calls(n, J.studentt.sample, dof); }
R.dt=function(x, dof) { return J.studentt.pdf(x, dof); }
R.pt=function(q, dof) { return J.studentt.cdf(q, dof); }
R.qt=function(p, dof) { return J.studentt.inv(p, dof); }
// Beta 分布
R.rbeta=function(n, alpha, beta) { return R.calls(n, J.beta.sample, alpha, beta); }
R.dbeta=function(x, alpha, beta) { return J.beta.pdf(x, alpha, beta); }
R.pbeta=function(q, alpha, beta) { return J.beta.cdf(q, alpha, beta); }
R.qbeta=function(p, alpha, beta) { return J.beta.inv(p, alpha, beta); }
// 柯西分布
R.rcauchy=function(n, local, scale) { return R.calls(n, J.cauchy.sample, local, scale); }
R.dcauchy=function(x, local, scale) { return J.cauchy.pdf(x, local, scale); }
R.pcauchy=function(q, local, scale) { return J.cauchy.cdf(q, local, scale); }
R.qcauchy=function(p, local, scale) { return J.cauchy.inv(p, local, scale); }
// chisquare 分布
R.rchisq=function(n, dof) { return R.calls(n, J.chisquare.sample, dof); }
R.dchisq=function(x, dof) { return J.chisquare.pdf(x, dof); }
R.pchisq=function(q, dof) { return J.chisquare.cdf(q, dof); }
R.qchisq=function(p, dof) { return J.chisquare.inv(p, dof); }
// 指數分布
R.rexp=function(n, rate) { return R.calls(n, J.exponential.sample, rate); }
R.dexp=function(x, rate) { return J.exponential.pdf(x, rate); }
R.pexp=function(q, rate) { return J.exponential.cdf(q, rate); }
R.qexp=function(p, rate) { return J.exponential.inv(p, rate); }
// Gamma 分布
R.rgamma=function(n, shape, scale) { return R.calls(n, J.gamma.sample, shape, scale); }
R.dgamma=function(x, shape, scale) { return J.gamma.pdf(x, shape, scale); }
R.pgamma=function(q, shape, scale) { return J.gamma.cdf(q, shape, scale); }
R.qgamma=function(p, shape, scale) { return J.gamma.inv(p, shape, scale); }
// 反 Gamma 分布
R.rinvgamma=function(n, shape, scale) { return R.calls(n, J.invgamma.sample, shape, scale); }
R.dinvgamma=function(x, shape, scale) { return J.invgamma.pdf(x, shape, scale); }
R.pinvgamma=function(q, shape, scale) { return J.invgamma.cdf(q, shape, scale); }
R.qinvgamma=function(p, shape, scale) { return J.invgamma.inv(p, shape, scale); }
// 對數常態分布
R.rlognormal=function(n, mu, sigma) { return R.calls(n, J.lognormal.sample, mu, sigma); }
R.dlognormal=function(x, mu, sigma) { return J.lognormal.pdf(x, mu, sigma); }
R.plognormal=function(q, mu, sigma) { return J.lognormal.cdf(q, mu, sigma); }
R.qlognormal=function(p, mu, sigma) { return J.lognormal.inv(p, mu, sigma); }
// Pareto 分布
R.rpareto=function(n, scale, shape) { return R.calls(n, J.pareto.sample, scale, shape); }
R.dpareto=function(x, scale, shape) { return J.pareto.pdf(x, scale, shape); }
R.ppareto=function(q, scale, shape) { return J.pareto.cdf(q, scale, shape); }
R.qpareto=function(p, scale, shape) { return J.pareto.inv(p, scale, shape); }
// Weibull 分布
R.rweibull=function(n, scale, shape) { return R.calls(n, J.weibull.sample, scale, shape); }
R.dweibull=function(x, scale, shape) { return J.weibull.pdf(x, scale, shape); }
R.pweibull=function(q, scale, shape) { return J.weibull.cdf(q, scale, shape); }
R.qweibull=function(p, scale, shape) { return J.weibull.inv(p, scale, shape); }
// 三角分布
R.rtriangular=function(n, a, b, c) { return R.calls(n, J.triangular.sample, a, b, c); }
R.dtriangular=function(x, a, b, c) { return J.triangular.pdf(x, a, b, c); }
R.ptriangular=function(q, a, b, c) { return J.triangular.cdf(q, a, b, c); }
R.qtriangular=function(p, a, b, c) { return J.triangular.inv(p, a, b, c); }
// 類似 Beta 分布，但計算更簡單
R.rkumaraswamy=function(n, alpha, beta) { return R.calls(n, J.kumaraswamy.sample, alpha, beta); }
R.dkumaraswamy=function(x, alpha, beta) { return J.kumaraswamy.pdf(x, alpha, beta); }
R.pkumaraswamy=function(q, alpha, beta) { return J.kumaraswamy.cdf(q, alpha, beta); }
R.qkumaraswamy=function(p, alpha, beta) { return J.kumaraswamy.inv(p, alpha, beta); }

module.exports = R;
