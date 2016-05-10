var B = require("./base");
var J = require("jStat").jStat;

var S = {
// 均等分布 : jStat.uniform( a, b )
runif:function(n, a, b) { return B.calls(n, J.uniform.sample, a, b); },
dunif:J.uniform.pdf,
punif:J.uniform.cdf,
qunif:J.uniform.inv,
// 常態分布 : jStat.normal( mean, std )
rnorm:function(n, mean, sd) { return B.calls(n, J.normal.sample, mean, sd); },
dnorm:J.normal.pdf,
pnorm:J.normal.cdf,
qnorm:J.normal.inv,
// 布瓦松分布 : jStat.poisson
rpois:function(n, l) { return B.calls(n, J.poisson.sample, l); },
dpois:J.poisson.pdf,
ppois:J.poisson.cdf,
qpois:J.poisson.inv,
// F 分布 : jStat.centralF( df1, df2 )
rf:function(n, df1, df2) { return B.calls(n, J.centralF.sample, df1, df2); },
df:J.centralF.pdf,
pf:J.centralF.cdf,
qf:J.centralF.inv,
// T 分布 : jStat.studentt( dof )
rt:function(n, dof) { return B.calls(n, J.studentt.sample, dof); },
dt:J.studentt.pdf,
pt:J.studentt.cdf,
qt:J.studentt.inv,
// Beta 分布 : jStat.beta( alpha, beta )
rbeta:function(n, alpha, beta) { return B.calls(n, J.beta.sample, alpha, beta); },
beta:J.beta.pdf,
pbeta:J.beta.cdf,
qbeta:J.beta.inv,
// 柯西分布 : jStat.cauchy( local, scale )
rcauchy:function(n, local, scale) { return B.calls(n, J.cauchy.sample, local, scale); },
dcauchy:J.cauchy.pdf,
pcauchy:J.cauchy.cdf,
qcauchy:J.cauchy.inv,
// chisquare 分布 : jStat.chisquare( dof )
rchisq:function(n, dof) { return B.calls(n, J.chisquare.sample, dof); },
dchisq:J.chisquare.pdf,
pchisq:J.chisquare.cdf,
qchisq:J.chisquare.inv,
// 指數分布 : jStat.exponential( rate )
rexp:function(n, rate) { return B.calls(n, J.exponential.sample, rate); },
dexp:J.exponential.pdf,
pexp:J.exponential.cdf,
qexp:J.exponential.inv,
// Gamma 分布 : jStat.gamma( shape, scale )
rgamma:function(n, shape, scale) { return B.calls(n, J.gamma.sample, shape, scale); },
dgamma:J.gamma.pdf,
pgamma:J.gamma.cdf,
qgamma:J.gamma.inv,
// 反 Gamma 分布 : jStat.invgamma( shape, scale )
rinvgamma:function(n, shape, scale) { return B.calls(n, J.invgamma.sample, shape, scale); },
dinvgamma:J.invgamma.pdf,
pinvgamma:J.invgamma.cdf,
qinvgamma:J.invgamma.inv,
// 對數常態分布 : jStat.lognormal( mu, sigma )
rlognormal:function(n, mu, sigma) { return B.calls(n, J.dlognormal.sample, mu, sigma); },
lognormal:J.lognormal.pdf,
plognormal:J.lognormal.cdf,
qlognormal:J.lognormal.inv,
// Pareto 分布 : jStat.pareto( scale, shape )
rpareto:function(n, scale, shape) { return B.calls(n, J.pareto.sample, scale, shape); },
dpareto:J.pareto.pdf,
ppareto:J.pareto.cdf,
qpareto:J.pareto.inv,
// Weibull 分布
rweibull:function(n, scale, shape) { return B.calls(n, J.weibull.sample, scale, shape); },
dweibull:J.weibull.pdf,
pweibull:J.weibull.cdf,
qweibull:J.weibull.inv,
// 三角分布 : jStat.triangular
rtriangular:function(n, a, b, c) { return B.calls(n, J.triangular.sample, a, b, c); },
dtriangular:J.triangular.pdf,
ptriangular:J.triangular.cdf,
qtriangular:J.triangular.inv,
// 類似 Beta 分布，但計算更簡單 : jStat.kumaraswamy( alpha, beta )
rkumaraswamy:function(n, alpha, beta) { return B.calls(n, J.kumaraswamy.sample, alpha, beta); },
dkumaraswamy:J.kumaraswamy.pdf,
pkumaraswamy:J.kumaraswamy.cdf,
qkumaraswamy:J.kumaraswamy.inv,

// ========== 離散分佈的 r, q 函數 ============
qcdf:function(cdf, q, N, p) {
  for (var i=0; i<=N; i++) {
    if (cdf(i, N, p) > q) return i;
  }
  return N;
},
rcdf:function(cdf, n, N, p) {
  var a = [];
  for (var i=0; i<n; i++) {
    var q = Math.random();
    a.push(cdf(q, N, p));
  }
  return a;
},
// 二項分布 : jStat.binomial
dbinom:J.binomial.pdf,
pbinom:J.binomial.cdf,
qbinom:function(q, N, p) { return S.qcdf(S.pbinom, q, N, p); },
rbinom:function(n, N, p) { return S.rcdf(S.qbinom, n, N, p); },
// 負二項分布 : jStat.negbin
dnbinom:J.negbin.pdf,
pnbinom:J.negbin.cdf,
qnbinom:function(q, N, p) { return S.qcdf(S.pnbinom, q, N, p); },
rnbinom:function(n, N, p) { return S.rcdf(S.qnbinom, n, N, p); },
// 超幾何分布 : jStat.hypgeom
dhyper:J.hypgeom.pdf,
phyper:J.hypgeom.cdf,
qhyper:function(q, N, m, n) { return S.qcdf(S.phyper, q, N, p); },
rhyper:function(n, N, m, k) { return S.rcdf(S.qhyper, n, N, p); },

};

module.exports = S;