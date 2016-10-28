var B = require('./base');
var J = require('jStat').jStat;
var S = require('./calculus');
var ncall = B.ncall;
// var T, R;

// ========== 離散分佈的 r, q 函數 ============
S.qcdf=function(cdf, q, N, p) {
  for (var i=0; i<=N; i++) {
    if (cdf(i, N, p) > q) return i;
  }
  return N;
}

S.rcdf=function(cdf, n, N, p) {
  var a = [];
  for (var i=0; i<n; i++) {
    var q = Math.random();
    a.push(cdf(q, N, p));
  }
  return a;
}

S.EPSILON=0.0000000001;
// 均等分布 : jStat.uniform( a, b )
S.dunif=(x,a=0,b=1)=>J.uniform.pdf(x,a,b);
S.punif=(q,a=0,b=1)=>J.uniform.cdf(q,a,b);
S.qunif=(p,a=0,b=1)=>J.uniform.inv(p,a,b);
S.runif=(n,a=0,b=1)=>ncall(n, J.uniform, 'sample', a, b);
// 常態分布 : jStat.normal( mean, sd )
S.dnorm=(x,mean=0,sd=1)=>J.normal.pdf(x,mean,sd);
S.pnorm=(q,mean=0,sd=1)=>J.normal.cdf(q,mean,sd);
S.qnorm=(p,mean=0,sd=1)=>J.normal.inv(p,mean,sd);
S.rnorm=(n,mean=0,sd=1)=>ncall(n, J.normal, 'sample', mean, sd);
// F 分布 : jStat.centralF( df1, df2 )
S.df=(x,df1,df2)=>J.centralF.pdf(x,df1,df2);
S.pf=(q,df1,df2)=>J.centralF.cdf(q,df1,df2);
S.qf=(p,df1,df2)=>J.centralF.inv(p,df1,df2);
S.rf=(n,df1,df2)=>ncall(n, J.centralF, 'sample', df1, df2);
// T 分布 : jStat.studentt( dof )
S.dt=(x,dof)=>J.studentt.pdf(x,dof);
S.pt=(q,dof)=>J.studentt.cdf(q,dof);
S.qt=(p,dof)=>J.studentt.inv(p,dof);
S.rt=(n,dof)=>ncall(n, J.studentt, 'sample', dof);
// Beta 分布 : jStat.beta( alpha, beta )
S.dbeta=(x,alpha,beta)=>J.beta.pdf(x,alpha,beta);
S.pbeta=(q,alpha,beta)=>J.beta.cdf(q,alpha,beta);
S.qbeta=(p,alpha,beta)=>J.beta.inv(p,alpha,beta);
S.rbeta=(n,alpha,beta)=>ncalls(n, J.beta, 'sample', alpha, beta);
// 柯西分布 : jStat.cauchy( local, scale )
S.dcauchy=(x,local,scale)=>J.cauchy.pdf(x,local,scale);
S.pcauchy=(q,local,scale)=>J.cauchy.cdf(q,local,scale);
S.qcauchy=(p,local,scale)=>J.cauchy.inv(q,local,scale);
S.rcauchy=(n,local,scale)=>ncall(n, J.cauchy, 'sample', local, scale);
// chisquare 分布 : jStat.chisquare( dof )
S.dchisq=(x,dof)=>J.chisquare.pdf(x,dof);
S.pchisq=(q,dof)=>J.chisquare.cdf(q,dof);
S.qchisq=(p,dof)=>J.chisquare.inv(p,dof);
S.rchisq=(n,dof)=>ncall(n, J.chisquare, 'sample', dof);
// 指數分布 : jStat.exponential( rate )
S.dexp=(x,rate)=>J.exponential.pdf(x,rate);
S.pexp=(q,rate)=>J.exponential.cdf(q,rate);
S.qexp=(p,rate)=>J.exponential.inv(p,rate);
S.rexp=(n,rate)=>ncall(n, J.exponential, 'sample', rate);
// Gamma 分布 : jStat.gamma( shape, scale )
S.dgamma=(x,shape,scale)=>J.gamma.pdf(x,shape,scale);
S.pgamma=(q,shape,scale)=>J.gamma.cdf(q,shape,scale);
S.qgamma=(p,shape,scale)=>J.gamma.inv(p,shape,scale);
S.rgamma=(n,shape,scale)=>ncall(n, J.gamma, 'sample', shape, scale);
// 反 Gamma 分布 : jStat.invgamma( shape, scale )
S.rinvgamma=(n,shape,scale)=>ncall(n, J.invgamma, 'sample', shape, scale);
S.dinvgamma=(x,shape,scale)=>J.invgamma.pdf(x,shape,scale);
S.pinvgamma=(q,shape,scale)=>J.invgamma.cdf(q,shape,scale);
S.qinvgamma=(p,shape,scale)=>J.invgamma.inv(p,shape,scale);
// 對數常態分布 : jStat.lognormal( mu, sigma )
S.dlognormal=(n, mu, sigma)=>J.lognormal.pdf(x,sigma);
S.plognormal=(n, mu, sigma)=>J.lognormal.cdf(q,sigma);
S.qlognormal=(n, mu, sigma)=>J.lognormal.inv(p,sigma);
S.rlognormal=(n, mu, sigma)=>ncall(n, J.dlognormal, 'sample', mu, sigma);
// Pareto 分布 : jStat.pareto( scale, shape )
S.dpareto=(n, scale, shape)=>J.pareto.pdf(x,scale,shape);
S.ppareto=(n, scale, shape)=>J.pareto.cdf(q,scale,shape);
S.qpareto=(n, scale, shape)=>J.pareto.inv(p,scale,shape);
S.rpareto=(n, scale, shape)=>ncall(n, J.pareto, 'sample', scale, shape);
// Weibull 分布 jStat.weibull(scale, shape)
S.dweibull=(n, scale, shape)=>J.weibull.pdf(x,scale,shape);
S.pweibull=(n, scale, shape)=>J.weibull.cdf(q,scale,shape);
S.qweibull=(n, scale, shape)=>J.weibull.inv(p,scale,shape);
S.rweibull=(n, scale, shape)=>ncall(n, J.weibull, 'sample', scale, shape);
// 三角分布 : jStat.triangular(a, b, c)
S.dtriangular=(n, a, b, c)=>J.triangular.pdf(x,a,b,c);
S.ptriangular=(n, a, b, c)=>J.triangular.cdf(q,a,b,c);
S.qtriangular=(n, a, b, c)=>J.triangular.inv(p,a,b,c);
S.rtriangular=(n, a, b, c)=>ncall(n, J.triangular, 'sample', a, b, c);
// 類似 Beta 分布，但計算更簡單 : jStat.kumaraswamy(alpha, beta)
S.dkumaraswamy=(n, alpha, beta)=>J.kumaraswamy.pdf(x,alpha,beta);
S.pkumaraswamy=(n, alpha, beta)=>J.kumaraswamy.cdf(q,alpha,beta);
S.qkumaraswamy=(n, alpha, beta)=>J.kumaraswamy.inv(p,alpha,beta);
S.rkumaraswamy=(n, alpha, beta)=>ncalls(n, J.kumaraswamy, 'sample', alpha, beta);

// ========== 離散分佈的 r, q 函數 ============
// 二項分布 : jStat.binomial(n, p0)
S.dbinom=(x, size, prob)=>J.binomial.pdf(x, size, prob);
S.pbinom=(q, size, prob)=>J.binomial.cdf(q, size, prob);
S.qbinom=(p, size, prob)=>S.qcdf(S.pbinom, p, size, prob);
S.rbinom=(n, size, prob)=>S.rcdf(S.qbinom, n, size, prob);
// 負二項分布 : jStat.negbin(r, p)
S.dnbinom=(x, size, prob)=>J.negbin.pdf(x, size, prob);
S.pnbinom=(q, size, prob)=>J.negbin.cdf(q, size, prob);
S.qnbinom=(p, size, prob)=>S.qcdf(S.pnbinom, p, size, prob);
S.rnbinom=(n, size, prob)=>S.rcdf(S.qnbinom, n, size, prob);
// 超幾何分布 : jStat.hypgeom(N, m, n)
S.dhyper=(x, m, n, k)=>J.hypgeom.pdf(k, m, n, k);
S.phyper=(q, m, n, k)=>J.hypgeom.cdf(q, m, n, k);
S.qhyper=(p, m, n, k)=>S.qcdf(S.phyper, p, m, n, k);
S.rhyper=(nn,m, n, k)=>S.rcdf(S.qhyper, nn, m, n, k);
// 布瓦松分布 : jStat.poisson(l)
S.dpois=(x, lambda)=>J.poisson.pdf(x, lambda);
S.ppois=(q, lambda)=>J.poisson.cdf(q, lambda);
S.qpois=(p, lambda)=>S.qcdf(S.ppois, p, lambda);
S.rpois=(n, lambda)=>S.rcdf(S.qpois, n, lambda);

// ====================== statistics =================================
// extend function

S.normalize=function(a) {
	var sum = S.sum(a);
	return a.map(function(x) { return x/sum});
}

// Vector Functionality
B.copyFunctions(S, J, "sumsqrt,sumsqerr,sumrow,mean,meansqerr,geomean,median,cumsum,cumprod,diff,mode,range,variance,stdev,meandev,meddev,skewness,kurtosis,coeffvar,quartiles,quantiles,percentile,percentileOfScore,histogram,covariance,corrcoeff,calcRdx,betafn,betacf,ibetainv,ibeta,gammafn,gammaln,gammap,lowRegGamma,gammapinv,factorialln,factorial,combination,combinationln,permutation,erf,erfc,erfcinv,randn,randg".split(","));

B.mapFunctions(S, J, {
	C:'combination',// C(n,m)
	choose:'combination',// C(n,m)
	lchoose:'combinationln',// log C(n,m)
	P:'permutation', // P(n,m)
	sd:'stdev',
	cov:'covariance',
	cor:'corrcoeff',
});

// =============== 檢定 ==============================
B.mix(S, B);

var T = S;

T.test = function(o) { // name, D, x, mu, sd, y, alpha, op
  Object.assign(o, {alpha:0.05, op:"="});
  var alpha = o.alpha;
  var pvalue, interval;
  var D      = o.D;
  var q1     = D.o2q(o); // 單尾檢定的 pvalue
  
  if (o.op === "=") {
    if (q1>0.5) q1 = 1-q1; // (q1>0.5) 取右尾，否則取左尾。
    pvalue= 2*q1; // 對稱情況：雙尾檢定的 p 值是單尾的兩倍。
    interval = [D.q2p(alpha/2, o, "L"), D.q2p(1-alpha/2, o, "R")];
  } else {
    if (o.op === "<") { // 右尾檢定 H0: q < 1-alpha, 
      interval = [ D.q2p(alpha, o, "L"), Infinity ]; 
      pvalue = 1-q1;
    }
    if (o.op === ">") { // 左尾檢定 H0: q > alpha
      interval=[-Infinity, D.q2p(1-alpha, o, "R")];
      pvalue = q1;
    }
  }
  return { 
    name: o.name,
    h: D.h(o),
    alpha: alpha,
    op: o.op, 
    pvalue: pvalue, 
    ci : interval, 
    df : D.df(o),
    report: function() { S.report(this) }
  };
}

T.report = function(o) {
  console.log("=========== report ==========");
  for (var k in o) {
    if (typeof o[k] !== "function")
      console.log(k+"\t: "+S.str(o[k]));
  }
}

var t1 = { // 單樣本 T 檢定 t = (X-mu)/(S/sqrt(n))
  h:function(o) { return "H0:mu"+o.op+o.mu; }, 
  o2q:function(o) {
    var x = o.x, n = x.length;
    var t = (S.mean(x)-o.mu)/(S.sd(x)/Math.sqrt(n)); 
    return S.pt(t, n-1);
  },
  // P(X-mu/(S/sqrt(n))<t) = q ; 信賴區間 P(T<q)
  // P(mu > X-t*S/sqrt(n)) = q ; 這反而成了右尾檢定，所以左尾與右尾確實會反過來
  q2p:function(q, o) {
    var x = o.x, n = x.length;
    return S.mean(x) + S.qt(q, n-1) * S.sd(x) / Math.sqrt(n);
  },
  df:function(o) { return o.x.length-1; }
}

var t2vareq = { // σ1=σ2, 合併 T 檢定 (雙樣本)
  h:function(o) { return "H0:mu1"+o.op+"mu2" }, 
  // S^2 = (n1-1)*S1^2+(n2-1)*S2^2)/(n1-1+n2-1)
  sd:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1= S.sd(x), S2 = S.sd(y);
    var S = Math.sqrt(((n1-1)*S1*S1+(n2-1)*S2*S2)/(n1-1+n2-1)); 
    return S;
  },
  // T = ((X-Y)-(mu1-mu2))/(sqrt(1/n1+1/n2)*S)
  o2q:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S = this.sd(o);
    var t = (S.mean(x)-S.mean(y)-o.mu)/(Math.sqrt(1/n1+1/n2)*S);
    return S.pt(t, n1+n2-2);
  },
  // t=((X-Y)-mu)/(sqrt(1/n1+1/n2)*S), (X-Y)-t*sqrt(1/n1+1/n2)*S = mu
  q2p:function(q, o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S = this.sd(o);
    return S.mean(x)-S.mean(y)+ S.qt(q, n1+n2-2)*Math.sqrt(1/n1+1/n2)*S;
  },
  df:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;  
    return n1+n2-2; 
  }
}

var t2paired = { // 成對 T 檢定 T = (X-Y-mu)/(S/sqrt(n)) (雙樣本)
  h:function(o) { return "H0:mu1"+o.op+"mu2" }, 
  sd:function(o) { // S = sd(x-y)
    var x = o.x, n = x.length, y=o.y;
    var S= S.sd(S.sub(x,y));
    return S;
  },
  o2q:function(o) { 
    var x = o.x, n = x.length, y=o.y;
    var S = this.sd(o);
    var t = (S.mean(S.sub(x,y))-o.mu)/(S/Math.sqrt(n));
    return S.pt(t, n-1);
  },
  // mean(x-y)+t*S/sqrt(n)
  q2p:function(q, o) {
    var x = o.x, n = x.length, y=o.y;
    var S = this.sd(o);
    return S.mean(S.sub(x,y))+ S.qt(q, n-1)*S/Math.sqrt(n);
  },
  df:function(o) {
    return o.x.length-1; 
  }
}

var t2varneq = { // σ1≠σ2, Welch's t test (雙樣本) (又稱 Smith-Satterwaite 程序)
// 參考：http://en.wikipedia.org/wiki/Welch's_t_test
  h:function(o) { return "H0:mu1"+o.op+"mu2" }, 
  // T = ((X-Y)-(mu1-mu2))/sqrt(S1^2/n1+S2^2/n2)
  o2q:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1 = S.sd(x), S2=S.sd(y);
    var t = (S.mean(x)-S.mean(y)-o.mu)/Math.sqrt(S1*S1/n1+S2*S2/n2);
    return S.pt(t, this.df(o));
  },
  // t=((X-Y)-mu)/sqrt(S1^2/n1+S2^2/n2), (X-Y)-t*sqrt(S1^2/n1+S2^2/n2) = mu
  q2p:function(q, o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1 = S.sd(x), S2=S.sd(y);
    return S.mean(x)-S.mean(y)+ S.qt(q, this.df(o))*Math.sqrt(S1*S1/n1+S2*S2/n2);
  },
  df:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;  
    var S1 = S.sd(x), S2=S.sd(y);
    var Sn1 = S1*S1/n1, Sn2 = S2*S2/n2, Sn12 = Sn1+Sn2;
    var df = (Sn12*Sn12)/((Sn1*Sn1)/(n1-1)+(Sn2*Sn2)/(n2-1));
    return df;
  }
}

T.ttest = function(o) { 
  var t;
  if (typeof o.y === "undefined") {
    o.name = "ttest(X)";
    o.D = t1;
    t = T.test(o);
    t.mean = S.mean(o.x);
    t.sd   = S.sd(o.x);
  } else {
    var varequal = o.varequal || false;
    var paired = o.paired || false;
    if (varequal) {
      o.name = "ttest(X,Y,mu="+o.mu+",varequal=true) (pooled)";
      o.D = t2vareq;
      t = T.test(o);
    } else if (paired) {
      o.name = "ttest(x,y,mu="+o.mu+",paired=true)";
      o.D = t2paired;
      t = T.test(o);
      t.mean = "mean(x-y)="+S.str(S.mean(S.sub(o.x, o.y)));
      t.sd   = "sd(x-y)="+S.str(S.sd(S.sub(o.x, o.y)));
    } else {
      o.name = "ttest(x,y,mu="+o.mu+",varequal=false), Welch t-test";
      o.D = t2varneq;
      t = T.test(o);
    }
    if (typeof t.mean === "undefined") {
      t.mean = "mean(x)="+S.str(S.mean(o.x))+" mean(y)="+S.str(S.mean(o.y));
      t.sd   = "sd(x)="+S.str(S.sd(o.x))+" sd(y)="+S.str(S.sd(o.y));
    }
  }
  return t;
}

var f2 = { // 檢定 σ1/σ2 = 1? 
  h:function(o) { return "H0:σ1/σ2"+o.op+"1"; }, 
  // F = S1^2/S2^2
  o2q:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1 = S.sd(x), S2=S.sd(y);
    var f = (S1*S1)/(S2*S2);
    var pf = S.pf(f, n1-1, n2-1);
    return pf;
  },
  // 信賴區間公式： S1^2/(S2^2*F(1-α/2), S1^2/(S2^2*F(α/2))
  // 也就是要用 S1^2/(S2^2*f(1-q)) ，參考 R 軟體、應用統計方法 (陳景祥) 389 頁。
  q2p:function(q, o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1 = S.sd(x), S2=S.sd(y);
    var qf = S.qf(1-q, n1-1, n2-1);
    return (S1*S1)/(S2*S2*qf);
  },
  df:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    return [n1-1, n2-1];
  }
}

T.ftest = function(o) { 
  o.name = "ftest(X, Y)";
  o.D = f2;
  var t = T.test(o);
  var rsd = S.sd(o.x)/S.sd(o.y);
  t.ratio = (rsd*rsd);
  return t;
}

// R 軟體沒有此函數，測試請看湯銀才 143 頁
var chisq1 = { // 檢定 σ1 = σ ?
  h:function(o) { return "H0:σ1"+o.op+"σ"; }, 
  // χ(n-1) = (n-1)S^2/σ^2
  o2q:function(o) {
    var x = o.x, n = x.length, S=S.sd(x);
    var v = (n-1)*S*S/(o.sd*o.sd);
    return S.pchisq(v, n-1);
  },
  // 信賴區間公式： (n-1)S^2/χ^2(1-q)，參考 R 語言與統計分析 (湯銀才) 142 頁。
  q2p:function(q, o) {
    var x = o.x, n = x.length, S=S.sd(x);
    return (n-1)*S*S/S.qchisq(1-q, n-1);
  },
  df:function(o) {
    var x = o.x, n = x.length;
    return n-1;
  }
}

T.chisqtest = function(o) { 
  o.name = "chisqtest(X)";
  o.D = chisq1;
  return T.test(o);
}

T.vartest = function(o) {
  if (typeof o.y === "undefined")
    return S.chisqtest(o);
  else
    return S.ftest(o);
}

var z1 = { // 單樣本 Z 檢定
  h:function(o) { return "H0:mu"+o.op+o.mu+" when sd="+o.sd; }, 
  o2q:function(o) {
    var x = o.x, n = x.length;
    var z = (S.mean(x)-o.mu)/(o.sd/Math.sqrt(n)); // z=(X-mu)/(sd/sqrt(n))
    return S.pnorm(z, 0, 1);
  },
  q2p:function(q, o) {
    var x = o.x, n = x.length;
    return S.mean(x) + S.qnorm(q, 0, 1) * S.sd(x) / Math.sqrt(n);
  },
  df:function(o) { return o.x.length; }
}

T.ztest = function(o) { 
  o.name = "ztest(X)";
  o.D = z1;
  return T.test(o);
}

var zprop1 = { // 比例的檢定， n 較大時的近似解 o={ x, n, p } // x 為數值，n 個中出現 x 個 1
  h:function(o) { return "H0:p"+o.op+o.p; }, 
  // Z = (p1-p)/sqrt(p(1-p)/n)
  o2q:function(o) {
    var x=o.x, n=o.n, p1=x/n, p=o.p||p1;
    var z = (p1-p)/Math.sqrt(p*(1-p)/n);
    return S.pnorm(z, 0, 1);
  },
  // 信賴區間公式： p1+z*sqrt(p1*(1-p1)/n)，參考 R 語言與統計分析 (湯銀才) 149 頁。
  q2p:function(q, o) {
    var x=o.x, n=o.n, p1=x/n, p=p1;
    var z = S.qnorm(q, 0, 1);
    var z22n = z*z/(2*n);
    return (p1+z22n+z*Math.sqrt( p*(1-p)/n + z22n/(2*n) ))/(1+2*z22n); // R 的版本，比較複雜的估計公式
//  return p1+z*Math.sqrt(p*(1-p)/n); //  語言與統計分析 (湯銀才) 149 頁的版本。
  },
  df:function(o) { return 1; }
}

var zprop2 = { // 比例的檢定， n 較大時的近似解 o={ x, y, n1, n2 }
  h:function(o) { return "H0:p1-p2"+o.op+o.p; }, 
  // Z = (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2)), p = (n1p1+n2p2)/(n1+n2)，參考 R 語言與統計分析 (湯銀才) 175 頁。
  o2q:function(o) {
    var x=o.x, y=o.y, n1=o.n1, n2=o.n2, p1=x/n1, p2=y/n2, p=(n1*p1+n2*p2)/(n1+n2);
    var z = (p1-p2)/Math.sqrt(p*(1-p)*(1/n1+1/n2));
    return S.pnorm(z, 0, 1);
  },
  // 信賴區間公式： p1-p2+z*sqrt(p*(1-p)*(1/n1+1/n2));
  q2p:function(q, o) {
    var x=o.x, y=o.y, n1=o.n1, n2=o.n2, p1=x/n1, p2=y/n2, p=(n1*p1+n2*p2)/(n1+n2);
    var z = S.qnorm(q, 0, 1);
    return p1-p2+z*Math.sqrt(p*(1-p)*(1/n1+1/n2));
  },
  df:function(o) { return 1; }
}

/* 在 prop.test.R 當中，雙邊檢定的 pvalue 是用 pchisq, 單邊才是用 z ，為何呢？ ( 但是信賴區間則是全部用 z)
 if (alternative == "two.sided")
  PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    else {
  if (k == 1)
      z <- sign(ESTIMATE - p) * sqrt(STATISTIC)
  else
      z <- sign(DELTA) * sqrt(STATISTIC)
  PVAL <- pnorm(z, lower.tail = (alternative == "less"))
    }
*/

T.proptest = function(o) {
  o.p = o.p || 0.5;
  o.name = "proptest("+S.str(o)+")";
  o.correct = o.correct|| false;
  if (o.correct) {
    o.name += ", binomtest";
    o.D += binom1;
  } else {
    if (typeof o.y === "undefined") {
      o.name += ", zprop1";
      o.D = zprop1;
    } else {
      o.p = 0; // p1-p2 = 0
      o.name += ", zprop2";
      o.D = zprop2;
    }
  }
  var t=T.test(o);
  if (typeof o.y === "undefined")
    t.p = o.x/o.n;
  else
    t.p = [o.x/o.n1, o.y/o.n2];
  return t;
}

// 參考： https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/binom.test.R
var binom1 = { // 比例的檢定， n 較大時的近似解 o={ x, n, p } // x 為數值，n 個中出現 x 個 1
  h:function(o) { return "H0:p"+o.op+o.p; }, 
  
  // Z = C(n, k)*p1^k*(1-p1)^(n-k), CDF(z: from 1 to x)
  o2q:function(o) {
    var x=o.x, n=o.n, p = o.p, q;
    var dx = S.dbinom(x, n, p);
    if (o.op === "=") { // 雙尾檢定，去雙尾後 / 2
      var q = 0;
      for (var i=0; i<=n; i++) {
        var di = S.dbinom(i, n, p);
        if (di > dx+1e-5) q += di; // 為何 x 本身不算，如果算應該用 di > dx-1e-5 才對。
      }
      q=1-((1-q)/2); // 因為 test 會 * 2 所進行的減半調整
    } else { // 單尾檢定
      if (Math.abs(x - n*p)<1e-5) // 正確預測， q=1
        q = 1;
      else {
        if (o.op === ">")
          q = S.pbinom(x, n, p); // 去右尾
        else // op=== "<"
          q = S.pbinom(x-1, n, p); // 去右尾
      }
    }
    return q;
  },
/* 注意上述 R 原始碼另一尾的計算方式，是用 < pbinom(最接近 x 者) 的算法，而不是直接 * 2。 問題是我們在 test 中是直接用*2 的方式。
  d <- dbinom(x, n, p)
  ...
  else if (x < m) {
      i <- seq.int(from = ceiling(m), to = n)
      y <- sum(dbinom(i, n, p) <= d * relErr)
      pbinom(x, n, p) 左尾 + pbinom(n - y, n, p, lower.tail = FALSE) 右尾
  } else {
      i <- seq.int(from = 0, to = floor(m))
      y <- sum(dbinom(i, n, p) <= d * relErr)
      pbinom(y - 1, n, p) 左尾 + pbinom(x - 1, n, p, lower.tail = FALSE) 右尾
  }
*/               
  // 信賴區間公式： P(T>c) = Σ (n, i) C(n, i) p1^i (1-p1)^(n-i) for i>c < q
  
  q2p:function(q, o, side) { 
    var x=o.x, n=o.n, p=o.p, op=o.op;
    if (side === "L")
      return S.qbeta(q, x, n - x + 1); // 這裏採用 qbeta 是 R 的作法; 直覺上應該採用 S.qbinom(q, n, p);
    else
      return S.qbeta(q, x + 1, n - x);
  },
  df:function(o) { return 1; }
}

T.binomtest = function(o) {
  o.p = o.p || 0.5;
  o.name = "binomtest("+S.str(o)+")";
  o.D = binom1;
  var t=T.test(o);
  t.p = o.x/o.n;
  t.ci[0]=(o.op === ">")?0:t.ci[0];
  t.ci[1]=(o.op === "<")?1:t.ci[1];
  return t;
}

// anova f-test : array1, array2, array3, ...
T.anovaftest = function() {
  return { 
    h0 : "σ1=σ2=...=σ"+arguments.length, 
    pvalue: J.anovaftest(), 
    score: J.anovafscore(),
  };
}

module.exports = S;
