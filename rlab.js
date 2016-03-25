var R=require("lodash");
var J=require("jStat").jStat;

R.M=require("numeric");

var repeats=function(f, x, n) {
  var a=[];
  for (var i=0; i<n; i++) {
    a.push(f(x));
  }
  return R(a);
}

R.mixin({repeats:repeats});

var calls=function() {
  var args = Array.prototype.slice.call(arguments);
  var n = args[0];
  var f = args[1];
  var params = args.slice(2, args.length);
  var a=[];
  for (var i=0; i<n; i++)
    a.push(f.apply(null, params));
  return R(a);  
}

R.mixin({calls:calls});

var samples=function(space, size, arg) {
  var arg = R.defaults(arg, {replace:true});
  if (arg.replace)
    return R.chain(R.sampleSize(space, size));
  else
    return R.repeats(R.sample, space, size);
}

R.mixin({samples:samples});

R.steps=function(start, end, step) {
  step = step || 1;
  return R.range(start,end+0.00001*step, step); 
}

// ============== format ==================
R.precision = 2;
// R.largeArray = 50;

R.num2str=function(x) {
  if (isNaN(x))
    return "NaN";
  else if (Math.floor(x) === x) // 整數
    return x.toString();
  else if (isFinite(x)) // 有小數點
    return R.ctrim(x.toFixed(R.precision), "0", "right");
  else if (x < 0)
    return "-Infinity";
  else
    return "Infinity";
}

R.ctrim=function(s, set, side) {
  side = side||"both";
  for (var b=0; b<s.length; b++)
    if (set.indexOf(s[b])<0) break;
  for (var e=s.length-1; e>=0; e--)
    if (set.indexOf(s[e])<0) break;
  if (side === "right") b=0;
  if (side === "left") e=s.length-1;
  return s.substring(b, e+1);
}

R.array2str=function(x) {
  var s = "";
  for (var i in x)
    s+= R.str(x[i])+", ";
  return "["+R.ctrim(s, ", ")+"]";
}

R.obj2str=function(x, sp) {
  var s = "";
  for (var k in x)
    s+= k+":"+R.str(x[k])+sp;
  return "{"+R.ctrim(s,sp)+"}";
}

var str = function(x) {
  if (typeof x === "number") return R.num2str(x);
  else if(typeof x === "string") return '\"'+x.toString()+'\"';
  else if (x instanceof Array)  return R.array2str(x);
  else if (typeof x === "object") return R.obj2str(x, ",");
  else return x.toString();
}

R.mixin({str:str}, {chain:false});

// ============== /format ==================

var sd = function(a, flag) { 
  flag = flag || 1;
  return J.stdev(a, flag); 
}

R.mixin({sd:sd}, {chain:false});

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

// ============= 離散分佈 ============
var qf=function(cdf, q, N, p) {
  for (var i=0; i<=N; i++) {
    if (cdf(i, N, p) > q) return i;
  }
  return N;
} 
var rf=function(cdf, n, N, p) {
  var samples = [];
  for (var i=0; i<n; i++) {
    var q = Math.random();
    samples.push(cdf(q, N, p));
  }
  return samples;
} 
// 二項分布
R.dbinom=function(x, N, p) { return J.binomial.pdf(x, N, p); }
R.pbinom=function(k, N, p) { return J.binomial.cdf(k, N, p); }
R.qbinom=function(q, N, p) { return qf(R.pbinom, q, N, p); } 
R.rbinom=function(n, N, p) { return rf(R.qbinom, n, N, p); }
// 負二項分布
R.dnbinom=function(x, N, p) { return jStat.negbin.pdf(x, N, p); }
R.pnbinom=function(k, N, p) { return jStat.negbin.cdf(k, N, p); }
R.qnbinom=function(q, N, p) { return qf(R.pnbinom, q, N, p); } 
R.rnbinom=function(n, N, p) { return rf(R.qnbinom, n, N, p); }
// 超幾何分布
R.dhyper=function(x, N, m, n) { return jStat.hypgeom.pdf(x, N, m, n); }
R.phyper=function(k, N, m, n) { return jStat.hypgeom.cdf(k, N, m, n); }
R.qhyper=function(q, N, m, n) { return qf(R.phyper, q, N, p); } 
R.rhyper=function(n, N, m, k) { return rf(R.qhyper, n, N, p); }

// =============== 檢定 ==============================

// ---------------- 統計與檢定 ----------------------
function opAlt(op) {
  if (op === "=") return "!=";
  if (op === "<") return ">=";
  if (op === ">") return "<=";
  return null;
}

R.test = function(o) { // name, D, x, mu, sd, y, alpha, op
  o = R.defaults(o, {alpha:0.05, op:"="});
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
    report: function() { R.report(this) }
  };
}

var t1 = { // 單樣本 T 檢定 t = (X-mu)/(S/sqrt(n))
  h:function(o) { return "H0:mu"+o.op+o.mu; }, 
  o2q:function(o) {
    var x = o.x, n = x.length;
    var t = (R.mean(x)-o.mu)/(R.sd(x)/Math.sqrt(n)); 
    return R.pt(t, n-1);
  },
  // P(X-mu/(S/sqrt(n))<t) = q ; 信賴區間 P(T<q)
  // P(mu > X-t*S/sqrt(n)) = q ; 這反而成了右尾檢定，所以左尾與右尾確實會反過來
  q2p:function(q, o) {
    var x = o.x, n = x.length;
    return R.mean(x) + R.qt(q, n-1) * R.sd(x) / Math.sqrt(n);
  },
  df:function(o) { return o.x.length-1; }
}

var t2vareq = { // σ1=σ2, 合併 T 檢定 (雙樣本)
  h:function(o) { return "H0:mu1"+o.op+"mu2" }, 
  // S^2 = (n1-1)*S1^2+(n2-1)*S2^2)/(n1-1+n2-1)
  sd:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1= R.sd(x), S2 = R.sd(y);
    var S = Math.sqrt(((n1-1)*S1*S1+(n2-1)*S2*S2)/(n1-1+n2-1)); 
    return S;
  },
  // T = ((X-Y)-(mu1-mu2))/(sqrt(1/n1+1/n2)*S)
  o2q:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S = this.sd(o);
    var t = (R.mean(x)-R.mean(y)-o.mu)/(Math.sqrt(1/n1+1/n2)*S);
    return R.pt(t, n1+n2-2);
  },
  // t=((X-Y)-mu)/(sqrt(1/n1+1/n2)*S), (X-Y)-t*sqrt(1/n1+1/n2)*S = mu
  q2p:function(q, o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S = this.sd(o);
    return R.mean(x)-R.mean(y)+ R.qt(q, n1+n2-2)*Math.sqrt(1/n1+1/n2)*S;
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
    var S= R.sd(R.sub(x,y));
    return S;
  },
  o2q:function(o) { 
    var x = o.x, n = x.length, y=o.y;
    var S = this.sd(o);
    var t = (R.mean(R.sub(x,y))-o.mu)/(S/Math.sqrt(n));
    return R.pt(t, n-1);
  },
  // mean(x-y)+t*S/sqrt(n)
  q2p:function(q, o) {
    var x = o.x, n = x.length, y=o.y;
    var S = this.sd(o);
    return R.mean(R.sub(x,y))+ R.qt(q, n-1)*S/Math.sqrt(n);
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
    var S1 = R.sd(x), S2=R.sd(y);
    var t = (R.mean(x)-R.mean(y)-o.mu)/Math.sqrt(S1*S1/n1+S2*S2/n2);
    return R.pt(t, this.df(o));
  },
  // t=((X-Y)-mu)/sqrt(S1^2/n1+S2^2/n2), (X-Y)-t*sqrt(S1^2/n1+S2^2/n2) = mu
  q2p:function(q, o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1 = R.sd(x), S2=R.sd(y);
    return R.mean(x)-R.mean(y)+ R.qt(q, this.df(o))*Math.sqrt(S1*S1/n1+S2*S2/n2);
  },
  df:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;  
    var S1 = R.sd(x), S2=R.sd(y);
    var Sn1 = S1*S1/n1, Sn2 = S2*S2/n2, Sn12 = Sn1+Sn2;
    var df = (Sn12*Sn12)/((Sn1*Sn1)/(n1-1)+(Sn2*Sn2)/(n2-1));
    return df;
  }
}

R.ttest = function(o) { 
  var t;
  if (typeof o.y === "undefined") {
    o.name = "ttest(X)";
    o.D = t1;
    t = R.test(o);
    t.mean = R.mean(o.x);
    t.sd   = R.sd(o.x);
  } else {
    var varequal = R.opt(o, "varequal", false);
    var paired = R.opt(o, "paired", false);
    if (varequal) {
      o.name = "ttest(X,Y,mu="+o.mu+",varequal=true) (pooled)";
      o.D = t2vareq;
      t = R.test(o);
    } else if (paired) {
      o.name = "ttest(x,y,mu="+o.mu+",paired=true)";
      o.D = t2paired;
      t = R.test(o);
      t.mean = "mean(x-y)="+R.str(R.mean(R.sub(o.x, o.y)));
      t.sd   = "sd(x-y)="+R.str(R.sd(R.sub(o.x, o.y)));
    } else {
      o.name = "ttest(x,y,mu="+o.mu+",varequal=false), Welch t-test";
      o.D = t2varneq;
      t = R.test(o);
    }
    if (typeof t.mean === "undefined") {
      t.mean = "mean(x)="+R.str(R.mean(o.x))+" mean(y)="+R.str(R.mean(o.y));
      t.sd   = "sd(x)="+R.str(R.sd(o.x))+" sd(y)="+R.str(R.sd(o.y));
    }
  }
  return t;
}

var f2 = { // 檢定 σ1/σ2 = 1? 
  h:function(o) { return "H0:σ1/σ2"+o.op+"1"; }, 
  // F = S1^2/S2^2
  o2q:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1 = R.sd(x), S2=R.sd(y);
    var f = (S1*S1)/(S2*S2);
    var pf = R.pf(f, n1-1, n2-1);
    return pf;
  },
  // 信賴區間公式： S1^2/(S2^2*F(1-α/2), S1^2/(S2^2*F(α/2))
  // 也就是要用 S1^2/(S2^2*f(1-q)) ，參考 R 軟體、應用統計方法 (陳景祥) 389 頁。
  q2p:function(q, o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    var S1 = R.sd(x), S2=R.sd(y);
    var qf = R.qf(1-q, n1-1, n2-1);
    return (S1*S1)/(S2*S2*qf);
  },
  df:function(o) {
    var x = o.x, n1 = x.length, y=o.y, n2=y.length;
    return [n1-1, n2-1];
  }
}

R.ftest = function(o) { 
  o.name = "ftest(X, Y)";
  o.D = f2;
  var t = R.test(o);
  var rsd = R.sd(o.x)/R.sd(o.y);
  t.ratio = (rsd*rsd);
  return t;
}

// R 軟體沒有此函數，測試請看湯銀才 143 頁
var chisq1 = { // 檢定 σ1 = σ ?
  h:function(o) { return "H0:σ1"+o.op+"σ"; }, 
  // χ(n-1) = (n-1)S^2/σ^2
  o2q:function(o) {
    var x = o.x, n = x.length, S=R.sd(x);
    var v = (n-1)*S*S/(o.sd*o.sd);
    return R.pchisq(v, n-1);
  },
  // 信賴區間公式： (n-1)S^2/χ^2(1-q)，參考 R 語言與統計分析 (湯銀才) 142 頁。
  q2p:function(q, o) {
    var x = o.x, n = x.length, S=R.sd(x);
    return (n-1)*S*S/R.qchisq(1-q, n-1);
  },
  df:function(o) {
    var x = o.x, n = x.length;
    return n-1;
  }
}

R.chisqtest = function(o) { 
  o.name = "chisqtest(X)";
  o.D = chisq1;
  return R.test(o);
}

R.vartest = function(o) {
  if (typeof o.y === "undefined")
    return R.chisqtest(o);
  else
    return R.ftest(o);
}

var z1 = { // 單樣本 Z 檢定
  h:function(o) { return "H0:mu"+o.op+o.mu+" when sd="+o.sd; }, 
  o2q:function(o) {
    var x = o.x, n = x.length;
    var z = (R.mean(x)-o.mu)/(o.sd/Math.sqrt(n)); // z=(X-mu)/(sd/sqrt(n))
    return R.pnorm(z, 0, 1);
  },
  q2p:function(q, o) {
    var x = o.x, n = x.length;
    return R.mean(x) + R.qnorm(q, 0, 1) * R.sd(x) / Math.sqrt(n);
  },
  df:function(o) { return o.x.length; }
}

R.ztest = function(o) { 
  o.name = "ztest(X)";
  o.D = z1;
  return R.test(o);
}

var zprop1 = { // 比例的檢定， n 較大時的近似解 o={ x, n, p } // x 為數值，n 個中出現 x 個 1
  h:function(o) { return "H0:p"+o.op+o.p; }, 
  // Z = (p1-p)/sqrt(p(1-p)/n)
  o2q:function(o) {
    var x=o.x, n=o.n, p1=x/n, p=R.def(o.p, p1);
    var z = (p1-p)/Math.sqrt(p*(1-p)/n);
    return R.pnorm(z, 0, 1);
  },
  // 信賴區間公式： p1+z*sqrt(p1*(1-p1)/n)，參考 R 語言與統計分析 (湯銀才) 149 頁。
  q2p:function(q, o) {
    var x=o.x, n=o.n, p1=x/n, p=p1;
    var z = R.qnorm(q, 0, 1);
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
    return R.pnorm(z, 0, 1);
  },
  // 信賴區間公式： p1-p2+z*sqrt(p*(1-p)*(1/n1+1/n2));
  q2p:function(q, o) {
    var x=o.x, y=o.y, n1=o.n1, n2=o.n2, p1=x/n1, p2=y/n2, p=(n1*p1+n2*p2)/(n1+n2);
    var z = R.qnorm(q, 0, 1);
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

R.proptest = function(o) {
  o.p = R.opt(o, "p", 0.5);
  o.name = "proptest("+R.str(o)+")";
  o.correct = R.opt(o.correct, false);
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
  var t=R.test(o);
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
    var dx = R.dbinom(x, n, p);
    if (o.op === "=") { // 雙尾檢定，去雙尾後 / 2
      var q = 0;
      for (var i=0; i<=n; i++) {
        var di = R.dbinom(i, n, p);
        if (di > dx+1e-5) q += di; // 為何 x 本身不算，如果算應該用 di > dx-1e-5 才對。
      }
      q=1-((1-q)/2); // 因為 test 會 * 2 所進行的減半調整
    } else { // 單尾檢定
      if (Math.abs(x - n*p)<1e-5) // 正確預測， q=1
        q = 1;
      else {
        if (o.op === ">")
          q = R.pbinom(x, n, p); // 去右尾
        else // op=== "<"
          q = R.pbinom(x-1, n, p); // 去右尾
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
      return qbeta(q, x, n - x + 1); // 這裏採用 qbeta 是 R 的作法; 直覺上應該採用 R.qbinom(q, n, p);
    else
      return qbeta(q, x + 1, n - x);
  },
  df:function(o) { return 1; }
}

R.binomtest = function(o) {
  o.p = R.opt(o, "p", 0.5);
  o.name = "binomtest("+R.str(o)+")";
  o.D = binom1;
  var t=R.test(o);
  t.p = o.x/o.n;
  t.ci[0]=(o.op === ">")?0:t.ci[0];
  t.ci[1]=(o.op === "<")?1:t.ci[1];
  return t;
}

R.report = function(o) {
  console.log("=========== report ==========");
  for (var k in o) {
    if (typeof o[k] !== "function")
      console.log(k+"\t: "+str(o[k]));
  }
}

// anova f-test : array1, array2, array3, ...
R.anovaftest = function() {
  return { 
    h0 : "σ1=σ2=...=σ"+arguments.length, 
    pvalue: J.anovaftest(), 
    score: J.anovafscore(),
  };
}

module.exports = R;
