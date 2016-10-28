var G = require("./statistics");

// 各種 space : https://en.wikipedia.org/wiki/Space_(mathematics)#/media/File:Mathematical_implication_diagram_eng.jpg

G.Space={} // Space = HllbertSpace + BanachSpace + Manifold (流形)

G.HilbertSpace={} // HilbertSpace => InnerProductSpace => LocallyConvaxSpace => VectorSpace (LinearSpace)

G.InnerProductSpace={}

G.LocallyConvaxSpace={}

G.LinearSpace=G.VectorSpace={} // (Algebraic)

G.BanachSpace={} // BanachSpace => NormedVectorSpace => MetricSpace => TopologicalSpace

G.NormedVectorSpace={}

G.MetricSpace={}

G.TopologicalSpace={} // (Analytic)

G.Manifold={}

G.Rn = {} // (R,R,....)

G.steps = function(from, to, step) {
	step=step || 1;
	var a=[];
	for (var t=from; t<=to; t+=step)
		a.push(t);
	return a;
}

G.curve=function(f, from=-10, to=10, step=0.1) {
	var x=G.steps(from, to, step);
	var y=x.map(f);
	return { type:"curve", x:x,	y:y	};
}

G.hist=function(a, from, to, step=1) {
//	console.log("from=%d to=%d step=%d", from, to, step);
	from = from||a.min();
	to = to||a.max();
  var n = Math.ceil((to-from+G.EPSILON)/step);
  var xc = G.steps(from+step/2.0, to, step);
//	console.log("from=%d to=%d step=%d xc=%j", from, to, step, xc);
  var bins = G.newV(n, 0);
  for (var i in a) {
    var slot=Math.floor((a[i]-from)/step);
    if (slot>=0 && slot < n)
      bins[slot]++;
  }
	return { type:'histogram', xc:xc, bins:bins, from:from, to:to, step:step};
}
/*
G.ihist=function(a) {
	console.log("a.min()=%d a.max()=%d", a.min(), a.max());
	return G.hist(a, a.min()-0.5, a.max()+0.5, 1);
}
*/
module.exports=G;