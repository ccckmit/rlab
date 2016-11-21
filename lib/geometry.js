var G = require("./statistics");
var extend = G.extend;

// 各種 space : https://en.wikipedia.org/wiki/Space_(mathematics)#/media/File:Mathematical_implication_diagram_eng.jpg

// Space = HllbertSpace + BanachSpace + Manifold (流形)
G.Space = extend({}, G.Set);

// ============= Topology ====================
// https://en.wikipedia.org/wiki/Topological_space
// p : point
G.TopologicalSpace={ // 拓撲空間 : a Space with neighbor()
	neighbor:function(p1, p2) {}, // TopologicalSpace is a Space with neighbor() function.
  coarser:function(T1, T2) {}, // 更粗糙的 (coarser, weaker or smaller) 的拓樸
	finer:function(T1, T2) { return !courser(T2,T1) }, // 更細緻的 (finer, stronger or larger) 的拓樸
	continuous:function() {}, // if for all x in X and all neighbourhoods N of f(x) there is a neighbourhood M of x such that f(M) is subset of N.
	homeomorphism:function() {}, // 同胚: 同胚是拓撲空間範疇中的同構(存在 f 雙射，連續，且 f-1 也連續) (注意和代數的 Homomorphism (同態) 不同，差一個 e)
}

extend(G.TopologicalSpace, G.Space);

// Kolmogorov classification T0, T1, T2, ....
// https://en.wikipedia.org/wiki/Separation_axiom

G.T0 = G.KolmogorovSpace = {
// every pair of distinct points of X, at least one of them has a neighborhood not containing the other	
}

extend(G.T0, G.TopologicalSpace);

G.T1 = G.FréchetSpace = {}

// 任兩點都可鄰域分離者，為郝斯多夫空間。
G.T2 = G.HausdorffSpace = {}

G.T2Complete = {}

G.T2p5	= G.Urysohn = {} // T2½

G.T3 = G.RegularSpace = {
	// every closed subset C of X and a point p not contained in C admit non-overlapping open neighborhoods
}

G.T3p5	= G.TychonoffSpace = {} // T3½

G.T4 = G.NormalHausdorff = {}

G.T5 = G.CompletelyNormalHausdorff = {}

G.T6 = G.PerfectlyNormalHausdorff = {}

// https://en.wikipedia.org/wiki/Discrete_space
G.DiscreteSpace = {} // 

extend(G.DiscreteSpace, G.TopologicalSpace);

// 1. d(x,y)>=0, 2. d(x,y)=0 iff x=y 3. d(x,y)=d(y,x) 4. d(x,z)<=d(x,y)+d(y,z)
G.Metric = {
	d:function(x,y) {}, 
	test0:function() {
		be(d(x,y)>=0);
		be(d(x,x)==0);
		be(d(x,y)==d(y,x));
	},
	test:function() {
		test0();
		be(d(x,z)<=d(x,y)+d(y,z));
	},
}

G.UltraMetric = {
	test:function() {
		test0();
		be(d(x,z)<=max(d(x,y),d(y,z)));
	},
}

extend(G.UltraMetric, G.Metric);

// https://en.wikipedia.org/wiki/Metric_space
G.MetricSpace = { // distances between all members are defined
	d:function(p1,p2) {},
	test:function() {
		be(d(x,y)>=0);
		be(d(x,x)===0);
		be(d(x,y)===d(y,x));
		be(d(x,z)<=d(x,y)+d(y,z));
	}
}

G.CompleteMetricSpace = { // 空間中的任何柯西序列都收斂在該空間之內。
}

extend(G.CompleteMetricSpace, G.MetricSpace);

// https://en.wikipedia.org/wiki/Complete_metric_space
G.CompleteMetricSpace = {
	// M is called complete (or a Cauchy space) if every Cauchy sequence of points in M has a limit that is also in M or, alternatively, if every Cauchy sequence in M converges in M.
}

// https://en.wikipedia.org/wiki/Uniform_space
G.UniformSpace = { // 帶有一致結構的集合，「x 鄰近於a 勝過y 鄰近於b」之類的概念，在均勻空間中是有意義的。
	
}

G.Rn = {
	
} // (R,R,....)

// 向量空間： VectorSpace = Rn + Linearity
// V:AbelGroup(交換群), F:Field
// 1. F × V → V，(r,v)→ rv
// 2. r(x+y) = rx+ry
// 3. (r+s)x = rx+sx
// 4. (rs)x = r(sx)
// 5. 1x = x
G.VectorSpace = G.LinearSpace = {
	add(x,y) { return x.add(y) },
	mul(a,x) { return a.mul(x) },
	bilinear(b) { // https://en.wikipedia.org/wiki/Bilinear_form
		linear((u)=>b(u,w));
		linear((w)=>b(u,w));
	},
	positiveDefinite(b) { }, // 正定
	symmetric(b) {}, // 對稱
	dualSpace() {}, // https://en.wikipedia.org/wiki/Dual_space
}

extend(G.VectorSpace, G.Rn);

G.FiniteVectorSpace = {} // 有限體向量空間

extend(G.FiniteVectorSpace, G.VectorSpace);

G.NormedVectorSpace={ // 賦範空間
	norm:function(v) { return x.norm() },
}

extend(G.NormedVectorSpace, G.VectorSpace);

G.InnerProductSpace={ // 內積空間
	dot(x,y) { return x.vdot(y) },
}

extend(G.InnerProductSpace, G.NormedVectorSpace);

// 仿射空間：a1 v1 + ... + an vn 中 sum(ai)=1
// https://en.wikipedia.org/wiki/Affine_space
G.AffineSpace = {
	sub(x,y) { }, // 點與點的差是一向量
	add(x,v) { }, // 點加向量得一點，但點與點不可作加法
}

// ============= Euclidean Geometry 歐氏幾何 ====================
G.EuclideanSpace = {
	d(x,y) { return x.sub(y).norm() }, // v= x.sub(y); v.dot(v).sqrt()
	dot(x,y) { return x.vdot(y) },
}

G.LocallyConvaxSpace={}
extend(G.LocallyConvaxSpace, G.VectorSpace);

G.HilbertSpace={} // HilbertSpace => InnerProductSpace => LocallyConvaxSpace => VectorSpace (LinearSpace)
extend(G.HilbertSpace, G.InnerProductSpace);

// BanachSpace => NormedVectorSpace => MetricSpace => TopologicalSpace
// NormedVectorSpace that Cauchy sequence of vectors always converges to a well defined limit
G.BanachSpace={
	
}

extend(G.BanachSpace, G.NormedVectorSpace);

// ============= Elliptic Geometry 橢圓幾何 ====================
G.EllipticGeometry = {} // 橢圓幾何

G.SphericalGeometry = {} // 球面幾何

G.SphericalTrigonometry = {} // 球面三角學

G.Geodesy = {} // 大地測量學

G.GreatCircleDistance = {} // 大球距離

// 圓： x^2+y^2 = 1 =>  (x,y)=(sin t, cos t)
// 雙曲線： x^2-y^2 = 1 => (x,y) = (sinh t, cosh t)
// sinh = (e^x-e^{-x})/2, cosh = (e^x+e^{-x})/2
// https://en.wikipedia.org/wiki/Hyperbolic_function
// sinh x = -i sin ix ; cosh x = cos ix
// sin x  泰勒展開 = x - 1/3! x^3 + 1/5! x^5 ....
// sinh x 泰勒展開 = x + 1/3! x^3 + 1/5! x^5 ....
// int sinh cx dx = 1/c cosh cx + C 
  

// ============= Hyperbolic Geometry 雙曲幾何 ====================
G.HyperbolicGeometry = {}

// ============= Riemannian Geometry ====================
// https://en.wikipedia.org/wiki/Riemannian_geometry
G.RiemannianGeometry = {} // 黎曼幾何

G.RiemannianMetrics = {} // 黎曼度規

G.MetricTensor = {} // 度規張量

G.GaussBonnetTheorem = {} // 高斯博內定理

// ============= Manifold : 流形 ====================
G.Manifold={}

G.C0 = {}

G.Coo = {}

// https://en.wikipedia.org/wiki/Topological_vector_space
G.TopologicalVectorSpace = {}  

// https://en.wikipedia.org/wiki/Locally_convex_topological_vector_space
G.LocallyConvexSpace = {}

extend(G.LocallyConvexSpace, G.TopologicalVectorSpace);
extend(G.LocallyConvexSpace, G.NormedVectorSpace);

// m 維拓撲流形
G.TopologicalManifold={
// M是豪斯多夫空間，x in M 有鄰域 U 同胚於 m 維歐幾里得空間 R^{m}的一個開集
}

G.BanachManifold = {}


// 多邊形：Polygon
// 多面體：Polyhedron : V-E+F = 2
// 多胞形：Polytope
// ============= Spherical Geometry ====================
G.SphericalGeometry = {}


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

// 碎形幾何 : Fractal
// http://andrew-hoyer.com/
// http://andrew-hoyer.com/experiments/fractals/
// http://flam3.com/flame.pdf
// https://en.wikipedia.org/wiki/List_of_fractals_by_Hausdorff_dimension

// http://rembound.com/articles/drawing-mandelbrot-fractals-with-html5-canvas-and-javascript
// https://github.com/rembound/Mandelbrot-Fractal-HTML5/blob/master/mandelbrot-fractal.js

G.MandelbrotSet = {}
/*
For each pixel (Px, Py) on the screen, do:
{
  x0 = scaled x coordinate of pixel (scaled to lie in the Mandelbrot X scale (-2.5, 1))
  y0 = scaled y coordinate of pixel (scaled to lie in the Mandelbrot Y scale (-1, 1))
  x = 0.0
  y = 0.0
  iteration = 0
  max_iteration = 1000
  // Here N=2^8 is chosen as a reasonable bailout radius.
  while ( x*x + y*y < (1 << 16)  AND  iteration < max_iteration ) {
    xtemp = x*x - y*y + x0
    y = 2*x*y + y0
    x = xtemp
    iteration = iteration + 1
  }
  // Used to avoid floating point issues with points inside the set.
  if ( iteration < max_iteration ) {
    // sqrt of inner term removed using log simplification rules.
    log_zn = log( x*x + y*y ) / 2
    nu = log( log_zn / log(2) ) / log(2)
    // Rearranging the potential function.
    // Dividing log_zn by log(2) instead of log(N = 1<<8)
    // because we want the entire palette to range from the
    // center to radius 2, NOT our bailout radius.
    iteration = iteration + 1 - nu
  }
  color1 = palette[floor(iteration)]
  color2 = palette[floor(iteration) + 1]
  // iteration % 1 = fractional part of iteration.
  color = linear_interpolate(color1, color2, iteration % 1)
  plot(Px, Py, color)
}
*/

// 操控繪圖
// http://fabricjs.com/ (讚！)
// https://github.com/kangax/fabric.js

// 3D 動畫
// https://threejs.org/
// http://haptic-data.com/toxiclibsjs/

// 地理資訊
// ArcGIS : https://developers.arcgis.com/javascript/3/

// 動畫
// http://paperjs.org/features/
// https://processing.org/

// 向量
// http://victorjs.org/
// 3d: https://evanw.github.io/lightgl.js/docs/vector.html

module.exports=G;