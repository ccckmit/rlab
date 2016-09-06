var _ = require("lodash");
var R = require("./lib/statistics");
var M = require("./lib/matrix");
var Symbol = require("./lib/symbolic");
var D = require("./lib/calculus");
R.NN = require("./plugin/neural");
R.NN.RBM = require("./plugin/neural/rbm");

R.Matrix = R.M = M;
R._ = _;
R.S = R.Sym = R.Symbol = Symbol;

R.PI = Math.PI;
R.E  = Math.E;

// space 沒有加上機率參數 , 不能指定機率
R.samples = function(space, size, arg) {
	var arg = _.defaults(arg, {replace:true});
	if (arg.replace)
		return _.times(size, ()=>_.sample(space));
	else
		return _.sampleSize(space, size);
}

// Graph

R.G = G = {}

G.curve=function(f, from=-10, to=10, step=0.1) {
	var x=R.steps(from, to, step);
	var y=x.map(f);
	return { type:"curve", x:x,	y:y	};
}

G.hist=function(a, from, to, step=1) {
	from = from||a.min();
	to = to||a.max();
  var n = Math.ceil((to-from+R.EPSILON)/step);
  var xc = R.steps(from+step/2.0, to, step);
  var bins = R.M.newV(n, 0);
  for (var i in a) {
    var slot=Math.floor((a[i]-from)/step);
    if (slot>=0 && slot < n)
      bins[slot]++;
  }
	return { type:'histogram', xc:xc, bins:bins, from:from, to:to, step:step};
}

G.ihist=function(a) {
	return G.hist(a, a.min()-0.5, a.max()+0.5, 1);
}

// Global
debug = function() {
	var arg = _.slice(arguments);
	console.debug.apply(console, arg);
}

print = function() {
	var arg = _.slice(arguments);
	console.log.apply(console, arg);
}

R.debug = debug;
R.print = print;

// ==== Copy functions to R ======
R.copyFunctions(R, D, "differential,integral".split(","));
R.copyFunctions(R, Math, "abs,acos,asin,atan,ceil,cos,exp,floor,log,pow,random,round,sin,sqrt,tan".split(","));

R.copyFunctions(R, M, "solveLP,solveMP,ODE,minimize,complex,spline,linspace".split(","));


// not include : bench, xxxeq, ccsXXX, cgrid, 
R.mixThis(Array.prototype, M, [
"cLU",
"cdelsq",
"clone",
"rows",
"cols",
"row",
"col",
"tr",
// "str",
"not",
"bnot",
"neg",
"abs",
"sin",
"cos",
"tan",
"asin",
"acos",
"atan",
"atan2",
"inv",
"all",
"any",
"same",
"isFinite",
"isNaN",
"sqrt",
"ceil",
"floor",
"round",
"log",
"exp",
"pow",
"mapreduce",
"lshifteq",
"rshifteq",
"add",
"sub",
"mul",
"div",
"mod",
"and",
"or",
"xor",
"band",
"bor",
"bxor",
"eq",
"neq",
"geq",
"leq",
"lt",
"gt",
"complex",
"det",
"norm2",
"norm2Squared",
"norm2inf",
"dot",
"dim",
"eig",
"LU",
"svd",
"sum",
"rowSum",
"colSum",
"rowMean",
"colMean",
"addMV",
"mapM",
"mapMM",
"flatM",
"fillVM",
"fillMM",
"getBlock",
"setBlock",
"getDiag",
"diag",
"parseFloat",
"parseDate",
"parseCSV",
"toCSV",
"strM",
"sumM",
]);

// not include : bench, xxxeq, ccsXXX, cgrid, 
R.mixThis(Array.prototype, G, [
"hist",
"ihist",
]);

R.mixThis(Array.prototype, R, [
// statistics
"max",
"min",
// "sum",
"product",
"mean",
"range",
"median",
"variance",
"deviation",
"sd",
"cov",
"cor",
"normalize",
]);

R.mixThis(Array.prototype, _, [
// lodash
"chunk",
"compact",
// concat:_.concat
"difference",
"differenceBy",
"differenceWith",
"drop",
"dropRight",
"dropRightWhile",
"dropWhile",
// fill:_.fill,
// findIndex:_.findIndex,
"findLastIndex",
"flatten",
"flattenDeep",
"flattenDepth",
"fromPairs",
"head",
// indexOf:_.indexOf,
"initial",
"intersection",
"intersectionBy",
"intersectionWith",
// _.join
"last",
// _.lastIndexOf
"nth",
"pull",
"pullAll",
"pullAllBy",
"pullAllWith",
"pullAt",
"remove",
// _.reverse
// _.slice
"sortedIndex",
"sortedIndexBy",
"sortedIndexOf",
"sortedLastIndex",
"sortedLastIndexBy",
"sortedLastIndexOf",
"sortedUniq",
"sortedUniqBy",
"tail",
"take",
"takeRight",
"takeRightWhile",
"takeWhile",
"union",
"unionBy",
"unionWith",
"uniq",
"uniqBy",
"uniqWith",
"unzip",
"unzipWith",
"without",
// _.xor
// _.xorBy
// _.xorWith
"zip",
"zipObject",
"zipObjectDeep",
"zipWith",
// Collection
"countBy",
// _.each → forEach
// _.eachRight → forEachRight
// every:_.every
// filter:_.filter
// find:_.find
"findLast",
"flatMap",
"flatMapDeep",
"flatMapDepth",
// _.forEach
"forEachRight",
"groupBy",
// includes:_.includes
"invokeMap",
"keyBy",
// _.map
"orderBy",
"partition",
// _.reduce
// reduceRight:_.reduceRight,
"reject",
"sample",
"sampleSize",
"shuffle",
"size",
// some:_.some
"sortBy",
]);

// R.mixThis(Array.prototype,  {str:R.astr}, ['str']);
R.mixThisMap(Array.prototype,  R, {astr:'str',print:'print'});
R.mixThisMap(Number.prototype, R, {nstr:'str',print:'print'});
R.mixThisMap(String.prototype, R, {sstr:'str',print:'print'});
R.mixThisMap(Object.prototype, R, {ostr:'str',print:'print'});
R.mixThisMap(Object.prototype, M, {strM:'strM'});

module.exports = R;


