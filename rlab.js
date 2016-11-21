var R    = require("./lib/math");
var _ = R._ = require("lodash");
R.Symbol = require("./lib/symbolic");
R.NN     = require("./plugin/neural");
R.NN.RBM = require("./plugin/neural/rbm");

// space 沒有加上機率參數 , 不能指定機率
R.samples = function(space, size, arg) {
	var arg = _.defaults(arg, {replace:true});
	if (arg.replace)
		return _.times(size, ()=>_.sample(space));
	else
		return _.sampleSize(space, size);
}

// Global
R.debug = debug = function() {
	var arg = _.slice(arguments);
	console.debug.apply(console, arg);
}

R.print = print = function() {
	var arg = _.slice(arguments);
	console.log.apply(console, arg);
}

p = R.parse;
be = R.be;

R.mixThisMap(Array.prototype, R, {
lu:"lu",
luSolve:"luSolve",
svd:"svd",
// "cdelsq",
// "clone",
rows:"rows",
cols:"cols",
row:"row",
col:"col",
tr:"tr",
inv:"inv",
// "all",
// "any",
// "same",
// "isFinite",
// "isNaN",
// "mapreduce",
// "complex",
det:"det",
// "norm2",
// "norm2Squared",
// "norm2inf",
dot:"dot",
// "dim",
eig:"eig",
// "sum",
rowSum:"rowSum",
colSum:"colSum",
rowMean:"rowMean",
colMean:"colMean",
addMV:"addMV",
mapM:"mapM",
mapMM:"mapMM",
flatM:"flatM",
fillVM:"fillVM",
fillMM:"fillMM",
getBlock:"getBlock",
setBlock:"setBlock",
getDiag:"getDiag",
diag:"diag",
// "parseFloat",
// "parseDate",
// "parseCSV",
// "toCSV",
mstr:"mstr",
// "sumM",
str:'astr',
print:'print',
});

R.mixThis(Array.prototype, R, [
"max",
"min",
"sum",
"product",
"norm",
"mean",
"range",
"median",
"variance",
"deviation",
"sd",
"cov",
"cor",
"normalize",
"curve",
"hist",
"ihist",
"eval",
// "map",
// +-*/%
"add",
"sub",
"mul",
"div",
"mod",
"neg",
// "inv", 和矩陣相衝
// logical
"and",
"or",
"xor",
"not",
// bits operation
"bnot",
"band",
"bor",
"bxor",
// function
"power",
// "dot", 和矩陣相衝
"sqrt",
"log",
"exp",
"abs",
"sin",
"cos",
"tan",
"asin",
"acos",
"atan",
"ceil",
"floor",
"round",
]);

R.mixThisMap(Number.prototype, R, {
	str:'nstr',
	print:'print',
});

R.mixThis(Number.prototype, R, [
	'map',
	'add',
	'sub',
	'mul',
	'div',
	'mod',
	'power',
	'neg',
	'inv',
	'sqrt',
	'log',
	'exp',
	'abs',
	'sin',
	'cos',
	'tan',
	'asin',
	'acos',
	'atan',
	'ceil',
	'floor',
	'round',
]);
R.mixThisMap(Function.prototype, R, {
	add:'fadd',
	sub:'fsub',
	mul:'fmul',
	div:'fdiv',
	compose:'fcompose',
	eval:'feval',
	diff:'fdiff',
	integral:'fintegral',
});
R.mixThisMap(String.prototype, R, {str:'sstr',print:'print'});
R.mixThisMap(Object.prototype, R, {str:'ostr',print:'print'});
R.mixThisMap(Object.prototype, R, {
	proto:'proto',
	eq:'eq',
	neq:'neq',
	geq:'geq',
	leq:'leq',
	gt:'gt',
	lt:'lt',
});

R.mixThisMap(Array.prototype, R._, {
// lodash
_chunk:'chunk',
_compact:'compact',
_concat:'concat',
_difference:'difference',
_differenceBy:'differenceBy',
_differenceWith:'differenceWith',
_drop:'drop',
_dropRight:'dropRight',
_dropRightWhile:'dropRightWhile',
_dropWhile:'dropWhile',
_fill:'fill',
_findIndex:'findIndex',
_findLastIndex:'findLastIndex',
_flatten:'flatten',
_flattenDeep:'flattenDeep',
_flattenDepth:'flattenDepth',
_fromPairs:'flattenPairs',
_head:'head',
_indexOf:'indexOf',
_initial:'initial',
_intersection:'intersection',
_intersectionBy:'intersectonBy',
_intersectionWith:'intersectionWith',
_join:'join',
_last:'last',
_lastIndexOf:'lastIndexOf',
_nth:'nth',
_pull:'pull',
_pullAll:'pullAll',
_pullAllBy:'pullAllBy',
_pullAllWith:'pullAllWith',
_pullAt:'pullAt',
_remove:'remove',
_reverse:'reverse',
_slice:'slice',
_sortedIndex:'sortedIndex',
_sortedIndexBy:'sortedIndexBy',
_sortedIndexOf:'sortedIndexOf',
_sortedLastIndex:'sortedLastIndex',
_sortedLastIndexBy:'sortedLastIndexBy',
_sortedLastIndexOf:'sortedLastIndexOf',
_sortedUniq:'sortedUniq',
_sortedUniqBy:'sortedUniqBy',
_tail:'tail',
_take:'take',
_takeRight:'takeRight',
_takeRightWhile:'takeRightWhile',
_takeWhile:'takeWhile',
_union:'union',
_unionBy:'unionBy',
_unionWith:'unionWith',
_uniq:'uniq',
_uniqBy:'uniqBy',
_uniqWith:'uniqWith',
_unzip:'unzip',
_unzipWith:'unzipWith',
_without:'without',
_xor:'xor',
_xorBy:'xorBy',
_xorWith:'xorWith',
_zip:'zip',
_zipObject:'zipObject',
_zipObjectDeep:'zipObjectDeep',
_zipWith:'zipWith',
// Collection
_countBy:'countBy',
// each→ forEach
// _eachRight → forEachRight
_every:'every',
_filter:'filter',
_find:'find',
_findLast:'findLast',
_flatMap:'flatMap',
_flatMapDeep:'flatMapDeep',
_flatMapDepth:'flatMapDepth',
_forEach:'forEach',
_forEachRight:'forEachRight',
_groupBy:'groupBy',
_includes:'includes',
_invokeMap:'invokeMap',
_keyBy:'keyBy',
_map:'map',
_orderBy:'orderBy',
_partition:'partition',
_reduce:'reduce',
_reduceRight:'reduceRight',
_reject:'reject',
_sample:'sample',
_sampleSize:'sampleSize',
_shuffle:'shuffle',
_size:'size',
_some:'some',
_sortBy:'sortBy',
});

module.exports = R;

// R.mixThis(Array.prototype,  {str:R.astr}, ['str']);
// R.mixThisMap(Array.prototype, R, {str:'astr',print:'print'});
// R.mixThisMap(Object.prototype, R, {strM:'strM'});
// ==== Copy functions to R ======
// R.copyFunctions(R, D, "differential,integral".split(","));
// R.copyFunctions(R, Math, "abs,acos,asin,atan,ceil,cos,exp,floor,log,pow,random,round,sin,sqrt,tan".split(","));

// R.copyFunctions(R, M, "solveLP,solveMP,ODE,minimize,complex,spline,linspace".split(","));

// not include : bench, xxxeq, ccsXXX, cgrid, 
// Complex = R.Complex;
// Ratio = R.Ratio;
