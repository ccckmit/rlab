var _ = require("lodash");
var R = require("./lib/R");
var M = require("./lib/matrix");
var Math = require("./lib/math");
var Symbol = require("./lib/symbolic");
var T = require("./lib/test");
var D = require("./lib/calculus");
R.NN = require("./plugin/neural");
R.NN.RBM = require("./plugin/neural/rbm");

R.Matrix = R.M = M;
R.Differential = R.D = D;
R._ = _;
R.S = R.Sym = R.Symbol = Symbol;

R.mix(R, Math);
R.mix(R, T);

// Global
debug = function() {
	var arg = _.slice(arguments);
	console.debug.apply(console, arg);
}

log = function() {
	var arg = _.slice(arguments);
	console.log.apply(console, arg);
}

// space 沒有加上機率參數 , 不能指定機率
R.samples = function(space, size, arg) {
	var arg = _.defaults(arg, {replace:true});
	if (arg.replace)
		return R.calls(size, function() { return _.sample(space); });
	else
		return _.sampleSize(space, size);
}

R.ODE=M.dopri; // dopri(x0,x1,y0,f,tol,maxit,event) #Ordinary Diff Eq
R.solveLP =M.solveLP; // solveLP(c,A,b,Aeq,beq,tol,maxit) #Linear programming
R.solveQP =M.solveQP; // solveQP(Dmat, dvec, Amat, bvec, meq, factorized); // Quadratic Programming
R.minimize=M.uncmin; // uncmin(f,x0,tol,gradient,maxit,callback,options); // Unconstrained optimization
R.sparse=M.sparse; // Matrix => Sparse
R.sparse2full=M.sparse2full; // Sparse => Matrix
R.complex=M.t;
R.spline=M.spline;
R.linspace=M.linspace;

// not include : bench, xxxeq, ccsXXX, cgrid, 
R.mixThis(Array.prototype, {
cLU:M.cLU,
cdelsq:M.cdelsq, // Laplacian
clone:M.clone,
// matrix
rows:M.rows,
cols:M.cols,
row:M.row,
col:M.col,
tr:M.tr,
strM:M.str,
not:M.not,
bnot:M.bnot,
neg:M.neg,
abs:M.abs,
sin:M.sin,
cos:M.cos,
tan:M.tan,
asin:M.asin,
acos:M.acos,
atan:M.atan,
//https://zh.wikipedia.org/wiki/Atan2
atan2:M.atan2,
inv:M.inv,
all:M.all,
any:M.any,
same:M.same,
isFinite:M.isFinite,
isNaN:M.isNaN,
sqrt:M.sqrt,
ceil:M.ceil,
floor:M.floor,
round:M.round,
log:M.log,
exp:M.exp,
pow:M.pow,
mapreduce:M.mapreduce,
lshifteq:M.lshifteq,
rshifteq:M.rshifteq,
add:M.add,
sub:M.sub,
mul:M.mul,
div:M.div,
mod:M.mod,
and:M.and,
or:M.or,
xor:M.xor,
band:M.band,
bor:M.bor,
bxor:M.bxor,
eq:M.eq,
neq:M.neq,
geq:M.geq,
leq:M.leq,
lt:M.lt,
gt:M.gt,
t:M.t,
det:M.det,
norm2:M.norm2,
norm2Squared:M.norm2Squared,
norm2inf:M.norm2inf,
dot:M.dot,
det:M.det,
dim:M.dim,
eig:M.eig,
LU:M.LU,
svd:M.svd,
sumM:M.sum,
rowSum:M.rowSum,
colSum:M.colSum,
rowMean:M.rowMean,
colMean:M.colMean,
addMV:M.addMV,
mapM:M.mapM, 
mapMM:M.mapMM,
flatM:M.flatM,
getBlock:M.getBlock,
setBlock:M.setBlock,
getDiag:M.getDiag,
diag:M.diag,
parseFloat:M.parseFloat,
parseDate:M.parseDate,
parseCSV:M.parseCSV,
toCSV:M.toCSV,

// statistics
max:R.max,
min:R.min,
sum:R.sum,
product:R.product,
mean:R.mean,
range:R.range,
median:R.median,
variance:R.variance,
deviation:R.deviation,
sd:R.sd,
cov:R.cov,
cor:R.cor,
normalize:R.normalize,
hist:R.hist,

// lodash
chunk:_.chunk,
compact:_.compact,
// concat:_.concat
difference:_.difference,
differenceBy:_.differenceBy,
differenceWith:_.differenceWith,
drop:_.drop,
dropRight:_.dropRight,
dropRightWhile:_.dropRightWhile,
dropWhile:_.dropWhile,
// fill:_.fill,
// findIndex:_.findIndex,
findLastIndex:_.findLastIndex,
flatten:_.flatten,
flattenDeep:_.flattenDeep,
flattenDepth:_.flattenDepth,
fromPairs:_.fromPairs,
head:_.head,
// indexOf:_.indexOf,
initial:_.initial,
intersection:_.intersection,
intersectionBy:_.intersectionBy,
intersectionWith:_.intersectionWith,
// _.join
last:_.last,
// _.lastIndexOf
nth:_.nth,
pull:_.pull,
pullAll:_.pullAll,
pullAllBy:_.pullAllBy,
pullAllWith:_.pullAllWith,
pullAt:_.pullAt,
remove:_.remove,
// _.reverse
// _.slice
sortedIndex:_.sortedIndex,
sortedIndexBy:_.sortedIndexBy,
sortedIndexOf:_.sortedIndexOf,
sortedLastIndex:_.sortedLastIndex,
sortedLastIndexBy:_.sortedLastIndexBy,
sortedLastIndexOf:_.sortedLastIndexOf,
sortedUniq:_.sortedUniq,
sortedUniqBy:_.sortedUniqBy,
tail:_.tail,
take:_.take,
takeRight:_.takeRight,
takeRightWhile:_.takeRightWhile,
takeWhile:_.takeWhile,
union:_.union,
unionBy:_.unionBy,
unionWith:_.unionWith,
uniq:_.uniq,
uniqBy:_.uniqBy,
uniqWith:_.uniqWith,
unzip:_.unzip,
unzipWith:_.unzipWith,
without:_.without,
// _.xor
// _.xorBy
// _.xorWith
zip:_.zip,
zipObject:_.zipObject,
zipObjectDeep:_.zipObjectDeep,
zipWith:_.zipWith,
// Collection
countBy:_.countBy,
// _.each → forEach
// _.eachRight → forEachRight
// every:_.every
// filter:_.filter
// find:_.find
findLast:_.findLast,
flatMap:_.flatMap,
flatMapDeep:_.flatMapDeep,
flatMapDepth:_.flatMapDepth,
// _.forEach
forEachRight:_.forEachRight,
groupBy:_.groupBy,
// includes:_.includes
invokeMap:_.invokeMap,
keyBy:_.keyBy,
// _.map
orderBy:_.orderBy,
partition:_.partition,
// _.reduce
// reduceRight:_.reduceRight,
reject:_.reject,
sample:_.sample,
sampleSize:_.sampleSize,
shuffle:_.shuffle,
size:_.size,
// some:_.some
sortBy:_.sortBy,
});

B = R;
B.precision=4;

B.mixThis(Number.prototype, { 
str:function(n, precision=B.precision) {
	return n.toFixed(precision);
},
});

B.mixThis(Array.prototype, { str:function(a, precision=B.precision) {
	var s="";
	for (var i in a) {
		s+=a[i].str(precision)+",";
	}
	return "["+B.ctrim(s, ',')+"]";
}});

B.mixThis(String.prototype, { 
str:function(s) {
	return s.toString();
},
lpad:B.lpad,
});

B.mixThis(Object.prototype, { str:function(o, precision) {
  var s = "";
  for (var k in o)
    s+= k+":"+B.str(o[k], precision)+",";
  return "{"+B.ctrim(s,",")+"}";
}});

B.mixThis(Object.prototype, { strM:M.str });

B.str=function(o) { return o.str(); }

module.exports = R;