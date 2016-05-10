// jStat API : https://jstat.github.io/all.html
var B = require("./base");
var P = require("./probability");
var S = require("./statistics");

var R = B;
B.mix(R, P);
B.mix(R, S);

module.exports = R;

// var T = require("./lib/test");
// var M = require("./lib/matrix");
// var _ = require("lodash");
/*
console.log("S=", S);

R.M = M;
R.NN = require("./lib/neural");
R.NN.RBM = require("./lib/neural/rbm");

M.precision = R.precision;
*/
/*
console.log("mix B");
B.mix(R, B);
console.log("mix M");
R.mix(R, M);
console.log("mix S");
R.mix(R, S);
console.log("mix P");
R.mix(R, P);
console.log("mix T");
R.mix(R, T);
console.log("mix _");
R.mix(R, _);
// R.mix(R, NN);
R.NN = NN;
R.NN.RBM = RBM;
*/
/*
console.log("mix ...");
B.mix(R, {
samples:function(space, size, arg) { 
		var arg = _.defaults(arg, {replace:true});
		if (arg.replace)
			return R.calls(size, function() { return _.sample(space); });
		else
			return _.sampleSize(space, size);
},
});
// console.log("R.randomM=", R.randomM);

console.log("amix");
B.amix({
// matrix
rows:M.rows,
cols:M.cols,
row:M.row,
col:M.col,
tr:M.tr,
strM:M.strM,
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
sumM:M.sumM,
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
max:S.max,
min:S.min,
sum:S.sum,
product:S.product,
mean:S.mean,
range:S.range,
unique:S.unique,
median:S.median,
variance:S.variance,
deviation:S.deviation,
sd:S.sd,
cov:S.cov,
cor:S.cor,
normalize:S.normalize,
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
fill:_.fill,
findIndex:_.findIndex,
findLastIndex:_.findLastIndex,
flatten:_.flatten,
flattenDeep:_.flattenDeep,
flattenDepth:_.flattenDepth,
fromPairs:_.fromPairs,
head:_.head,
indexOf:_.indexOf,
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
reduceRight:_.reduceRight,
reject:_.reject,
sample:_.sample,
sampleSize:_.sampleSize,
shuffle:_.shuffle,
size:_.size,
// some:_.some
sortBy:_.sortBy,
});
*/
