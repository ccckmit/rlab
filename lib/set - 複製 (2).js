var R = require("./rule");
var I = require("./integer");
var S = {}
// ========== Set =================
var create=S.create=function(has) {
	return { has:has }
}

S.union=function(x,y) {
  return create(function(a) { return x.has(a) || y.has(a) })
}

S.intersection=function(x,y) { 
	return create(function(a) { return x.has(a) && y.has(a) })
}

S.difference=function(x,y) { 
	return create(function(a) { return x.has(a) && !y.has(a) })
}

S.symmetricDifference=function(x,y) { 
	return S.difference(S.union(x,y), S.intersection(x,y));
}

S.cartesianProduct=function(x,y) { 
	return create(function(a) { return x.has(a[0]) && y.has(a[1]) })
}

// 因為羅素悖論 A={x|x∉x}，對於無限集合，powerSet 很難定義，除非改採用公理化集合論，但那個很麻煩。
S.powerSet=function(set) {
	return create(function(subset) { 
		for (var e of subset) {
			if (!set.has(e)) return false;
		}
		return true;
	})
}

S.subSet=function(subset, set) { // 此方法僅限用於有限集合。
	for (var a of subset) {
		if (!set.has(a)) return false;
	}
	return true;
}

// 所有東西的集合
S.All=create(R.yes)
// 空集合
S.Empty=create(R.no)
// 浮點數集合
S.Float=create(R.isFloat)
// 整數集合
S.Z=S.Integer=create(R.isInteger);
// 自然數集合
S.N0=S.N=create(function(a) { R.isInteger(a)&&a>=0 })
S.N1=create(function(a) { R.isInteger(a)&&a>=1 })
// 偶數集合
S.Even=create(function(a) { return (R.isInteger(a) && a%2===0) } )
// 奇數集合
S.Odd=create(function(a) { return (R.isInteger(a) && a%2===1) })
// 質數集合
S.Prime=create(function(a) { return I.isPrime(a) })
// 函數集合
S.Function=create(function(a) { return typeof a==='function' })
// 字串集合
S.String=create(function(a) { return typeof a==='string' })
// 物件集合
S.Object=create(function(a) { return typeof a==='object' })
// 陣列集合
S.Array=create(function(a) { return a.constructor===Array })
// 有限元素的集合
S.Set=Set; // Finite=Default JavaScript Set Object

// =============== measure =======================
// ref : https://en.wikipedia.org/wiki/Measure_(mathematics)
// ref : https://en.wikipedia.org/wiki/Sigma-algebra
S.monotonicity=function(s1,s2,measure) {
	if (S.subset(s1,s2)) leq(measure(s1), measure(s2));
}

S.count=function(set) { }

module.exports=S;