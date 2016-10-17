var S = {}
var R = S.Rule = require("./rule");
var I = require("./integer");
// ========== Set =================
var create=S.create=function(has) {
  return { has:has }
}

S.setUnion=function(set) {
  return create(function(e) {
    for (var s of set) {
      if (s.has(e)) return true;
    }
    return false;
  })
}

S.setIntersection=function(set) {
  return create(function(e) {
    for (var s of set) {
      if (!s.has(e)) return false;
    }
    return true;
  })
}

S.union=function(x,y) { return S.setUnion([x,y]) }

S.intersection=function(x,y) { return S.setIntersection([x,y]) }

S.difference=function(x,y) { 
  return create(function(e) { return x.has(e) && !y.has(e) })
}

S.symmetricDifference=function(x,y) { 
  return S.difference(S.union(x,y), S.intersection(x,y));
}

S.cartesianProduct=function(x,y) { 
  return create(function(e) { return x.has(e[0]) && y.has(e[1]) })
}

// 因為羅素悖論 e={x|x∉x}，對於無限集合，powerSet 很難定義，除非改採用公理化集合論，但那個很麻煩。
S.powerSet=function(set) {
  return create(function(subset) { 
    for (var e of subset) {
      if (!set.has(e)) return false;
    }
    return true;
  })
}

S.subSet=function(subset, set) { // 此方法僅限用於有限集合。
  for (var e of subset) {
    if (!set.has(e)) return false;
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
S.N0=S.N=create(function(e) { return R.isInteger(e)&&e>=0 })
S.N1=create(function(e) { return R.isInteger(e)&&e>=1 })
// 偶數集合
S.Even=create(function(e) { return (R.isInteger(e) && e%2===0) } )
// 奇數集合
S.Odd=create(function(e) { return (R.isInteger(e) && e%2===1) })
// 質數集合
S.Prime=create(function(e) { return I.isPrime(e) })
// 函數集合
S.Function=create(function(e) { return typeof e==='function' })
// 字串集合
S.String=create(function(e) { return typeof e==='string' })
// 物件集合
S.Object=create(function(e) { return typeof e==='object' })
// 陣列集合
S.Array=create(function(e) { return e.constructor===Array })
// 有限元素的集合
S.Set=Set; // Finite=Default JavaScript Set Object
S.RussellSet=create(function(e) { // 羅素集合, 悖論 A={ x | !x.has(x) }
  return !e.has(e); // 問題： Is A in A ?
	// 當您呼叫 A = S.RussellSet; A.has(A) 時會當機。
});

// =============== measure =======================
// ref : https://en.wikipedia.org/wiki/Measure_(mathematics)
// ref : https://en.wikipedia.org/wiki/Sigma-algebra
S.monotonicity=function(s1,s2,measure) {
  if (S.subset(s1,s2)) return measure(s1)<measure(s2);
}

S.count=function(set) { }

// https://en.wikipedia.org/wiki/Partially_ordered_set
// S.PartiallyOrderedSet=function()

module.exports=S;