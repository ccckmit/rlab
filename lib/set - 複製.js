var I = require("./integer");

class _Set {
	constructor(f) { this.f = f }
	has(a) { return this.f(a) }
}

// ========== Set =================
var S={	Set:_Set, Finite:Set } // Finite=Default JavaScript Set Object

S.union=function(x,y) { 
  return new S.Set(function(a) { return x.has(a) || y.has(a) })
}

S.intersection=function(x,y) { 
	return new S.Set(function(a) { return x.has(a) && y.has(a) })
}

S.difference=function(x,y) { 
	return new S.Set(function(a) { return x.has(a) && !y.has(a) })
}

S.symmetricDifference=function(x,y) { 
	return S.difference(S.union(x,y), S.intersection(x,y));
}

S.cartesianProduct=function(x,y) { 
	return new S.Set(function(a) { return x.has(a[0]) && y.has(a[1]) })
}

// 因為羅素悖論 A={x|x∉x}，對於無限集合，powerSet 很難定義，除非改採用公理化集合論，但那個很麻煩。
S.powerSet=function(set) {
	return new S.Set(function(subset) { 
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

S.Float=new S.Set(function(a) { return typeof a==='number' })

S.Integer=new S.Set(function(a) { return (typeof a==='number' && a%1===0) })

S.Even=new S.Set(function(a) { return (S.Integer.has(a) && a%2===0) } )

S.Odd=new S.Set(function(a) { return (S.Integer.has(a) && a%2===1) })

S.Function=new S.Set(function(a) { return typeof a==='function' })

S.String=new S.Set(function(a) { return typeof a==='string' })

S.Object=new S.Set(function(a) { return typeof a==='object' })

S.Array=new S.Set(function(a) { return a.constructor===Array })

S.Prime=new S.Set(function(a) { return I.isPrime(a) })

module.exports=S;