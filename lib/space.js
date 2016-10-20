var S={}
// var _Set = require("./set");

// 各種 space : https://en.wikipedia.org/wiki/Space_(mathematics)#/media/File:Mathematical_implication_diagram_eng.jpg

S.Space={} // Space = HllbertSpace + BanachSpace + Manifold (流形)

S.HilbertSpace={} // HilbertSpace => InnerProductSpace => LocallyConvaxSpace => VectorSpace (LinearSpace)

S.InnerProductSpace={}

S.LocallyConvaxSpace={}

S.LinearSpace=S.VectorSpace={} // (Algebraic)

S.BanachSpace={} // BanachSpace => NormedVectorSpace => MetricSpace => TopologicalSpace

S.NormedVectorSpace={}

S.MetricSpace={}

S.TopologicalSpace={} // (Analytic)

S.Manifold={}

S.Rn = {} // (R,R,....)


module.exports=S;