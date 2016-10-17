var R = require("../rlab");
var S = R.Set;
var be = R.Rule.be;
be('3.5 in Float', S.Float.has(3.5));
be('3.5 not in Integer', !S.Integer.has(3.5));
be('3 not in Even', !S.Even.has(3));
be('3 in Odd', S.Odd.has(3));
be('6 not in Prime', !S.Prime.has(6));
be('5 in Prime', S.Prime.has(5));

var A = new S.Set([1,2,3,4,5]);
var B=new S.Set([1,3,5,7]);
print("A=", A);
print("B=", B);
be("3 in A", A.has(3));
be("6 not in A", !A.has(6));

var UAB = S.union(A,B);
be('3 in A∪B:', UAB.has(3));
be('8 not in A∪B:', !UAB.has(8));
var IAB = S.intersection(A,B);
be('3 in A∩B:', IAB.has(3));
be('2 not in A∩B', !IAB.has(2));
var DAB = S.difference(A,B);
be('3 not in A-B', !DAB.has(3));
be('2 in A-B', DAB.has(2));
var SAB = S.symmetricDifference(A,B);
be('3 not in A(-)B', !SAB.has(3));
be('2 in A(-)B', SAB.has(2));
var EP = S.setUnion([S.Even, S.Prime]);
be('3 in Even∪Prime', EP.has(3));
be('8 in Even∪Prime', EP.has(8));
be('9 not in Even∪Prime', !EP.has(9));

be('2 not in Odd', !S.Odd.has(2));
be('3 in Odd', S.Odd.has(3));
var CEO = S.cartesianProduct(S.Even,S.Odd);
be('[2,3] in Even×Odd', CEO.has([2,3]));
be('[2,4] not in Even×Odd', !CEO.has([2,4]));
var P = S.powerSet(A);
be('{1,4} in P5', P.has(new S.Set([1,4])));

be('subset({1,3},A)', S.subSet(new S.Set([1,3]), A));
be('not subset({1,7},A)', !S.subSet(new S.Set([1,7]), A));

// Russell's Paradox, RP={ x | !x.has(x) }
var A=new Set([1,2,3]);
var RP=S.RussellSet;
be('not A.has(A)', !A.has(A));
be('RP.has(A)', RP.has(A));
var B = A.add(A);
be('!RP.has(B)', !RP.has(B));
be('RP.has(RP)', RP.has(RP)); // not halt, stack overflow

