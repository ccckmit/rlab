var assert = require('chai').assert;
var S = require("../rlab");
var eq = assert.equal;

describe('Set', function() {
	var A = new S.Set([1,2,3,4,5]);
	var B = new S.Set([1,3,5,7]);
  describe('Set.has(e)', function () {
    it('Default Set', function () {
      eq(true, S.Float.has(3.5));
			eq(false, S.Integer.has(3.5));
			eq(false, S.Even.has(3));
			eq(true, S.Odd.has(3));
			eq(false, S.Prime.has(6));
			eq(true, S.Prime.has(5));
		})
    it('Set Operation', function () {
			var UAB = S.union(A,B);
      eq(true, UAB.has(3));
			eq(false, UAB.has(8));
			var IAB = S.intersection(A,B);
			eq(true, IAB.has(3));
			eq(false, IAB.has(2));
			var DAB = S.difference(A,B);
			print('3 in A-B:', DAB.has(3));
			print('2 in A-B:', DAB.has(2));
			
			eq(false, S.Even.has(3));
			eq(true, S.Odd.has(3));
			eq(false, S.Prime.has(6));
			eq(true, S.Prime.has(5));
		})
  })
})
/*
var UAB = A.union(B);
print('3 in A∪B:', UAB.has(3));
print('8 in A∪B:', UAB.has(8));
var IAB = S.intersection(A,B);
print('3 in A∩B:', IAB.has(3));
print('2 in A∩B:', IAB.has(2));
var DAB = S.difference(A,B);
print('3 in A-B:', DAB.has(3));
print('2 in A-B:', DAB.has(2));
var SAB = S.symmetricDifference(A,B);
print('3 in A(-)B:', SAB.has(3));
print('2 in A(-)B:', SAB.has(2));
var EP = S.setUnion([S.Even, S.Prime]);
print('3 in Even∪Prime', EP.has(3));
print('8 in Even∪Prime', EP.has(8));
print('9 in Even∪Prime', EP.has(9));

print('2 in Odd:', S.Odd.has(2));
print('3 in Odd:', S.Odd.has(3));
var CEO = S.cartesianProduct(S.Even,S.Odd);
print('[2,3] in Even×Odd:', CEO.has([2,3]));
print('[2,4] in Even×Odd:', CEO.has([2,4]));
var P = S.powerSet(A);
print('{1,4} in P5:', P.has(new S.Set([1,4])));

print('subset({1,3},A)=', S.subSet(new S.Set([1,3]), A));
print('subset({1,7},A)=', S.subSet(new S.Set([1,7]), A));

*/