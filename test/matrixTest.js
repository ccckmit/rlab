var assert = require('chai').assert;
var M = require("../rlab").M;

describe('Matrix', function() {
  var A = [[1,2,3],[4,5,6],[7,3,9]];
	var I = M.identity(3);
  describe('test tr(), eq(), all()', function () {
    it('A=A.tr().tr()', function () {
      assert.equal(true, A.tr().tr().eq(A).all());
    });
  });
  describe('test inv(), dot(), sub(), lt()', function () {
    it('A-1*A=I', function () {
      assert.isTrue(A.inv().dot(A).sub(I).abs().sumM()<0.001);
    });
  });
  describe('test eig()', function () {
    it('E,λ=A.eig(); E*λ.diag()-A=0', function () {
			var Eλ = M.eigR(A);
			var E = Eλ.E, λ=Eλ.lambda;
			assert.isTrue(E.dot(λ.diag()).dot(E.inv()).sub(A).abs().sumM()<0.001);
    });
  });
});
