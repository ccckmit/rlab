var D = {	dx:0.001 };

// diff 和 df 都不能用
//   df 在 statistics 中定義過，f 分布的機率密度函數。
//   diff 則在 jStat 當中定義過
// print('R.diff=', R.diff); 
// print('R.df=', R.df);

D.differential = function(f, x) {
  var dy = f(x+D.dx) - f(x);
  return dy/D.dx;
}

D.integral = function(f, a, b) {
  var dx = D.dx;
  var area = 0.0;
  for (var x=a; x<b; x=x+D.dx) {
    area = area + f(x)*D.dx;
  }
  return area;
}

module.exports = D;