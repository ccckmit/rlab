var D = {	dx:0.001 };

D.differential = function(f, x) {
  var dy = f(x+D.dx) - f(x);
  return dy/D.dx;
}

D.d = D.differential;

D.integral = function(f, a, b) {
  var dx = D.dx;
  var area = 0.0;
  for (var x=a; x<b; x=x+D.dx) {
    area = area + f(x)*D.dx;
  }
  return area;
}

D.i = D.integral;

module.exports = D;