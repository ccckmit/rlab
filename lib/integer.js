var I = {};

I.gcd = function(a, b) {
  if (!b) return a;
  return I.gcd(b, a % b);
}

I.lcm = function(a, b) {
  return (a * b) / gcd(a, b);   
}

I.isPrime=function(n) {
	for (var i=2; i<=Math.sqrt(n); i++)
		if (n%i===0) return false;
	return n%1===0;
}

module.exports = I;