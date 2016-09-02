// http://blog.smartbear.com/testing/four-serious-math-libraries-for-javascript/

// Four Serious Math Libraries for JavaScript

var B = {};
/*
B.slice = function(a) {
	return Array.prototype.slice.call(a);
}

B.curry = function(f,o) {
  return function() {
		var args = Array.prototype.slice.call(arguments);
    return f.apply(null, [o].concat(args));
  }
}
*/
B.thisAsArg1 = function(f) {
  return function() {
		var args = Array.prototype.slice.call(arguments);
    return f.apply(null, [this].concat(args));
  }
}

B.calls = function() {
  var args = Array.prototype.slice.call(arguments);
  var n = args[0];
  var f = args[1];
  var params = args.slice(2);
  var a=[];
  for (var i=0; i<n; i++)
    a.push(f.apply(null, params));
  return a;
}

B.mix=function(self, members) {
	for (var name in members) {
		var member = members[name];
		if (typeof self[name] === 'undefined') {
			Object.defineProperty(self, name, {
				enumerable: true,
				writable: true,
				value: member,
			});
		} else {
	    console.log("B.mix fail:", name, " already exists!");
		}
	}
}

B.arg1this = function(f) {
  return function() {
		var args = Array.prototype.slice.call(arguments);
    return f.apply(null, [this].concat(args));
  }
}

B.mixThis=function(proto, fmap) {
	for (var name in fmap) {
		var f = fmap[name];
		if (typeof proto[name] === 'undefined') {
			Object.defineProperty(proto, name, {
				enumerable: false,
				writable: true,
				value: B.arg1this(f), // proto.f(args) = f(this, args)
			});
		} else {
	    console.log("B.mixThis:", name, " fail!");
		}
	}
}

B.ctrim=function(s, set, side) {
  side = side||"both";
  for (var b=0; b<s.length; b++)
    if (set.indexOf(s[b])<0) break;
  for (var e=s.length-1; e>=0; e--)
    if (set.indexOf(s[e])<0) break;
  if (side === "right") b=0;
  if (side === "left") e=s.length-1;
  return s.substring(b, e+1);
}

B.lpad=function(s, width, ch) {
  return s.length >= width ? s : new Array(width - s.length + 1).join(ch) + s;
}

B.steps = function(from, to, step) {
	step=step || 1;
	var a=[];
	for (var t=from; t<=to; t+=step)
		a.push(t);
	return a;
}

B.def = function(x, value) {
	return (typeof x !== 'undefined')?x:value;
}

module.exports = B;

