// http://blog.smartbear.com/testing/four-serious-math-libraries-for-javascript/

// Four Serious Math Libraries for JavaScript

var B = {};

B.slice = function(a) {
	return Array.prototype.slice.call(a);
}

B.curry = function(f,o) {
  return function() {
		var args = Array.prototype.slice.call(arguments);
    return f.apply(null, [o].concat(args));
  }
}

B.curryThis = function(f) {
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

// options={thisify:, enumerable:, override:}
B.mix=function(o, map, options) {
	options = options || {};
	for (var name in map) {
		f = map[name];
		if (options.thisify) {
			f = B.curryThis(f);
		}
		if (typeof o[name] === 'undefined' || options.override) {
      if (options.override) {
		    B.print("B.mix:", name, " overrided!");
			}
			Object.defineProperty(o, name, {
				enumerable: options.enumerable || false,
				writable: true,
				value: f,
			});
		} else if (typeof o[name] !== 'undefined') {
	    B.print("B.mix:", name, " skipped! f=", f);
		}
	}
}

B.mixThis=function(o, map) { B.mix(o, map, {thisify:true}); }

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

module.exports = B;

