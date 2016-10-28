// http://blog.smartbear.com/testing/four-serious-math-libraries-for-javascript/

// Four Serious Math Libraries for JavaScript

var B = {};

B.slice = function(a) {
	return Array.prototype.slice.call(a);
}

B.bind=function(o, member) {
	if (typeof o[member]==='Function')
		return o[member].bind(o); 
	else
		return o[member];
}

B.ncall = function() {
  var args = B.slice(arguments);
  var n = args[0];
  var o = args[1];
  var fname = args[2];
  var params = args.slice(3);
  var a=[];
  for (var i=0; i<n; i++)
    a.push(o[fname].apply(o, params));
  return a;
}

B.mapFunctions=function(host, obj, pairs) {
	for (var h in pairs) {
		var o = pairs[h];
		if (typeof host[h] !=='undefined')
			console.log('mapBind: error!', h, ' has been defined!');
		host[h]=B.bind(obj, o);
	}
}

B.copyFunctions=function(host, obj, names) {
	for (var name of names) {
		if (typeof host[name] !=='undefined')
			console.log('namesBind: error!', name, ' has been defined!');
		host[name]=B.bind(obj, name);
	}
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

B.arg1this = function(f,obj) { // 傳回一個已經綁定 f, obj 的函數
  return function() { 
		var args = B.slice(arguments);
    return f.apply(obj, [this].concat(args)); // 效果相當於 obj.f(this, args)
  }
}

B.mixThis=function(proto, obj, fmembers) {
	for (var fname of fmembers) {
		var f = obj[fname];
		if (typeof proto[fname] === 'undefined') {
			Object.defineProperty(proto, fname, {
				enumerable: false,
				writable: true,
				value: B.arg1this(f,obj), // proto.f(args) => obj.f(this, args) , 這行盡量不要動，除非想得很清楚了！
			});
		} else {
	    console.log("B.mixThis:", fname, " fail!");
		}
	}
}

B.mixThisMap=function(proto, obj, poMap) {
	for (var pname in poMap) {
		var oname = poMap[pname];
		var f = obj[oname];
		if (typeof proto[pname] === 'undefined') {
			Object.defineProperty(proto, pname, {
				enumerable: false,
				writable: true,
				value: B.arg1this(f,obj), // proto.f(args) = f(this, args) , 這行盡量不要動，除非想得很清楚了！
			});
		} else {
			console.log('pname=', pname, 'proto[pname]=', proto[pname]);
	    console.log("B.mixThisMap:", oname, " fail!");
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

B.def = function(x, value) {
	return (typeof x !== 'undefined')?x:value;
}

B.precision=2;

B.nstr = function(n, precision=B.precision) {
	if (n % 1 === 0) return n.toString();
	return n.toFixed(precision);
}

B.astr = function(a, precision=B.precision) {
	var s=[];
	for (var i in a) {
		s.push(a[i].str(precision));
	}
	return "["+s.join(', ')+"]";
}

B.sstr = function(s) { return s.toString(); }

B.ostr = function(o, precision=B.precision) {
  var s = [];
  for (var k in o)
    s.push(k+":"+B.str(o[k], precision));
  return "{"+s.join(", ")+"}";
}

B.str = function(o, precision=B.precision) {
	if (typeof o ==='undefined')
		return 'undefined';
	else
		return o.str(precision);
}

module.exports = B;

/*
B.calls = function() {
  var args = B.slice(arguments);
  var n = args[0];
  var f = args[1];
  var params = args.slice(2);
  var a=[];
  for (var i=0; i<n; i++)
    a.push(f.apply(null, params));
  return a;
}
*/
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
B.thisAsArg1 = function(f) {
  return function() {
		var args = B.slice(arguments);
    return f.apply(null, [this].concat(args));
  }
}

*/

/*
B.mixThis=function(proto, obj, fmembers) {
	for (var fname of fmembers) {
//		var f = B.bind(obj, fname);
		var f = obj[fname];
		if (typeof proto[fname] === 'undefined') {
			Object.defineProperty(proto, fname, {
				enumerable: false,
				writable: true,
				value: B.arg1this(f,obj), // proto.f(args) = f(this, args)
			});
		} else {
	    console.log("B.mixThis:", fname, " fail!");
		}
	}
}
*/