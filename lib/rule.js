var R = {
	throwError:true,
}

// check : assert
var check = R.check = function(cond, msg) { 
	if (cond)
		console.log("O:"+msg);
	else {
		console.log("X:"+msg);
		if (R.throwError) throw Error('check fail!');
	}
}

R.be =function(msg,cond) { return check(cond, msg) }

R.proto=function(o) { return Object.getPrototypeOf(o) }

// relation
eq=R.eq=function(a,b)  {
  return (typeof a === typeof b) && a.toString()===b.toString() 
}
R.neq=function(a,b)  { return !R.eq(a,b) }
R.leq=function(a,b)  { return a<=b }
R.geq=function(a,b)  { return a>=b }
R.lt =function(a,b)  { return a<b  }
R.gt =function(a,b)  { return a>b  }
// ========= type checking =================
R.yes=function(a) { return true }
R.no=function(a) {return false }
R.isBool=function(a) { 
  return typeof a === 'boolean' || a instanceof Boolean 
}
R.isFunction=function(a) { 
  return typeof a==='function' || a instanceof Function 
}
R.isString=function(a) { 
  return typeof a==='string' || a instanceof String 
}
R.isObject=function(a) { 
  return typeof a==='object' || a instanceof Object 
}
R.isArray=function(a) { 
  return a instanceof Array 
}
R.isUndefined=function(a) { 
  return typeof a === 'undefined' 
}
R.isSet=function(a) { 
  return a instanceof Set 
}
R.isFloat=R.isNumber=function(a) { 
  return typeof a === 'number' || a instanceof Number 
}
R.isInteger=function(a) { return R.isFloat(a) && a%1===0 }
R.isZero=function(a) { return a===0 }
R.isPositive=function(a) { return a>0 }
R.isNegative=function(a) { return a<0 }
R.isEven=function(a) { return (R.isInteger(a) && a%2===0) }
R.isOdd=function(a) { return (I.isInteger(a) && a%2===1) }

module.exports = R;
