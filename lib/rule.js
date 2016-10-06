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

// relation
eq=R.eq=function(a,b)  {
  return (typeof a === typeof b) && a.toString()===b.toString() 
}
R.leq=function(a,b)  { return a<=b }
R.geq=function(a,b)  { return a>=b }
R.le =function(a,b)  { return a<b  }
R.ge =function(a,b)  { return a>b  }
// ========= type checking =================
R.yes=function(a) { return true }
R.no=function(a) {return false }
R.isFloat=function(a) { return typeof a==='number'}
R.isInteger=function(a) { return typeof a==='number' && a%1===0 }
R.isZero=function(a) { return a===0 }
R.isPositive=function(a) { return a>0 }
R.isNegative=function(a) { return a<0 }
R.isEven=function(a) { return (R.isInteger(a) && a%2===0) }
R.isOdd=function(a) { return (I.isInteger(a) && a%2===1) }
R.isFunction=function(a) { return typeof a==='function' }
R.isString=function(a) { return typeof a==='string' }
R.isObject=function(a) { return typeof a==='object' }
R.isArray=function(a) { return a.constructor===Array }
R.isSet=function(a) { return a.constructor===Set }

module.exports = R;
