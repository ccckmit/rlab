var R = require("../rlab");
var s = R.S.sym;

var changeLabel = function(variable) {
    variable.label += "s";
    return variable;
};

var E = s("E = m*\\c**2");
log("E=MC2\nstr:", E.toString());
log(" latex:", E.toLaTeX());
log(" nmathml:", E.toMathML());
log("E.simplify()=", E.simplify().str()); 
log("E.approx()=", E.approx().str()); 
log("E.getAllVariables()=", E.getAllVariables().str()); 
log("E.mapOverVariables(changeLabel)=", E.mapOverVariables(changeLabel).str());  // 有錯
log("E.copy()=", E.copy().str());
log("E.getUncertainty()=", E.getUncertainty().str());
log("E.differentiate('m')=", E.differentiate("m").str());

var f = s("a*x**2 + b*x + c - y").toFunction("x", "y");
log("(a*x**2 + b*x + c - y).toFunction(x, y)(3,2)=", f(3,2).str());
