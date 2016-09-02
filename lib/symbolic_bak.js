// http://coffeequate.readthedocs.org/ (讚!)
// https://www.npmjs.com/package/coffeequate
// https://github.com/jiggzson/nerdamer
// https://github.com/aantthony/javascript-cas
var CQ = require("coffeequate");

var s1 = CQ("E = m*\\c**2");
s1.__proto__.str = s1.__proto__.toString;

var S = {
  sym:CQ,
};

module.exports = S;