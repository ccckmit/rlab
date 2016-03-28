// 關鍵：只要 chain 下去，就要用 value() 取出來。
// 否則可以直接呼叫！

/*

_.range(5)
[0, 1, 2, 3, 4]

_.chain(5).range().value()
[0, 1, 2, 3, 4]

_(5).range().value()
[0, 1, 2, 3, 4]

_(5).range().map((x)=>x*x).value()
[0, 1, 4, 9, 16]

_.range(5).map(function() { return _.sample(space) })
[5, 6, 2, 2, 3]
_.chain(5).range().map(function() { return _.sample(space) })
LodashWrapper {__wrapped__: 5, __actions__: Array[2], __chain__: true, __index__: 0, __values__: undefined}
_.chain(5).range().map(function() { return _.sample(space) }).value()
[3, 6, 1, 4, 4]

_.map(_.range(5),  _.ary((x)=>x*x, 1))
[0, 1, 4, 9, 16]
_.chain().map(_.range(5),  _.ary((x)=>x*x, 1))
LodashWrapper {__wrapped__: undefined, __actions__: Array[1], __chain__: true, __index__: 0, __values__: undefined}
_.chain(_.range(5)).map(_.ary((x)=>x*x, 1))
LodashWrapper {__wrapped__: Array[5], __actions__: Array[1], __chain__: true, __index__: 0, __values__: undefined}
_.chain(_.range(5)).map(_.ary((x)=>x*x, 1)).value()
[0, 1, 4, 9, 16]

*/

var _ = require("lodash");
var c = console;

var str = function(x, sp) {
	return _.join(x, sp);
}

// _.mixin({str:str}, {chain:true});
_.mixin({str:str}, {chain:false});

var a = [1,2,3];
var sq = _([1,2,3]).chain().map(function(x) {return x*x;});
c.log("join=", sq.join('-'));
c.log("join=", sq.join('-').value());
c.log("str=", sq.str('-').value());
c.log("a.join=", _.join(a, '-'));

// var s = a.str().value();
var s = sq.str();
c.log(s);

/*
function vowels(string) {
  return _.filter(string, function(v) {
    return /[aeiou]/i.test(v);
  });
}

_.mixin({ 'vowels': vowels }, { 'chain': false });
c.log(_('fred').vowels());
*/
// var s = a.log("a=").value();
// s.value();
// c.log(s);
// c.log(s.value());

/*
var log = function(x, prompt) {
	c.log(prompt, x);
//	return "prompt"+x;
}

// _.mixin({log:log}, {chain:true});
_.mixin({log:log}, {chain:false});
*/
