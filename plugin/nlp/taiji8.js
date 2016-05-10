var R = require("./R");
var c = console;

c.log("lpad=", String.prototype.lpad);
c.log(R.lpad("123", 5, '0'));
c.log("123".lpad(5, '0'));

c.log("R.rep([2,3],1)=", R.rep([2,3],1));

function gen(a,b,answerMap) {
  var m=[];
	for (var i=0; i<a.length;i++) {
	  m[i] = [];
	  for (var j=0; j<b.length; j++) {
	    var ab = a[i]+b[j];
	    var answer = answerMap[ab];
	    var abit=(typeof answer === 'undefined')?0:answer; 
      m[i][j] = ab+abit;
	  }
	}
	return m;
}

function str2bits(s) {
	var bits = [];
  for (var i=0; i<s.length; i++) {
  	if (s[i].match(/[01]/)) 
  		bits.push(s[i]);
  	else
  		bits.push(s.charCodeAt(i).toString(2).lpad(16,'0'));
  }
  return bits.join('');
}


var answerMap = { "隻貓":1, "個人":1, "條魚":1};

var m = gen("隻個條", "人魚貓", answerMap);
console.log("m=", m);

var a = R.flatM(m);
console.log("a=", a);

var abits = a.map(str2bits);
console.log("abits=", abits);

