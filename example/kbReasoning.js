var fs = require('fs'); // 引用檔案物件
var kb = require('../lib/kb');
var kb1 = new kb();
var code = fs.readFileSync(process.argv[2], "utf8").replace(/\n/gi, "");
kb1.load(code);
kb1.forwardChaining();
