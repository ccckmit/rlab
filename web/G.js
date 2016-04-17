// 考慮： https://plot.ly/javascript/3d-surface-plots/ 互動性好，有 3D
G = {};

// ------------------ 繪圖物件與函數 C3.js 部份  --------------------------------
// 注意： C3 的繪圖對像好像是放在 $$.config 裏，預設可能是 chart
var C3G=function() {
  this.g = {
    data: {
/*    x: "x",  */
      xs: {},
      columns: [ /*["x", 1, 2, 3, 4 ]*/ ],
      type: "line", 
      types : {}
    },
    axis: {
      x: {
        label: 'X',
        tick: { fit: false, format:d3.format(".2f") }
      },
      y: { label: 'Y', 
        tick : { format: d3.format(".2f") }
      }
    }, 
    bar: { width: { ratio: 0.9 } }, 
  };
  this.varcount = 0;
  this.xrange(-10, 10);
  this.step = 1;
  this.setrange = false;
}

C3G.prototype.tempvar = function() { return "T"+this.varcount++; }

C3G.prototype.xrange = function(xmin, xmax) {
  this.xmin = xmin;
  this.xmax = xmax;
  if (arguments.length == 0)
    this.setrange = false;
  else
    this.setrange = true;
}

C3G.prototype.plot = function(x,y, options) {
  var name = options.name || this.tempvar();
  this.g.data.types[name] = "scatter";
  this.g.data.xs[name] = name+"x";
  this.g.data.columns.push([name+"x"].concat(x));
  this.g.data.columns.push([name].concat(y));
}

C3G.prototype.curve = function(f, options) {
  var name = options.name|| this.tempvar();
  var step = options.step || this.step;
  var from = options.from || this.xmin;
  var to   = options.to || this.xmax;
  this.g.data.types[name] = "line";
  this.g.data.xs[name] = name+"x";
  var x = R.steps(from, to, step), y=[];
	var y = x.map(f);
/*	
	console.log("x=", x, "y=", y);
	console.log("f=", f, "curve:x=", x);
  for (var i in x) {
    if (typeof(f)==="string")
      y.push(eval(f.replace("x", x[i])));
    else if (typeof f === "function")
      y.push(f(x[i]));
  }
	console.log("y=", y);
*/	
  this.g.data.columns.push([name+"x"].concat(x));
  this.g.data.columns.push([name].concat(y));
}

C3G.prototype.hist = function(x, options) {
  var name = options.name || this.tempvar(); 
  var mode = options.mode || ""; 
  var step = options.step || this.step; 
  var from = options.from || this.xmin; 
  var to   = options.to   || this.xmax;
  this.g.data.types[name] = "bar";
  this.g.data.xs[name] = name+"x";
  var xc = R.steps(from+step/2.0, to, step);
  var n = (to-from)/step + 1;
  var count = R._.fill(Array(n), 0);
  for (var i in x) {
    var slot=Math.floor((x[i]-from)/step);
    if (slot>=0 && slot < n)
      count[slot]++;
  }
//	console.log("count=", count);
  this.g.data.columns.push([name+"x"].concat(xc));
  var total = R.sum(count);
  if (mode === "normalized")
    count = count.map(function(c) { return (1.0/step)*(c/total); });
//	console.log("count=", count);
  this.g.data.columns.push([name].concat(count));
}

C3G.prototype.show = function() {
  if (typeof(module)==="undefined")
    return c3.generate(this.g);
}

// ------------------ 繪圖物件與函數 vis.js 部份  --------------------------------
var VISG=function() {
      // specify options
  this.options = {
        width:  '95%',
        height: '95%',
        style: 'surface',
        showPerspective: true,
        showGrid: true,
        showShadow: false,
        keepAspectRatio: true,
        verticalRatio: 0.5
      };
  this.data = null; 
  this.graph= null;
}

VISG.prototype.curve3d = function(f, box) {
  // Create and populate a data table.
  var data = new vis.DataSet();
  // create some nice looking data with sin/cos
  var counter = 0;
  var steps = 50;  // number of datapoints will be steps*steps
  var axisMax = 314;
  var axisStep = axisMax / steps;
  for (var x = 0; x < axisMax; x+=axisStep) {
    for (var y = 0; y < axisMax; y+=axisStep) {
      var value = f(x,y);
      data.add({id:counter++,x:x,y:y,z:value,style:value});
    }
  }
  this.graph = new vis.Graph3d(box, data, this.options);
}

// ---------------------- 繪圖整合函數 ----------------------------------
G.newGraph = function(box) {
  G.c3g = new C3G();
  G.visg = new VISG();
  G.box = box;
//  G.c3g.show();
}

// C3 部份
G.xrange = function(xmin, xmax) { G.c3g.xrange(xmin, xmax); }

G.curve = function(f, options) {
  G.c3g.curve(f, options);
  G.c3g.show();
}

G.hist = function(x, options) {
  G.c3g.hist(x, options);
//	console.log("x=", x, "G=", G);
  G.c3g.show();
}

G.plot = function(x,y,options) {
  G.c3g.plot(x,y,options);
  G.c3g.show();
}

// Vis 部份
G.curve3d = function(f) {
  G.visg.curve3d(f, G.box);
}
