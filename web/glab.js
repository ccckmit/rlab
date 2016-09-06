G = glab = {
	dgMap:{}, // (dimension+graph) Map
}

G.new2D=function() {
  return { // c3 graph
    data: {
      xs: {},
      columns: [], 
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
  }
}

G.show = function(chartName, dg) {
	var g = dg.graph;
	if (dg.dim === '2D') {
		g.bindto=chartName;
		return c3.generate(g);		
	} else {
		var box=document.getElementById(chartName.replace('#',''));
		new vis.Graph3d(box, g.dataSet, g.options);		
	}
}

G.chart2D = function(chartName, f) {
	var g = G.new2D();
	f(g);
	var dg = {dim:'2D', graph:g};
	G.dgMap[chartName] = dg;
	return G.show(chartName, dg);
}

// type : line, spline, step, area, area-spline, area-step, bar, scatter, pie, donut, gauge
G.draw = function(g, name, x, y, type) {
  g.data.types[name] = type;
  g.data.xs[name] = name+"x";
  g.data.columns.push([name+"x"].concat(x));
  g.data.columns.push([name].concat(y));
}

G.curve = function(g, name, f, from=-10, to=10, step=0.1) {
	var rg  = R.G.curve(f, from, to, step);
	G.draw(g, name, rg.x, rg.y, 'line');
}

G.hist = function(g, name, x, type, from, to, step=1) {
	var rh = R.G.hist(x, from, to, step);
	G.draw(g, name, rh.xc, rh.bins, type || 'bar');
}

G.ihist=function(g, name, x, type) {
	G.hist(g, name, x, type, x.min()-0.5, x.max()+0.5, 1);
}

G.plot =(g, name, x, y)=>G.draw(g, name, x, y, 'scatter');

G.pie =function(g, countMap) {
	g.data.type = 'pie';
	for (var name in countMap) {
		var count = countMap[name];
		g.data.columns.push([name, count]);
	}
}

G.timeSeries =function(g, columns) {
	g.data.x = 'x';
	g.axis.x = { type:'timeseries', tick:{format:'%Y-%m-%d'}};
	g.data.columns = columns;
}

G.new3D = function() {
  return {
		dataSet:new vis.DataSet(),
		options:{
      width:  '95%',
      height: '95%',
      style: 'surface',
      showPerspective: true,
      showGrid: true,
      showShadow: false,
      keepAspectRatio: true,
      verticalRatio: 0.5
    }
	}
}
/*
G.chart3D = function(chartName, f) {
	var g = G.new3D();
	f(g);
	var box=document.getElementById(chartName.replace('#',''));
  new vis.Graph3d(box, g.dataSet, g.options);
}

G.curve3D = function(g, f) {
	var g = G.new3D();
  // create some nice looking data with sin/cos
  var counter = 0;
  var steps = 50;  // number of datapoints will be steps*steps
  var axisMax = 314;
  var axisStep = axisMax / steps;
  for (var x = 0; x < axisMax; x+=axisStep) {
    for (var y = 0; y < axisMax; y+=axisStep) {
      var value = f(x,y);
      g.dataSet.add({id:counter++,x:x,y:y,z:value,style:value});
    }
  }
	var box=document.getElementById(chartName.replace('#',''));
  new vis.Graph3d(box, g.dataSet, g.options);
}
*/

// style : surface, grid, bar, bar-color, bar-size, dot, dot-line, dot-color, dot-size, line

G.chart3D = function(chartName, style, f) {
	var g = G.new3D();
  // create some nice looking data with sin/cos
  var counter = 0;
  var steps = 50;  // number of datapoints will be steps*steps
  var axisMax = 314;
  var axisStep = axisMax / steps;
  for (var x = 0; x < axisMax; x+=axisStep) {
    for (var y = 0; y < axisMax; y+=axisStep) {
      var value = f(x,y);
      g.dataSet.add({id:counter++,x:x,y:y,z:value,style:value});
    }
  }
	g.options.style=style;
	G.dgMap[chartName] = {dim:'3D',graph:g};
	var box=document.getElementById(chartName.replace('#',''));
  new vis.Graph3d(box, g.dataSet, g.options);
}
