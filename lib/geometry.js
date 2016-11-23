var G = require("./statistics");

// 多邊形：Polygon
// 多面體：Polyhedron : V-E+F = 2
// 多胞形：Polytope
// ============= Spherical Geometry ====================
G.SphericalGeometry = {}

G.MandelbrotSet = {}
/*
For each pixel (Px, Py) on the screen, do:
{
  x0 = scaled x coordinate of pixel (scaled to lie in the Mandelbrot X scale (-2.5, 1))
  y0 = scaled y coordinate of pixel (scaled to lie in the Mandelbrot Y scale (-1, 1))
  x = 0.0
  y = 0.0
  iteration = 0
  max_iteration = 1000
  // Here N=2^8 is chosen as a reasonable bailout radius.
  while ( x*x + y*y < (1 << 16)  AND  iteration < max_iteration ) {
    xtemp = x*x - y*y + x0
    y = 2*x*y + y0
    x = xtemp
    iteration = iteration + 1
  }
  // Used to avoid floating point issues with points inside the set.
  if ( iteration < max_iteration ) {
    // sqrt of inner term removed using log simplification rules.
    log_zn = log( x*x + y*y ) / 2
    nu = log( log_zn / log(2) ) / log(2)
    // Rearranging the potential function.
    // Dividing log_zn by log(2) instead of log(N = 1<<8)
    // because we want the entire palette to range from the
    // center to radius 2, NOT our bailout radius.
    iteration = iteration + 1 - nu
  }
  color1 = palette[floor(iteration)]
  color2 = palette[floor(iteration) + 1]
  // iteration % 1 = fractional part of iteration.
  color = linear_interpolate(color1, color2, iteration % 1)
  plot(Px, Py, color)
}
*/


// 碎形幾何 : Fractal
// http://andrew-hoyer.com/
// http://andrew-hoyer.com/experiments/fractals/
// http://flam3.com/flame.pdf
// https://en.wikipedia.org/wiki/List_of_fractals_by_Hausdorff_dimension

// http://rembound.com/articles/drawing-mandelbrot-fractals-with-html5-canvas-and-javascript
// https://github.com/rembound/Mandelbrot-Fractal-HTML5/blob/master/mandelbrot-fractal.js


// 操控繪圖
// http://fabricjs.com/ (讚！)
// https://github.com/kangax/fabric.js

// 3D 動畫
// https://threejs.org/
// http://haptic-data.com/toxiclibsjs/

// 地理資訊
// ArcGIS : https://developers.arcgis.com/javascript/3/

// 動畫
// http://paperjs.org/features/
// https://processing.org/

// 向量
// http://victorjs.org/
// 3d: https://evanw.github.io/lightgl.js/docs/vector.html

module.exports=G;