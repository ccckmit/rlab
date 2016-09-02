var R = require("../rlab");

// FFT: fast fourier transform
var z = (new R.complex([1,2,3,4,5],[6,7,8,9,10])).fft()
log('z=', z.strM());
log('z.ifft=', z.ifft().strM());