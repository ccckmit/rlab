var map={
	a:null,
	b:'',
	c:'a',
	d:undefined,
}

for (var k in map) {
	console.log(k, map[k]);
}

// map.forEach((value,key)=>console.log(value, key));