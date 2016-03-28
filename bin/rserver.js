#!/usr/bin/env node
var express = require('express');
var path = require('path');
var serveIndex = require('serve-index');
var app = express();
var dir = "../web";
// var dir = path.join(__dirname, '/../web');
app.use('/', express.static(dir));
// app.use('/', express.static(dir));
app.use('/', serveIndex(dir, {'icons': true}));

app.listen(3000);

console.log('Server started: http://localhost:3000/');
console.log('__dirname=', __dirname);