// @ts-check

/**
 * @param {number[]|number[][]} matrixArray 
 * @param {boolean} columnMajor 
 */
 function Matrix(matrixArray, columnMajor=false) {
    // Must be a square matrix

    // Check if it is a multi dim array
    var multiDim = (typeof matrixArray[0] == "object");
    var mat = [];// saved in row major
    var dim;

    if (!multiDim) {
        dim = Math.sqrt(matrixArray.length);
        if (dim % 1 != 0) {
            throw "Matrix must be square";
        }

        for (var row = 0; row < dim; row++) {
            mat.push([]);
            for (var col = 0; col < dim; col++) {
                mat[row].push(matrixArray[row * dim + col]);
            }
        }

    } else {
        dim = matrixArray.length;
        //@ts-ignore // It says there's an error since matrixArray[0] can be a number or an array, but .length is not on numbers
        // but we ensure that we don't access .length unless it is an array
        if (dim != matrixArray[0].length) {
            throw "Matrix must be square";
        }
        for (var row = 0; row < dim; row++) {
            mat.push([]);
            for (var col = 0; col < dim; col++) {
                mat[row].push(matrixArray[row][col]);
            }
        }
    }
    // Rearange the matrix if it was given in column major order
    if (columnMajor) {
        // row and column refer to mat's format
        for (var row = 0; row < dim; row++) {
            for (var col = 0; col < dim; col++) {
                if (multiDim) {
                    mat[row][col] = matrixArray[col][row];
                } else {
                    mat[row][col] = matrixArray[col * dim + row];
                }
            }
        }
    }

    /**
     * @param {number[][]} matA 
     * @param {number[][]} matB 
     * @param {number} row 
     * @param {number} col 
     * @returns {number}
     */
    function dot(matA, matB, row, col) {
        // Calculates the dot product of a row and a col of 2 matricies
        // Used in the multiply function
        var result = 0;
        for (var i = 0; i < dim; i++) {
            result += matA[row][i] * matB[i][col];
        }
        return result;
    }


    /**
     * @param {Matrix} other 
     */
    this.multiply = function (other) {
        if (dim != other.getDim()) {
            throw "Matrix dimensions do not match";
        }

        /**@type {number[][]} */
        var result = [];
        var oMat = other.getMultiArray();
        for (var row = 0; row < dim; row++) {
            result.push([]);
            for (var col = 0; col < dim; col++) {
                result[row].push(dot(mat, oMat, row, col));
            }
        }

        return new Matrix(result);
    }


    this.getDim = function () {
        return dim;
    }

    /**
     * @param {Matrix} other 
     */
    this.add = function (other) {
        var oDim = other.getDim();
        if (dim != oDim) {
            throw "Cannot add matricies of differing dimensions";
        }

        var result = other.getMultiArray();
        for (var row = 0; row < dim; row++) {
            for (var col = 0; col < dim; col++) {
                result[row][col] += mat[row][col];
            }
        }

        return new Matrix(result);
    }

    /**
     * @param {number} value 
     * @returns Matrix
     */
    this.scale = function (value) {
        var result = Matrix.getIdenity(dim).getMultiArray();

        for (var row = 0; row < dim; row++) {
            for (var col = 0; col < dim; col++) {
                result[row][col] = value * mat[row][col];
            }
        }

        return new Matrix(result);
    }


    this.determinant = function () {
        // calculates the determinant
        // We will do this recursively
        if (dim == 2) {
            return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
        }

        if (dim < 2) {
            // Not sure why you need a 1x1 matrix but it is possible
            return 0;
        }

        // we will use the first row as our base
        var result = 0;

        // Build a 2d array that is 1 less than our current width
        var detMatArray = Matrix.getIdenity(dim - 1).getMultiArray();

        for (var col = 0; col < dim; col++) {
            for (var subCol = 0; subCol < dim; subCol++) {
                if (subCol == col) {
                    continue;
                }
                var detArrayCol = subCol < col ? subCol : subCol - 1;
                for (var subRow = 1; subRow < dim; subRow++) {
                    detMatArray[subRow - 1][detArrayCol] = mat[subRow][subCol];
                }
            }
            var detSign = col % 2 == 0 ? 1 : -1;
            result += detSign * mat[0][col] * new Matrix(detMatArray).determinant();

        }

        return result;
    }

    this.inverse = function () {
        var det = this.determinant();
        if (det != 0) {
            return new Matrix(this.cofactors().scale(1 / det).getMultiArray(true));
        } else {
            return null;
        }
    }

    this.minors = function () {
        var result = Matrix.getIdenity(dim).getMultiArray();
        var subMat = Matrix.getIdenity(dim - 1).getMultiArray();
        for (var row = 0; row < dim; row++) {
            for (var col = 0; col < dim; col++) {
                for (var subCol = 0; subCol < dim; subCol++) {
                    if (subCol == col) {
                        continue;
                    }
                    var subMatCol = subCol < col ? subCol : subCol - 1;
                    for (var subRow = 0; subRow < dim; subRow++) {
                        if (subRow == row) {
                            continue;
                        }
                        var subMatRow = subRow < row ? subRow : subRow - 1;
                        subMat[subMatRow][subMatCol] = mat[subRow][subCol];
                    }
                }
                result[row][col] = new Matrix(subMat).determinant();
            }
        }
        return new Matrix(result);
    }

    this.cofactors = function () {
        var result = Matrix.getIdenity(dim).getMultiArray();
        var subMat = Matrix.getIdenity(dim - 1).getMultiArray();
        for (var row = 0; row < dim; row++) {
            for (var col = 0; col < dim; col++) {
                for (var subCol = 0; subCol < dim; subCol++) {
                    if (subCol == col) {
                        continue;
                    }
                    var subMatCol = subCol < col ? subCol : subCol - 1;
                    for (var subRow = 0; subRow < dim; subRow++) {
                        if (subRow == row) {
                            continue;
                        }
                        var subMatRow = subRow < row ? subRow : subRow - 1;
                        subMat[subMatRow][subMatCol] = mat[subRow][subCol];
                    }
                }
                var detSign = (col + row) % 2 == 0 ? 1 : -1;
                result[row][col] = detSign * new Matrix(subMat).determinant();
            }
        }
        return new Matrix(result);
    }

    this.getMultiArray = function (columnMajor=false) {
        var result = [];
        if (columnMajor) {
            // row and column refer to mat's format
            for (var col = 0; col < dim; col++) {
                result.push([]);
                for (var row = 0; row < dim; row++) {
                    result[col].push(mat[row][col]);
                }
            }
        } else {
            // row and column refer to mat's format
            for (var row = 0; row < dim; row++) {
                result.push([]);
                for (var col = 0; col < dim; col++) {
                    result[row].push(mat[row][col]);
                }
            }
        }
        return result;
    }

    this.getArray = function (columnMajor=false) {
        var result = [];
        if (columnMajor) {
            // row and column refer to mat's format
            for (var col = 0; col < dim; col++) {
                for (var row = 0; row < dim; row++) {
                    result.push(mat[row][col]);
                }
            }
        } else {
            // row and column refer to mat's format
            for (var row = 0; row < dim; row++) {
                for (var col = 0; col < dim; col++) {
                    result.push(mat[row][col]);
                }
            }
        }
        return result;
    }
}

/**
 * @param {number} dim 
 * @returns {Matrix}
 */
Matrix.getIdenity = function (dim) {
    var result = [];
    for (var row = 0; row < dim; row++) {
        result.push([]);
        for (var col = 0; col < dim; col++) {
            if (col == row) {
                result[row].push(1);
            } else {
                result[row].push(0);
            }
        }
    }
    return new Matrix(result);
}

const MAX_SPARK_COUNT = 5000;
const MATRIX_BYTE_COUNT = 4 * 16;// 16 elements in a 4 by 4 matrix (4 bytes each element)

// arrays to be used for our vertex buffer
var matrixArray = new Float32Array(MAX_SPARK_COUNT * 16);

var lightArray = new Float32Array(MAX_SPARK_COUNT * 5);// red, green, blue, luminance, and intensity inputs per spark

var screenPosData = new Float32Array([
    -1.0, 1.0,
    1.0, -1.0,
    1.0, 1.0,
    -1.0, -1.0
]);
var screenFaceData = new Uint16Array([
    0, 1, 2,
    0, 3, 1
])

var modelArray = new Float32Array([
    0.5, 0.0,
    -0.25, 0.43,
    -0.25, -0.43,

    -0.5, 0.0,
    0.25, -0.43,
    0.25, 0.43
]);

var baseImage;
var lightImage;
var baseBuffer;

var blur1Image;
var blur2Image;
var blur1Buffer;
var blur2Buffer;


var textures = [];
var framebuffers = [];
/**
 * @type {WebGLVertexArrayObject}
 */
 var glowVAB;
 /**
 * @type {WebGLVertexArrayObject}
 */
  var blurVAB;
/**
 * @type {WebGLVertexArrayObject}
 */
    var combineVAB;
/**
 *  @type {WebGLBuffer}
 */
var matrixBuffer;
/**
*  @type {WebGLBuffer}
*/
var lightBuffer;
/**
 *  @type {WebGLBuffer}
 */
 var screenPosBuffer;
 /**
 *  @type {WebGLBuffer}
 */
  var screenFaceBuffer;
  /**
  *  @type {WebGLBuffer}
  */
var modelBuffer;
var glowProgram;
var glowVertSource = `#version 300 es
uniform mat4 world;
in mat4 model;
in vec3 color;
in vec2 pos;
in float intsy;
in float lum;
out vec3 _color;
out vec2 _pos;
out float _intsy;
out float _lum;
void main(){
    gl_Position = world * model * vec4(pos, 0.0, 1.0);
    _pos = pos;
    _lum = lum;
    _intsy = intsy;
    _color = color;
}`;
var glowFragSource = `#version 300 es

precision highp float;

in vec3 _color;
in vec2 _pos;
in float _intsy;
in float _lum;
layout(location = 0) out vec4 o_color0;
layout(location = 1) out vec4 o_color1;
float light;
void main(){
    light = pow(_lum, _intsy);
    o_color0 = vec4(_color * _lum, 1.0);
    o_color1 = vec4( _color * _color * _color * _color * light, 1.0);
}`;
var blurProgram;
var combineProgram;
var screenVertSource = `#version 300 es
in vec2 pos;
out vec2 f_uv;
void main(){
    f_uv = (pos + 1.0) / 2.0;// convert screen space to texture space
    gl_Position = vec4(pos, 0.0, 1.0);
}`
var blurFragSource = `#version 300 es
precision highp float;
in vec2 f_uv;
out vec4 o_color;
uniform sampler2D u_input;
uniform bool u_horizontal;
float weight[5] = float[] (0.227027, 0.1945946, 0.1216216, 0.054054, 0.016216);
vec2 offset;
void main(){
    float pixel_size = 2.0/float(textureSize(u_input, 0).x);
    if(u_horizontal){
        offset = vec2(pixel_size, 0.0);
    }else{
        offset = vec2(0.0, pixel_size);
    }

    vec3 color = texture(u_input, f_uv).rgb * weight[0];
    for(int i = 1; i < 5; i++){
        color += texture(u_input, f_uv + offset * float(i)).rgb * weight[i];
        color += texture(u_input, f_uv - offset * float(i)).rgb * weight[i];
    }
    o_color = vec4(color, 1.0);
}`;
var combineFragSource = `#version 300 es
precision highp float;
in vec2 f_uv;
out vec4 o_color;
uniform sampler2D u_base;
uniform sampler2D u_bloom;
void main(){
    vec3 color = texture(u_base, f_uv).rgb + texture(u_bloom, f_uv).rgb;
    o_color = vec4(color, 1.0);
}`;

var cameraMat = Matrix.getIdenity(4);
/**
 * @type {WebGLUniformLocation}
 */
 var glowWorldUniform;
 var blurInputUniform;
 var blurWeightUniform;
 var blurHorizontalUniform;
 var combineBaseUniform;
 var combineBloomUniform;


	// Camera Varibles
    var feildOfView = 90;// In degrees
    var camNear = 1;
    var camFar = 500;
    var camRot = 0;// Deg
    var camRotMat = Matrix.getIdenity(4);
    var antiRotMat = Matrix.getIdenity(4);
    
    function prepareCamera(){
        // Should only need to be called upon screen resize
        // Update the canvas's internal height and width
        // to match that of the screen
        
        var viewWidth = gl.drawingBufferWidth;
        var viewHeight = gl.drawingBufferHeight;

        // Calculate the slope of the view fulstrum
        var factor = Math.tan(Math.PI * 0.5 - 0.5 * feildOfView * Math.PI / 180);
        var rangeInv = 1 / (camNear - camFar);
        // Use the above values to compute a perspective matrix
        // Note we do not need to use the Matrix class to create this as
        // no operations will be done in javascript to this matrix
        var psm = [
            factor / (viewWidth / viewHeight), 0, 0, 0,
            0, factor, 0, 0,
            0, 0, (camNear + camFar) * rangeInv, -1,
            0, 0, camNear * camFar * rangeInv * 2, 0
        ];



        cameraMat = new Matrix(psm).multiply(new Matrix([1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, -200,
        0, 0 ,0, 1]));
        
    }

var canvas = document.createElement("canvas");

document.body.appendChild(canvas);

var width = canvas.clientWidth;
var height = canvas.clientHeight;

function updateCanvasSize(){
    width = canvas.clientWidth;
    height = canvas.clientHeight;
    canvas.width = width;
    canvas.height = height;
    gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
    // update buffer sizes
    prepareCamera();
}


window.addEventListener("resize", updateCanvasSize);

var gl = canvas.getContext("webgl2", {alpha:false});
if(gl == null || gl == undefined){
    throw "Could not start webgl2";
}

function startGL(){
    // Sets up webgl
    gl.clearColor(0.0, 0.0, 0.0, 1.0);
    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.LEQUAL);
    gl.clearDepth(1.0);
    gl.enable(gl.CULL_FACE);
    updateCanvasSize();

    glowVAB = gl.createVertexArray();
    blurVAB = gl.createVertexArray();
    combineVAB = gl.createVertexArray();
    matrixBuffer = gl.createBuffer();
    lightBuffer = gl.createBuffer();
    modelBuffer = gl.createBuffer();
    screenPosBuffer = gl.createBuffer();
    screenFaceBuffer = gl.createBuffer();


    glowProgram = buildProgram(gl, glowVertSource, glowFragSource);


    gl.useProgram(glowProgram);

    gl.bindVertexArray(glowVAB);
    glowWorldUniform = gl.getUniformLocation(glowProgram, "world");

    gl.bindBuffer(gl.ARRAY_BUFFER, modelBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, modelArray, gl.STATIC_DRAW);

    var gposLoc = gl.getAttribLocation(glowProgram, "pos");
    var gcolLoc = gl.getAttribLocation(glowProgram, "color");
    var gmatLoc = gl.getAttribLocation(glowProgram, "model");
    var glumLoc = gl.getAttribLocation(glowProgram, "lum");
    var gintsyLoc = gl.getAttribLocation(glowProgram, "intsy");

    gl.enableVertexAttribArray(gposLoc);
    gl.enableVertexAttribArray(gcolLoc);
    gl.enableVertexAttribArray(glumLoc);
    gl.enableVertexAttribArray(gintsyLoc);

    gl.vertexAttribPointer(gposLoc, 2, gl.FLOAT, false, 2 * 4, 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, matrixBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, MATRIX_BYTE_COUNT * MAX_SPARK_COUNT, gl.DYNAMIC_DRAW);

    for(var i = 0; i < 4; i++){
        var cLoc = gmatLoc + i;
        var offset = i * 16;
        
        gl.enableVertexAttribArray(cLoc);
        gl.vertexAttribPointer(cLoc, 4, gl.FLOAT, false, MATRIX_BYTE_COUNT, offset);
        gl.vertexAttribDivisor(cLoc, 1);
    }

    gl.bindBuffer(gl.ARRAY_BUFFER, lightBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, lightArray.byteLength, gl.DYNAMIC_DRAW);

    gl.vertexAttribPointer(gcolLoc, 3, gl.FLOAT, false, 5 * 4, 0);
    gl.vertexAttribDivisor(gcolLoc, 1);
    gl.vertexAttribPointer(glumLoc, 1, gl.FLOAT, false, 5 * 4, 3 * 4);
    gl.vertexAttribDivisor(glumLoc, 1);
    gl.vertexAttribPointer(gintsyLoc, 1, gl.FLOAT, false, 5 * 4, 4 * 4);
    gl.vertexAttribDivisor(gintsyLoc, 1);



    gl.bindVertexArray(blurVAB);
    blurProgram = buildProgram(gl, screenVertSource, blurFragSource);
    gl.useProgram(blurProgram);

    gl.bindBuffer(gl.ARRAY_BUFFER,screenPosBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, screenPosData, gl.STATIC_DRAW);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, screenFaceBuffer);
    gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, screenFaceData, gl.STATIC_DRAW);

    var bPos = gl.getAttribLocation(blurProgram, 'pos');

    gl.enableVertexAttribArray(bPos);
    gl.vertexAttribPointer(bPos, 2, gl.FLOAT, false, 8, 0);

    blurInputUniform = gl.getUniformLocation(blurProgram, 'u_input');
    blurHorizontalUniform = gl.getUniformLocation(blurProgram, 'u_horizontal');

    gl.bindVertexArray(combineVAB);
    combineProgram = buildProgram(gl, screenVertSource, combineFragSource);
    gl.useProgram(combineProgram);

    
    gl.bindBuffer(gl.ARRAY_BUFFER,screenPosBuffer);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, screenFaceBuffer);
    
    var cPos = gl.getAttribLocation(combineProgram, 'pos');

    gl.enableVertexAttribArray(cPos);
    gl.vertexAttribPointer(cPos, 2, gl.FLOAT, false, 8, 0);

    combineBaseUniform = gl.getUniformLocation(combineProgram, 'u_base');
    combineBloomUniform = gl.getUniformLocation(combineProgram, 'u_bloom');

    gl.bindVertexArray(null);


    gl.bindVertexArray(null);
    
    baseImage = buildTexture();
    blur1Image = buildTexture();
    blur2Image = buildTexture();
    lightImage = buildTexture();

    textures.push(baseImage, lightImage, blur1Image, blur2Image);

    baseBuffer = gl.createFramebuffer();
    blur1Buffer = gl.createFramebuffer();
    blur2Buffer = gl.createFramebuffer();

    gl.bindFramebuffer(gl.FRAMEBUFFER, baseBuffer);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, baseImage, 0);
    gl.framebufferTexture2D(gl.DRAW_FRAMEBUFFER, gl.COLOR_ATTACHMENT1, gl.TEXTURE_2D, lightImage, 0);
    gl.drawBuffers([gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1]);

    gl.bindFramebuffer(gl.FRAMEBUFFER, blur1Buffer);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, blur2Image, 0);
    
    gl.bindFramebuffer(gl.FRAMEBUFFER, blur2Buffer);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, blur1Image, 0);
    
    framebuffers.push(baseBuffer, blur1Buffer, blur2Buffer);

    gl.bindFramebuffer(gl.FRAMEBUFFER, null);


    prepareCamera();
}

function clear(){
    
    gl.clearColor(0, 0, 0, 0);
    framebuffers.forEach((val) => {
        gl.bindFramebuffer(gl.FRAMEBUFFER, val);
        gl.clear(gl.DEPTH_BUFFER_BIT | gl.COLOR_BUFFER_BIT);
    })
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.clearColor(0, 0, 1, 1)
    gl.clear(gl.DEPTH_BUFFER_BIT | gl.COLOR_BUFFER_BIT);
}

function buildTexture(){
    var texture = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, texture);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.drawingBufferWidth, gl.drawingBufferHeight, 0, gl.RGBA, gl.UNSIGNED_BYTE, null)
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return texture;
}

function buildBuffer(){
    var buffer = gl.createBuffer();
    var depthBuffer = gl.createRenderbuffer();

    gl.bindBuffer(gl.FRAMEBUFFER, buffer);
    gl.bindRenderbuffer(gl.RENDERBUFFER, depthBuffer);
    gl.renderbufferStorage(gl.RENDERBUFFER, gl.DEPTH_COMPONENT16, gl.drawingBufferWidth, gl.drawingBufferHeight);
    gl.framebufferRenderbuffer(gl.FRAMEBUFFER, gl.DEPTH_ATTACHMENT, gl.RENDERBUFFER, depthBuffer);
    return buffer;
}

function blur(times){
    
    gl.bindVertexArray(blurVAB);
    gl.useProgram(blurProgram);
    
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, blur1Image);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, blur2Image);
    
    gl.activeTexture(gl.TEXTURE2);
    gl.bindTexture(gl.TEXTURE_2D, lightImage);

    for(var i =0; i < times; i++){

        gl.uniform1i(blurHorizontalUniform, 0);
        if(i == 0){
            gl.uniform1i(blurInputUniform, 2);
        }else{
            gl.uniform1i(blurInputUniform, 1);
        }
        gl.bindFramebuffer(gl.FRAMEBUFFER, blur2Buffer);

        gl.drawElements(gl.TRIANGLES, 6, gl.UNSIGNED_SHORT, 0);

        gl.uniform1i(blurHorizontalUniform, 1);
        gl.uniform1i(blurInputUniform, 0);
        gl.bindFramebuffer(gl.FRAMEBUFFER, blur1Buffer);

        gl.drawElements(gl.TRIANGLES, 6, gl.UNSIGNED_SHORT, 0);
    }

    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
}

startGL();

canvas.addEventListener("webglcontextlost", function(event) {
    event.preventDefault();
    console.log("Webgl context lost");
}, false);

canvas.addEventListener(
    "webglcontextrestored", function(){startGL(); requestAnimationFrame(draw);}, false);

/**
 * @param {number} type
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} vx
 * @param {number} vy
 * @param {number} vz
 * @param {number} r
 * @param {number} g
 * @param {number} b
 * @param {number} lum
 * @param {number} insty
 * @param {number} time
 * @param {any} size
 * @param {any} tags
 */
function Spark(type, time, size, x, y, z, vx, vy, vz, r, g, b, lum, insty, tags, id=Spark.sparks.length){

    if(Spark.sparks.length >= MAX_SPARK_COUNT && id == Spark.sparks.length){
        // Don't create a new spark if we don't have room for it
        throw "To many sparks created";
    }
    // Add this spark to the list
    if(id == Spark.sparks.length){
        // detect if default value was used
        Spark.sparks.push(this);
    }else{
        // a specific id was given
        Spark.sparks[id] = this;
    }
    var lightArrayView = new Float32Array(lightArray.buffer, id * 5 * 4, 5);
    var matArrayView = new Float32Array(matrixArray.buffer, id * MATRIX_BYTE_COUNT, MATRIX_BYTE_COUNT / 4);
    var myMat = Matrix.getIdenity(4);

    function updateData(){
        var dataToLoad = myMat.getArray(true);
        for(var i = 0; i < 16; i++){
            matArrayView[i] = dataToLoad[i];
        }
        lightArrayView[0] = r;
        lightArrayView[1] = g;
        lightArrayView[2] = b;
        lightArrayView[3] = lum;
        lightArrayView[4] = insty;
    }

    this.isActive = function(){
        return isActive();
    }

    function isActive(){
        return time > 0;
    }

    this.isImportant = function(){
        switch(type){
            case Spark.types.ROCKET:
            case Spark.types.SHELL:
            case Spark.types.POP_SHELL:
                return true;
            default:
                return false;
        }
    }

    this.update = function(/** @type {number} */ deltaTime){
        if(!isActive()){
            myMat = new Matrix([0,0,0,0,
                0,0,0,0,
                0,0,0,0,
                0,0,0,0]);
                updateData();
            return;
        }
        time -= deltaTime;
        x += vx * deltaTime;
        y += vy * deltaTime;
        z += vz * deltaTime;

        //custom code for each type of spark
        switch(type){
            case Spark.types.ROCKET:
                if(Math.random() < 1){
                    var vec = angleToVector(getRandomAngle(Math.PI * 0.45, Math.PI * 0.55), 1);
                    var sr = vec.vy;
                    var sg = 0.5;
                    var sb = 0;
                    var nx = Math.random() * 30 - 15;
                    var ny = Math.random() * 30 - 15;
                    var nz = Math.random() * 30 - 15;
                    Spark.getNewSpark(Spark.types.EMBER, Math.random() * 2, 4, x, y, z, nx, ny, nz, sr, sg, sb, 0.9, 10, null);
                }
                
                vy -= deltaTime * 10;
                if(y > 0 && time > 2){
                    time = 0.5;
                }

                if(time < 2){
                    vx -= vx * deltaTime / 2;
                }

                if(time < 0.25){
                    r *= 0.9;
                    g *= 0.9;
                    b *= 0.9;
                    size *= 0.9;
                }

                if(time <= 0 && time + deltaTime > 0){
                    // did we just run out of time this frame?
                    // explode and spawn other sparks;
                    var newSparkCount = Math.floor(Math.random() * 50) + 100;
                    for(var i = 0; i < newSparkCount; i++){
                        var v1 = getRandomVector(Math.random() * 300);
                        var rand = Math.random() * 2 + 2;
                        var sr = r * rand;
                        var sg = g * rand;
                        var sb = b * rand;
                        sr %= 1;
                        sg %= 1;
                        sb %= 1;
                        var nx = v1.vx;
                        var ny = v1.vy;
                        var nz = v1.vz;
                        var nSpark = Spark.getNewSpark(Spark.types.SPARK, Math.random(), 10, x, y, z, nx, ny, nz, sr, sg, sb, 1, 15, null);
                        if(nSpark == null){
                            break;
                        } 
                    }
                }
                break;
            case Spark.types.SPARK:
                if(Math.random() < 1){
                    var sr = r;
                    var sg = g;
                    var sb = b;
                    var nx = Math.random() * 30 - 15;
                    var ny = Math.random() * 30 - 15;
                    var nz = Math.random() * 30 - 15;
                    Spark.getNewSpark(Spark.types.EMBER, Math.random(), 4, x, y, z, nx, ny, nz, sr, sg, sb, 0.9, 5, null);
                }
                vy -= deltaTime * 10;
                if(time < 0.25){
                    r *= 0.9;
                    g *= 0.9;
                    b *= 0.9;
                    size *= 0.9;
                }

                break;
            default:
                vy -= deltaTime * 10;
                if(time < 0.25){
                    r *= 0.9;
                    g *= 0.9;
                    b *= 0.9;
                    size *= 0.9;
                }
                break;
        
        }

        var transMat = new Matrix([
            size, 0, 0, x,
            0, size, 0, y,
            0, 0, size, z,
            0, 0, 0, 1]);
        
        myMat = transMat.multiply(antiRotMat);
    }

    this.render = function(){ // Temp function
        if(!isActive()){
            return;
        }
        updateData();
    }

}

function rgbToString(r, g, b){
    var rv = Math.floor(Math.min(Math.max(0, r), 1) * 255).toString(16);
    var gv = Math.floor(Math.min(Math.max(0, g), 1) * 255).toString(16);
    var bv = Math.floor(Math.min(Math.max(0, b), 1) * 255).toString(16);
    rv = rv.length < 2 ? "0" + rv : rv;
    gv = gv.length < 2 ? "0" + gv : gv;
    bv = bv.length < 2 ? "0" + bv : bv;
    return "#" + rv + gv + bv;  
}


/**
 * @type {Spark[]}
 */
Spark.sparks = [];

Spark.lastInactiveSpark = 0;

/**
 * @param {number} type
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} vx
 * @param {number} vy
 * @param {number} vz
 * @param {number} r
 * @param {number} g
 * @param {number} b
 * @param {number} lum
 * @param {number} insty
 * @param {number} time
 * @param {any} size
 * @param {any} tags
 */
Spark.getNewSpark = function(type, time, size, x, y, z, vx, vy, vz, r, g, b, lum, insty, tags){
    if(Spark.sparks.length < MAX_SPARK_COUNT){
        return new Spark(type, time, size, x, y, z, vx, vy, vz, r, g, b, lum, insty, tags);
    }else{
        var end = Spark.lastInactiveSpark;
        Spark.lastInactiveSpark ++;
        Spark.lastInactiveSpark %= MAX_SPARK_COUNT;

        // Search for an inactive spark to use
        while(Spark.lastInactiveSpark != end){
            if(Spark.sparks[Spark.lastInactiveSpark].isActive() == false){
                return new Spark(type, time, size, x, y, z, vx, vy, vz, r, g, b, lum, insty, tags, Spark.lastInactiveSpark);
            }
            Spark.lastInactiveSpark ++;
            Spark.lastInactiveSpark %= MAX_SPARK_COUNT;
        }

        // If we did not find one, find the next non important one
        Spark.lastInactiveSpark ++;
        Spark.lastInactiveSpark %= MAX_SPARK_COUNT;

        while(Spark.lastInactiveSpark != end){
            if(Spark.sparks[Spark.lastInactiveSpark].isImportant() == false){
                return new Spark(type, time, size, x, y, z, vx, vy, vz, r, g, b, lum, insty, tags, Spark.lastInactiveSpark);
            }
            Spark.lastInactiveSpark ++;
            Spark.lastInactiveSpark %= MAX_SPARK_COUNT;
        }

        // If we can't find one to replace, return nothing, there is not much we can do here
        return null;
    }
}

Spark.types = {
    ROCKET: 0,
    EMBER: 1,
    SHELL: 2,
    SPARK: 3,
    POP_SHELL: 4,
    POP_SPARK: 5
};

function getRandomAngle(start=0, end=2 * Math.PI){
    return (end - start) * Math.random() + start;
}


/**
 * @param {number} angle
 * @param {number} mag
 */
function angleToVector(angle, mag){
    return {vx:Math.cos(angle) * mag, vy:Math.sin(angle) * mag};
}

function anglesToVector(angle1, angle2, mag){
    return {vx:Math.cos(angle1) * mag, vy:Math.sin(angle1) * mag, vz:Math.cos(angle2) * mag};
}

/**
 * @param {number} mag
 */
function getRandomVector(mag){
    
    return anglesToVector(getRandomAngle(), getRandomAngle(), mag);
}

/**
		 * @param {WebGLRenderingContext} gl 
		 * @param {string} vertexShaderSource 
		 * @param {string} fragmentShaderSource
		 */
// builds our webgl program
 function buildProgram(gl, vertexShaderSource, fragmentShaderSource){
			
    var vertShader = gl.createShader(gl.VERTEX_SHADER);
    gl.shaderSource(vertShader, vertexShaderSource);
    gl.compileShader(vertShader);
    if (!gl.getShaderParameter(vertShader, gl.COMPILE_STATUS)) {
        throw "ERROR IN VERTEX SHADER : " + gl.getShaderInfoLog(vertShader);
    }

    var fragShader = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(fragShader, fragmentShaderSource);
    gl.compileShader(fragShader);
    if (!gl.getShaderParameter(fragShader, gl.COMPILE_STATUS)) {
        throw "ERROR IN FRAG SHADER : " + gl.getShaderInfoLog(fragShader);
    }

    var program = gl.createProgram();
    gl.attachShader(program, vertShader);
    gl.attachShader(program, fragShader);
    gl.linkProgram(program);

    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        throw "Unknown error in program";
    }

   return program;
}

var oldTime = 0;

/**
 * @param {number} time
 */
function draw(time){
    if(gl.isContextLost()){
        return;
    }
    //update buffers on gpu side

    //
    if(oldTime == 0){
        oldTime = time - 1;
    }

    camRot += (time - oldTime) / 1000 * 0.5;

    camRotMat = cameraMat.multiply(new Matrix([
        Math.cos(camRot), 0, Math.sin(camRot), 0,
        0, 1, 0, 0,
        -Math.sin(camRot), 0, Math.cos(camRot), 0,
        0, 0, 0, 1]));
    antiRotMat = new Matrix([
        Math.cos(-camRot), 0, Math.sin(-camRot), 0,
        0, 1, 0, 0,
        -Math.sin(-camRot), 0, Math.cos(-camRot), 0,
        0, 0, 0, 1]);

    
    gl.useProgram(glowProgram);
    gl.uniformMatrix4fv(glowWorldUniform, false, camRotMat.getArray(true));
    clear();
    Spark.sparks.forEach((spark) => {
        spark.update((time - oldTime) / 1000);
        spark.render();
    })
    //Update buffers
    

    gl.bindBuffer(gl.ARRAY_BUFFER, matrixBuffer);
    gl.bufferSubData(gl.ARRAY_BUFFER, 0, matrixArray);
    
    gl.bindBuffer(gl.ARRAY_BUFFER, lightBuffer);
    gl.bufferSubData(gl.ARRAY_BUFFER, 0, lightArray);

    // Draw glow
    
	gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);


    gl.bindVertexArray(glowVAB);
    gl.useProgram(glowProgram);

    gl.bindFramebuffer(gl.FRAMEBUFFER, baseBuffer);
    gl.drawArraysInstanced(gl.TRIANGLES, 0, 6, MAX_SPARK_COUNT);

    blur(8);

    gl.bindVertexArray(combineVAB);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.useProgram(combineProgram);

    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, baseImage);
    gl.activeTexture(gl.TEXTURE1);
    gl.bindTexture(gl.TEXTURE_2D, blur1Image);

    gl.uniform1i(combineBaseUniform, 0);
    gl.uniform1i(combineBloomUniform, 1);

    gl.drawElements(gl.TRIANGLES, 6, gl.UNSIGNED_SHORT, 0);

    gl.bindVertexArray(null);

    if(Math.random() < 0.025){
        var color = getRandomVector(1);
        
        Spark.getNewSpark(Spark.types.ROCKET, 5, 10, Math.random() * 100 - 50, -height / 2, Math.random() * 100 - 50, Math.random() * 50 - 25, Math.random() * 500 + 50, Math.random() * 50 - 25, Math.abs(color.vx), Math.abs(color.vy), Math.abs(color.vz), 1, 20, null);

    }
    oldTime = time;
    requestAnimationFrame(draw);
}
requestAnimationFrame(draw);



