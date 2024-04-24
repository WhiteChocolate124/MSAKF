#version 330
#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D src;
uniform float p1;

in vec2 tex_coords;
out vec4 frag_color;

//Modified Sobel Filter
vec4 StructTensor(sampler2D tex, vec2 Coord){
  vec2 uv = Coord;
  vec4 c = texture2D(tex, uv);
  vec2 resolution = textureSize(tex, 0);
  vec2 d = 1.0 / resolution;

  vec3 u = (
          +p1 * texture2D(tex, uv + vec2(-d.x, -d.y)).xyz +
          (1 - 2 * p1) * texture2D(tex, uv + vec2(-d.x,  0.0)).xyz +
          +p1 * texture2D(tex, uv + vec2(-d.x,  d.y)).xyz +
          -p1 * texture2D(tex, uv + vec2( d.x, -d.y)).xyz +
          (2 * p1 - 1) * texture2D(tex, uv + vec2( d.x,  0.0)).xyz +
          -p1 * texture2D(tex, uv + vec2( d.x,  d.y)).xyz
          ) / 2.0;
  vec3 v = (
          +p1 * texture2D(tex, uv + vec2(-d.x, -d.y)).xyz +
          (1 - 2 * p1) * texture2D(tex, uv + vec2( 0.0, -d.y)).xyz +
          +p1 * texture2D(tex, uv + vec2( d.x, -d.y)).xyz +
          -p1 * texture2D(tex, uv + vec2(-d.x,  d.y)).xyz +
          +(2 * p1 - 1) * texture2D(tex, uv + vec2( 0.0,  d.y)).xyz +
          -p1 * texture2D(tex, uv + vec2( d.x,  d.y)).xyz
          ) / 2.0;

  vec4 g = vec4(vec3(dot(u, u), dot(v, v), dot(u, v)), 1.);
  
  return g;
}

void main(){
    frag_color = StructTensor(src, tex_coords);
}