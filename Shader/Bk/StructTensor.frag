#version 330
#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D src;

in vec2 tex_coords;
out vec4 frag_color;


vec4 StructTensor(sampler2D tex, vec2 Coord){
  vec2 uv = Coord;
  vec4 c = texture2D(tex, uv);
  vec2 resolution = textureSize(tex, 0);
  vec2 d = 1.0 / resolution;
  
  vec3 u = (
          -1.0 * texture2D(tex, uv + vec2(-d.x, -d.y)).xyz +
          -2.0 * texture2D(tex, uv + vec2(-d.x,  0.0)).xyz +
          -1.0 * texture2D(tex, uv + vec2(-d.x,  d.y)).xyz +
          +1.0 * texture2D(tex, uv + vec2( d.x, -d.y)).xyz +
          +2.0 * texture2D(tex, uv + vec2( d.x,  0.0)).xyz +
          +1.0 * texture2D(tex, uv + vec2( d.x,  d.y)).xyz
          ) / 4.0;
  vec3 v = (
          -1.0 * texture2D(tex, uv + vec2(-d.x, -d.y)).xyz +
          -2.0 * texture2D(tex, uv + vec2( 0.0, -d.y)).xyz +
          -1.0 * texture2D(tex, uv + vec2( d.x, -d.y)).xyz +
          +1.0 * texture2D(tex, uv + vec2(-d.x,  d.y)).xyz +
          +2.0 * texture2D(tex, uv + vec2( 0.0,  d.y)).xyz +
          +1.0 * texture2D(tex, uv + vec2( d.x,  d.y)).xyz
          ) / 4.0;
  
  vec4 g = vec4(vec3(dot(u, u), dot(v, v), dot(u, v)), 1.);
  
  return g;
}

void main(){
    frag_color = StructTensor(src, tex_coords);
}