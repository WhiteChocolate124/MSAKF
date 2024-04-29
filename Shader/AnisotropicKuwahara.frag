#version 330
#ifdef GL_ES
precision mediump float;
#endif

in vec2 tex_coords;
out vec4 frag_color;

//Kuwahara
uniform sampler2D src;
//Anisotropy
uniform sampler2D Variance;

uniform float sigma;

struct lic_t { 
    vec2 p; 
    vec2 t;
    float w;
    float dw;
};

void step2(inout lic_t s) {
    vec2 resolution = textureSize(src, 0);
    vec2 t = texture(src, s.p).xy;
    if (dot(t, s.t) < 0.0) t = -t;
    s.t = t;

    s.dw = (abs(t.x) > abs(t.y))?
        abs((fract(s.p.x) - 0.5 - sign(t.x)) / t.x) :
        abs((fract(s.p.y) - 0.5 - sign(t.y)) / t.y);

    s.p += t * s.dw / resolution;
    s.w += s.dw;
}

vec4 mainImage(sampler2D Color, sampler2D Varia, vec2 uv) {
    vec2 resolution = textureSize(Color, 0);

    float twoSigma2 = 2.0 * sigma * sigma;
    float halfWidth = 2.0 * sigma;

    vec3 c = texture( Color, uv ).xyz;
    float w = 1.0;

    lic_t a, b;
    a.p = b.p = uv;
    a.t = texture( Varia, uv ).xy / resolution;
    b.t = -a.t;
    a.w = b.w = 0.0;

    while (a.w < halfWidth) {
        step2(a);
        float k = a.dw * exp(-a.w * a.w / twoSigma2);
        c += k * texture(Color, a.p).xyz;
        w += k;
    }
    while (b.w < halfWidth) {
        step2(b);
        float k = b.dw * exp(-b.w * b.w / twoSigma2);
        c += k * texture(Color, b.p).xyz;
        w += k;
    }
    c /= w;

    return vec4(c, 1.0);// * 1.05;
}

void main(){
  frag_color = mainImage(src, Variance, tex_coords);
}