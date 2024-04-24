#version 330
#ifdef GL_ES
precision mediump float;
#endif

in vec2 tex_coords;
out vec4 frag_color;

//Image
uniform sampler2D src;
uniform sampler2D prev;

//Variance//Anisotropy
uniform sampler2D srcVar;

uniform float radius;
uniform float q = 12.0;
uniform float alpha;

uniform float KSize;

const float PI = 3.14159265358979323846;

uniform float ps = 0.5;
uniform float pd = 1.25;
uniform float tauv = 0.1;

float LocalSD(sampler2D Color, sampler2D Variance, vec2 uv){
    vec2 resolution = textureSize(Color, 0);

    vec3 s[4];
    for (int k = 0; k < 4; ++k) {
        s[k] = vec3(0.0);
    }

    float piN = 2.0 * PI / float(4);
    mat2 X = mat2(cos(piN), sin(piN), -sin(piN), cos(piN));

    vec4 t = texture2D(Variance, uv);

    float a = radius * clamp((alpha + t.w) / alpha, 0.1, 2.0);
    float b = radius * clamp(alpha / (alpha + t.w), 0.1, 2.0);

    float cos_phi = cos(t.z);
    float sin_phi = sin(t.z);

    mat2 SR = mat2(cos_phi/a, -sin_phi/b, sin_phi/a, cos_phi/b);

    int max_x = int(sqrt(a*a * cos_phi*cos_phi +
                          b*b * sin_phi*sin_phi));
    int max_y = int(sqrt(a*a * sin_phi*sin_phi +
                          b*b * cos_phi*cos_phi));

    for (int j = 0; j <= max_y; ++j)  {
        for (int i = -max_x; i <= max_x; ++i) {
            if ((j !=0) || (i > 0)) {
                vec2 v = SR * vec2(float(i),float(j));
                float dot_v = dot(v,v);
                if (dot_v <= KSize){
                    vec4 c0_fix = texture2D(Color, uv + vec2(i,j) / resolution);
                    vec3 c0 = c0_fix.rgb;
                    vec4 c1_fix = texture2D(Color, uv - vec2(i,j) / resolution);
                    vec3 c1 = c1_fix.rgb;

                    vec3 cc0 = c0 * c0;
                    vec3 cc1 = c1 * c1;

                    float n = 0.0;
                    float wx[4];
                    {
                        float z;
                        float xx = 0.33 - 0.84 * v.x * v.x;
                        float yy = 0.33 - 0.84 * v.y * v.y;

                        z = max(0.0, v.y + xx);
                        n += wx[0] = z * z;

                        z = max(0.0, -v.x + yy);
                        n += wx[1] = z * z;

                        z = max(0.0, -v.y + xx);
                        n += wx[2] = z * z;

                        z = max(0.0, v.x + yy);
                        n += wx[3] = z * z;
                    }

                    float g = exp(-3.125 * dot_v) / n;

                    for (int k = 0; k < 4; ++k) {
                        float w = wx[k] * g;
                        s[k] += cc0 * w;
                        s[(k+2)&3] += cc1 * w;
                    }
                }
            }
        }
    }
    float sumS = 0.0;
    for (int k = 0; k < 4; ++k){
        sumS += pow(s[k].r, 2) + pow(s[k].g, 2) + pow(s[k].b, 2);
    }

    return sqrt(sumS / 12.0);
}

void main(){
    vec4 current = texture2D(src, tex_coords);
    vec4 previous = texture2D(prev, tex_coords);

    float scalefactor = textureSize(src, 0).x / textureSize(prev, 0).x;

    float currS = LocalSD(src, srcVar, tex_coords);

    float Beta = clamp(currS * ps * pow(pd, scalefactor) + sqrt(distance(current, previous)) - tauv, 0, 1);
    vec4 merged = Beta * current + (1 - Beta) * previous;
    frag_color = merged;
}