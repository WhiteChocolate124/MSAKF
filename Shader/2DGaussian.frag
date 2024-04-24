#version 330
#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D src;
uniform float variance = 0.02;

const int kernelSize = 3;

const float PI = 3.14159265358979323846;

in vec2 tex_coords;
out vec4 frag_color;

float gaussian(float x, float sigma) {
    return exp(-(x * x) / (2.0 * sigma * sigma)) / (sqrt(2.0 * PI) * sigma);
}

void main() {
    float sigma = sqrt(variance);

    vec2 resolution = textureSize(src, 0);

    float kernel[kernelSize];
    float totalWeight = 0.0;
    for (int i = 0; i < kernelSize; ++i) {
        float offset = float(i - (kernelSize - 1) / 2);
        kernel[i] = gaussian(offset, sigma);
        totalWeight += kernel[i];
    }

    for (int i = 0; i < kernelSize; ++i) {
        kernel[i] /= totalWeight;
    }

    vec3 sum = vec3(0.0);
    for (int i = 0; i < kernelSize; ++i) {
        for (int j = 0; j < kernelSize; ++j) {
            vec2 offset = vec2(float(i - (kernelSize - 1) / 2), float(j - (kernelSize - 1) / 2));
            vec2 texCoordOffset = tex_coords + offset / resolution;
            sum += texture2D(src, texCoordOffset).rgb * kernel[i] * kernel[j];
        }
    }
    frag_color = vec4(sum, 1.0);
}