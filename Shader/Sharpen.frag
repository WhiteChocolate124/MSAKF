#version 330
#ifdef GL_ES
precision mediump float;
#endif

uniform sampler2D src;

in vec2 tex_coords;
out vec4 frag_color;

uniform float SHARPEN_FACTOR = 16.0;

vec4 sharpenMask (sampler2D tex, vec2 fragCoord)
{
    vec2 resolution = textureSize(tex, 0);
    // Sharpen detection matrix [0,1,0],[1,-4,1],[0,1,0]
    // Colors
    vec4 up = texture2D(tex, (fragCoord + vec2 (0, 1))/resolution);
    vec4 left = texture2D(tex, (fragCoord + vec2 (-1, 0))/resolution);
    vec4 center = texture2D(tex, fragCoord/resolution);
    vec4 right = texture2D(tex, (fragCoord + vec2 (1, 0))/resolution);
    vec4 down = texture2D(tex, (fragCoord + vec2 (0, -1))/resolution);

    // Return edge detection
    return (1.0 + 4.0*SHARPEN_FACTOR)*center -SHARPEN_FACTOR*(up + left + right + down);
}

void main(){
    frag_color = sharpenMask(src, tex_coords);
}