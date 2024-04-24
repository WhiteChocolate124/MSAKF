#version 330
#ifdef GL_ES
precision mediump float;
#endif

//StructTensor
uniform sampler2D src;

in vec2 tex_coords;
out vec4 frag_color;

float normpdf(float x, float sig){
	return 0.39894*exp(-0.5*x*x/(sig*sig))/sig;
}

vec4 Variance(sampler2D tex, vec2 Coord){
  vec2 resolution = textureSize(tex, 0);

  vec3 sum = vec3(0.);
  const int mSize = 10;
  const int kSize = (mSize-1)/2;
  float kernel[mSize];
  
  float sig = 5.0;
  float Z = 0.0;
  for (int j = 0; j <= kSize; ++j){
    kernel[kSize+j] = kernel[kSize-j] = normpdf(float(j), sig);
  }
  
  for (int j = 0; j < mSize; ++j){
    Z += kernel[j];
  }
  
  for (int i=-kSize; i <= kSize; ++i){
    for (int j=-kSize; j <= kSize; ++j){
      sum += kernel[kSize+j]*kernel[kSize+i]*
        texture2D(tex, (Coord.xy+vec2(float(i),float(j))) 
                / resolution).rgb;
	}
  }
  
  sum = sum/(Z*Z);
  
  float lambda1 = 0.5 * (sum.y + sum.x +
                         sqrt(sum.y*sum.y - 2.0*sum.x*sum.y + sum.x*sum.x
                              + 4.0*sum.z*sum.z));
  float lambda2 = 0.5 * (sum.y + sum.x -
                         sqrt(sum.y*sum.y - 2.0*sum.x*sum.y + sum.x*sum.x 
                              + 4.0*sum.z*sum.z));
  
  vec2 v = vec2(lambda1 - sum.x, -sum.z);
  vec2 t;
  if (length(v) > 0.0) { 
    t = normalize(v);
  } else {
    t = vec2(0.0, 1.0);
  }
  
  float phi = atan(t.y, t.x);

  float A = (lambda1 + lambda2 > 0.0)?(lambda1 - lambda2) / (lambda1 + lambda2) : 0.0;
  
  return vec4(t, phi, A);
}

void main(){
  frag_color = Variance(src, tex_coords);
}