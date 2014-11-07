### 100:S test.fsh
#ifdef GL_ES
precision lowp float;
#endif

varying vec2 v_texCoord;

uniform sampler2D CC_Texture0;

void main() {
//  vec2 a = CC_Random01.xy*100.0+10000.0;
  gl_FragColor = texture2D(CC_Texture0, v_texCoord)*vec4(0.5,0.3,1.0,0.5);
}
