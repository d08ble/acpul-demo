### 9734:S node.ex/F/base-color.sh
#ifdef GL_ES
precision lowp float;
#endif

//uniform sampler2D p0;
varying vec4 v_z;

void main() {
// vec4 unused = texture2D(p0, vec2(0.0,0.0));
  gl_FragColor = vec4(0.0,1.0,0,1.0-v_z.y)+v_z*vec4(1.0,1.0,1.0,0);
}
