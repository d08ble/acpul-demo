### 9732:S node.ex/F/texXYDXDY_Move.sh
#ifdef GL_ES
precision lowp float;
#endif

varying vec2 v_texCoord;

uniform sampler2D p0;

void main() {
// vec4 a = vec4(v_texCoord.x, v_texCoord.y, 0.5, 1);
 vec4 b = texture2D(p0, v_texCoord);
// a.x += a.z*0.1;
// a.y += a.a*0.1;
// a = vec4(0);
 b.x += 0.01;
 b.y += 0.01;
 b.z -= 0.03;
// b.x = b.y = b.z;
// a.r=a.g=a.b = a.z;
 if (b.z < 0.01)
  b=vec4(0);
 gl_FragColor = b;

}
