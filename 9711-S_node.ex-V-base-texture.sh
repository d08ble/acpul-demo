### 9711:S node.ex/V/base-texture.sh
attribute vec4 a_position;
attribute vec2 a_texCoord;
varying mediump vec2 v_texCoord;

//attribute vec4 a_color;
//varying vec4 v_fragmentColor;

void main()
{
   gl_Position = CC_MVPMatrix * a_position+vec4(sin(a_position.x*0.05)*0.0,0.0,0.0,0.0);
//   gl_Position = CC_MVPMatrix * a_position+vec4(sin(a_position.x*0.05)*0.2,0.0,0.0,0.0);
//  gl_Position +=vec4(0.2);
//   v_texCoord = vec2(1.5)-a_texCoord*2.0;
  v_texCoord = a_texCoord;

   //v_fragmentColor = a_color;//vec4(1, 1, 1, 1);
}
