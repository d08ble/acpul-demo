### 10:S base-texture.vsh
attribute vec4 a_position;
attribute vec2 a_texCoord;
//attribute vec4 a_color;

varying mediump vec2 v_texCoord;
//varying vec4 v_fragmentColor;

void main()
{
   gl_Position = CC_MVPMatrix * a_position;
   //v_fragmentColor = a_color;//vec4(1, 1, 1, 1);
   v_texCoord = a_texCoord;
}
