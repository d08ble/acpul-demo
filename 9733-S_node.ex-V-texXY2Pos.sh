### 9733:S node.ex/V/texXY2Pos.sh
//attribute vec4 a_position;
attribute vec4 a0;

//varying vec2 v_texCoord;
//uniform sampler2D p0;

varying vec4 v_z;

void main()
{
v_z = a0;
// vec4 b = 
//texture2D(p0, a0);
// vec4 p = vec4(b.x*10.0+a0.x, b.y*10.0+a0.y, 0.0,1.0);
 vec4 p = vec4(a0.x*500.0+70.0, a0.y*500.0+50.0, 0.0, 1.0);
// vec4 b = a0;
 p.x+=(1.0+sin(a0.y*20.0))*100.0;
 p.x=320.0+(p.x-400.0)*sin(a0.y/2.0);
 p.y+=sin(a0.x*1009.0)*2.0;
// p = vec4(500.0,500.0, 0.0, 1.0)+a0;
 gl_Position = CC_MVPMatrix * p;//p;//a_position;
 gl_PointSize = 1.0;//(1.0-a0.x)*5.0;
}

//bilt-in:
//CC_PMatrix
//CC_MVMatrix
//+CC_MVPMatrix
//CC_Time
//CC_SinTime
//CC_CosTime
//CC_Random01
//CC_Texture0
//CC_alpha_value
//a_color
//+a_position
//a_texCoord
