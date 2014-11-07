### 9710:S node.ex/V/base.sh
attribute vec4 a_position;

void main()
{
   gl_Position = CC_MVPMatrix * a_position;
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
