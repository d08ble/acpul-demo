### 9712:S node.ex/V/base-attribute.sh
attribute vec4 a0;

void main()
{
   gl_Position = CC_MVPMatrix * a0;
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
