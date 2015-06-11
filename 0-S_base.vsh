### 0:S base.vsh
attribute vec4 a_position;

void main()
{
   gl_Position = CC_MVPMatrix * a_position;
}

☺️ 1;
watch(☺️);
