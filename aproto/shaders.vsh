
### 37000 shader-test

r0:=k==0;

draw(r0, 37001);
x:=-100;
y:=200;
k:=1;

### 37001 shader-test.draw

#color(u0, 0.1,0.1,0.1,0.1);
#rect(u0, 6, 0,0, 200,200);

_ @1010;
_ @sys.node;
_ @sys.display;     

#node.ex.shader.use(u0, 37100, 37101);
#node.ex.gl.rect(u0, 0,0, 100,100);

#node.ex.shader.use(u0, 37100, 37012);
#node.ex.gl.rect(u0, 0,0, 200,200);

#node.ex.shader.use(u0, 37100, 37016);
#node.ex.gl.rect(u0, 0,0, 200,200);

#node.ex.shader.use(u0, 37100, 37018);
#node.ex.gl.rect(u0, 0,0, 200,200);

#node.ex.shader.use(u0, 37100, 37020);
#node.ex.gl.rect(u0, 0,0, 200,200);
#node.ex.shader.use(u0, 37100, 37022);
#node.ex.gl.rect(u0, 0,0, 200,200);
#node.ex.shader.use(u0, 37100, 37024);
#node.ex.gl.rect(u0, 0,0, 100,100);

#node.ex.shader.use(u0, 37100, 37026);
#node.ex.gl.rect(u0, 0,0, 200,200);

#node.ex.shader.use(u0, 37100, 37028);
#node.ex.gl.rect(u0, 0,0, 200,200);

#node.ex.shader.use(u0, 37100, 37030);
#node.ex.gl.rect(u0, 0,0, 200,200);

#node.ex.shader.use(u0, 37100, 37032);
#node.ex.gl.rect(u0, 0,0, 200,200);

#node.ex.shader.use(u0, 37100, 37036); # gpu benchmark
#node.ex.gl.rect(u0, 0,0, 500,500);

#node.ex.shader.use(u0, 37100, 37034);
#node.ex.gl.rect(u0, 0,0, 500,500);

#node.ex.shader.use(u0, 37100, 37038);
#node.ex.gl.rect(u0, 0,0, 500,500);

#node.ex.shader.use(u0, 37100, 37040);
#node.ex.gl.rect(u0, 0,0, 200,200);

#node.ex.shader.use(u0, 37100, 37046);
#node.ex.gl.rect(u0, 0,0, 200,200);

b display.gl.blend;                         #   b is display.gl.blend
#b.func(u0,                                  #   
#        b.mode.src.alpha,                   #   src = SRC.ALPHA
#        b.mode.one.minus.src.alpha);        #   dst = 1-SRC.ALPHA

#b.func(u0, 6, 0);
#b.func(u0, 6, 2);


node.ex.shader.use(u0, 37100, 37048);
#node.ex.gl.rect(u0, 0,0, 200,200);
#node.ex.gl.rect(u0, 200,0, 200,200);

# 1. SHADER.NOISE --> TEXTURE.NOISE
# 2. TEXTURE.NOISE --> TEXTURE.[A,B...]
# 3. TEXTURE.A --> NORMALMAP.A
# 4. TEXTURE.A + NORMALMAP.A --> OUT

drawtex { _ node.ex;
 shader.use(u0, 37100, 37050);
 shader.uniform.texture(u0, 0, _0, 0);
 gl.rect(u0, _1,_2, _3,_4);
};

TEX.A l10;
TEX.B l11;
TEX.SIZE {w 200; h 200;};
step0 {
 if (r0)                                     
 {
  TEX.A:=node.ex.texture.create(r0, TEX.SIZE.w,TEX.SIZE.h);
  TEX.B:=node.ex.texture.create(r0, TEX.SIZE.w,TEX.SIZE.h);
 };
};

step1 { _ node.ex;
 gl.fbo(u0, TEX.A);                  # - FBO-->TEXTURE
 gl.clear(r0, 0,0,0,0);              # - TEXURE FILL COLOR

 shader.use(u0, 37100, 37048);
 gl.rect(u0, 0,0, TEX.SIZE.w,TEX.SIZE.h);
 
 gl.fbo(u0, -1);                     # - FBO-->DEFAULT
};

step3 { _ node.ex;
 gl.fbo(u0, TEX.B);
 gl.clear(r0, 0,0,0,0);

 shader.use(u0, 37100, 37052);
 shader.uniform.texture(u0, 0, TEX.A, 0);
 gl.rect(u0, 0,0, TEX.SIZE.w,TEX.SIZE.h);
 
 gl.fbo(u0, -1);
};

step4 { _ node.ex;
# gl.fbo(u0, TEX.B);
# gl.clear(r0, 0,0,0,0);

 shader.use(u0, 37100, 37054);
 shader.uniform.texture(u0, 0, TEX.A, 0);
 shader.uniform.texture(u0, 1, TEX.A, 1);
# shader.uniform.texture(u0, 0, TEX.B, 0);
 shader.uniform.texture(u0, 1, TEX.B, 1);

# shader.uniform.texture(u0, 1, TEX.B, 1);
 gl.rect(u0, 0,0, TEX.SIZE.w,TEX.SIZE.h);
 
# gl.fbo(u0, -1);
};

step0;
step1;
drawtex(TEX.A, 0,0, TEX.SIZE.w/2.,TEX.SIZE.h/2);
step3;
drawtex(TEX.B, 100,0, TEX.SIZE.w/2.,TEX.SIZE.h/2);
drawtex(TEX.B, 200,0, TEX.SIZE.w,TEX.SIZE.h);
step4;

### 37002 ---

// Texture+Normal [
### 37054:S Texture+Normal.fsh

uniform sampler2D p0;
uniform sampler2D p1;
varying mediump vec2 v_texCoord;

void main(void)
{
    vec3 normal = normalize(texture2D(p1, v_texCoord).rgb * 2.0 - 1.0); 
    vec3 light_pos = normalize(vec3(1.0, 1.0, sin(CC_Time[0])*0.5+1.5));
    float diffuse = max(dot(normal, light_pos), 0.0);
//    diffuse += 0.1;
//    float diffuse = 1.;
  
//    vec3 color = diffuse * texture2D(p0, v_texCoord).rgb *0. + 1.*texture2D(p1, v_texCoord).rgb;
    vec3 color = diffuse * texture2D(p0, v_texCoord).rgb;

    //-texture2D(p1, v_texCoord);
    gl_FragColor = vec4(color, 1.0); 
}

### 37055:S ---
// Texture+Normal ]


// NormalMap [
### 37052:S NormalMap.fsh

uniform sampler2D p0;
varying mediump vec2 v_texCoord;

vec4 texture_offset(sampler2D t, vec2 p, ivec2 o) {
    return texture2D(t, vec2(ivec2(p) + o));
}

void main(void)
{
    vec4 n = texture2D(p0, v_texCoord);
//    v=vec4(vec3(v.r+v.g+v.b), 1.);
//    v=vec4(vec3(v.r+v.g+v.b), 1.);
//    n = vec4(n.x*.5+.5, n.y*.5+.5, n.z*.5+.5, 1.);
/*
    const vec2 size = vec2(2.0,0.0);
    const ivec3 off = ivec3(-1,0,1);

    vec4 wave = n;
    float s11 = wave.x;

    float s01 = texture_offset(p0, v_texCoord, off.xy).x;
    float s21 = texture_offset(p0, v_texCoord, off.zy).x;
    float s10 = texture_offset(p0, v_texCoord, off.yx).x;
    float s12 = texture_offset(p0, v_texCoord, off.yz).x;
    vec3 va = normalize(vec3(size.xy,s21-s01));
    vec3 vb = normalize(vec3(size.yx,s12-s10));
    vec4 bump = vec4( cross(va,vb), s11 );

    gl_FragColor = vec4(bump);
    */
//    vec2 size = {2.0,0.0};
//    vec3 off = {-1.0,0.0,1.0};
    const vec2 size = vec2(2.0,0.0);
    const vec3 off = vec3(-1.,0.,1.);
    #define image p0
    vec2 uv = v_texCoord;
    vec4 color = texture2D(image, uv.xy);

    vec2 nTex = vec2(200., 200.);
    vec2 offxy = vec2(off.x/nTex.x , off.y/nTex.y);
    vec2 offzy = vec2(off.z/nTex.x , off.y/nTex.y);
    vec2 offyx = vec2(off.y/nTex.x , off.x/nTex.y);
    vec2 offyz = vec2(off.y/nTex.x , off.z/nTex.y);

    float s11 = color.x;
    float s01 = texture2D(image, uv.xy + offxy).x;
    float s21 = texture2D(image, uv.xy + offzy).x;
    float s10 = texture2D(image, uv.xy + offyx).x;
    float s12 = texture2D(image, uv.xy + offyz).x;
    vec3 va = vec3(size.x, size.y, s21-s01);
    vec3 vb = vec3(size.y, size.x, s12-s10);
    va = normalize(va);
    vb = normalize(vb);
    vec4 bump = vec4(vec3(cross(va,vb) / 2. + 0.5), 1.0);
    gl_FragColor = bump;
}


### 37053:S ---
// NormalMap ]

// TextureRAW [
### 37050:S TextureRAW.fsh

uniform sampler2D p0;
varying mediump vec2 v_texCoord;

void main(void)
{
    gl_FragColor = texture2D(p0, v_texCoord);
}

### 37051:S ---
// TextureRAW ]

// Fire [
### 37107:S Fire.fsh

#define iGlobalTime CC_Time[3]*10.0
#define iResolution vec2(400.0)
#define iMouse vec2(0.5)

float noise(vec3 p) //Thx to Las^Mercury
{
    vec3 i = floor(p);
    vec4 a = dot(i, vec3(1., 57., 21.)) + vec4(0., 57., 21., 78.);
    vec3 f = cos((p-i)*acos(-1.))*(-.5)+.5;
    a = mix(sin(cos(a)*a),sin(cos(1.+a)*(1.+a)), f.x);
    a.xy = mix(a.xz, a.yw, f.y);
    return mix(a.x, a.y, f.z);
}

float sphere(vec3 p, vec4 spr)
{
    return length(spr.xyz-p) - spr.w;
}

float flame(vec3 p)
{
    float d = sphere(p*vec3(1.,.5,1.), vec4(.0,-1.,.0,1.));
    return d + (noise(p+vec3(.0,iGlobalTime*2.,.0)) + noise(p*3.)*.5)*.25*(p.y) ;
}

float scene(vec3 p)
{
    return min(100.-length(p) , abs(flame(p)) );
}

vec4 raymarch(vec3 org, vec3 dir)
{
    float d = 0.0, glow = 0.0, eps = 0.02;
    vec3  p = org;
    bool glowed = false;
    
    for(int i=0; i<64; i++)
    {
        d = scene(p) + eps;
        p += d * dir;
        if( d>eps )
        {
            if(flame(p) < .0)
                glowed=true;
            if(glowed)
                glow = float(i)/64.;
        }
    }
    return vec4(p,glow);
}

void main()
{
    vec2 v = -1.0 + 2.0 * gl_FragCoord.xy / iResolution.xy;
    v.x *= iResolution.x/iResolution.y;
    
    vec3 org = vec3(0., -2., 4.);
    vec3 dir = normalize(vec3(v.x*1.6, -v.y, -1.5));
    
    vec4 p = raymarch(org, dir);
    float glow = p.w;
    
    vec4 col = mix(vec4(1.,.5,.1,1.), vec4(0.1,.5,1.,1.), p.y*.02+.4);
    
    gl_FragColor = mix(vec4(0.), col, pow(glow*2.,4.));
    //gl_FragColor = mix(vec4(1.), mix(vec4(1.,.5,.1,1.),vec4(0.1,.5,1.,1.),p.y*.02+.4), pow(glow*2.,4.));

}

### 37108:S ---

// Fire ]

# Clouds [

// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

//#define FULL_PROCEDURAL


#ifdef FULL_PROCEDURAL

// hash based 3d value noise
float hash( float n )
{
    return fract(sin(n)*43758.5453);
}
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f*f*(3.0-2.0*f);
    float n = p.x + p.y*57.0 + 113.0*p.z;
    return mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                   mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
               mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                   mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
}
#else

// LUT based 3d value noise
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    
    vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
    vec2 rg = texture2D( iChannel0, (uv+ 0.5)/256.0, -100.0 ).yx;
    return mix( rg.x, rg.y, f.z );
}
#endif

vec4 map( in vec3 p )
{
    float d = 0.2 - p.y;

    vec3 q = p - vec3(1.0,0.1,0.0)*iGlobalTime;
    float f;
    f  = 0.5000*noise( q ); q = q*2.02;
    f += 0.2500*noise( q ); q = q*2.03;
    f += 0.1250*noise( q ); q = q*2.01;
    f += 0.0625*noise( q );

    d += 3.0 * f;

    d = clamp( d, 0.0, 1.0 );
    
    vec4 res = vec4( d );

    res.xyz = mix( 1.15*vec3(1.0,0.95,0.8), vec3(0.7,0.7,0.7), res.x );
    
    return res;
}


vec3 sundir = vec3(-1.0,0.0,0.0);


vec4 raymarch( in vec3 ro, in vec3 rd )
{
    vec4 sum = vec4(0, 0, 0, 0);

    float t = 0.0;
    for(int i=0; i<64; i++)
    {
        if( sum.a > 0.99 ) continue;

        vec3 pos = ro + t*rd;
        vec4 col = map( pos );
        
        #if 1
        float dif =  clamp((col.w - map(pos+0.3*sundir).w)/0.6, 0.0, 1.0 );

        vec3 lin = vec3(0.65,0.68,0.7)*1.35 + 0.45*vec3(0.7, 0.5, 0.3)*dif;
        col.xyz *= lin;
        #endif
        
        col.a *= 0.35;
        col.rgb *= col.a;

        sum = sum + col*(1.0 - sum.a);  

        #if 0
        t += 0.1;
        #else
        t += max(0.1,0.025*t);
        #endif
    }

    sum.xyz /= (0.001+sum.w);

    return clamp( sum, 0.0, 1.0 );
}

void main(void)
{
    vec2 q = gl_FragCoord.xy / iResolution.xy;
    vec2 p = -1.0 + 2.0*q;
    p.x *= iResolution.x/ iResolution.y;
    vec2 mo = -1.0 + 2.0*iMouse.xy / iResolution.xy;
    
    // camera
    vec3 ro = 4.0*normalize(vec3(cos(2.75-3.0*mo.x), 0.7+(mo.y+1.0), sin(2.75-3.0*mo.x)));
    vec3 ta = vec3(0.0, 1.0, 0.0);
    vec3 ww = normalize( ta - ro);
    vec3 uu = normalize(cross( vec3(0.0,1.0,0.0), ww ));
    vec3 vv = normalize(cross(ww,uu));
    vec3 rd = normalize( p.x*uu + p.y*vv + 1.5*ww );

    
    vec4 res = raymarch( ro, rd );

    float sun = clamp( dot(sundir,rd), 0.0, 1.0 );
    vec3 col = vec3(0.6,0.71,0.75) - rd.y*0.2*vec3(1.0,0.5,1.0) + 0.15*0.5;
    col += 0.2*vec3(1.0,.6,0.1)*pow( sun, 8.0 );
    col *= 0.95;
    col = mix( col, res.xyz, res.w );
    col += 0.1*vec3(1.0,0.4,0.2)*pow( sun, 3.0 );
        
    gl_FragColor = vec4( col, 1.0 );
}

# Clouds ]

# Volcanic [
// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

//#define STEREO 

float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    
    vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
    vec2 rg = texture2D( iChannel0, (uv+ 0.5)/256.0, -100.0 ).yx;
    return mix( rg.x, rg.y, f.z );
}

float noise( in vec2 x )
{
    vec2 p = floor(x);
    vec2 f = fract(x);
    vec2 uv = p.xy + f.xy*f.xy*(3.0-2.0*f.xy);
    return texture2D( iChannel0, (uv+118.4)/256.0, -100.0 ).x;
}

#ifdef STEREO
#define lodbias -5.0
#else
#define lodbias 0.0
#endif

vec4 texcube( sampler2D sam, in vec3 p, in vec3 n )
{
    vec4 x = texture2D( sam, p.yz, lodbias );
    vec4 y = texture2D( sam, p.zx, lodbias );
    vec4 z = texture2D( sam, p.xy, lodbias );
    return x*abs(n.x) + y*abs(n.y) + z*abs(n.z);
}

//=====================================================================

float lava( vec2 p )
{
    p += vec2(2.0,4.0);
    float f;
    f  = 0.5000*noise( p ); p = p*2.02;
    f += 0.2500*noise( p ); p = p*2.03;
    f += 0.1250*noise( p ); p = p*2.01;
    f += 0.0625*noise( p );
    return f;
}

const mat3 m = mat3( 0.00,  0.80,  0.60,
                    -0.80,  0.36, -0.48,
                    -0.60, -0.48,  0.64 );

float displacement( vec3 p )
{
    p += vec3(1.0,0.0,0.8);
    
    float f;
    f  = 0.5000*noise( p ); p = m*p*2.02;
    f += 0.2500*noise( p ); p = m*p*2.03;
    f += 0.1250*noise( p ); p = m*p*2.01;
    f += 0.0625*noise( p ); 
    
    float n = noise( p*3.5 );
    f += 0.03*n*n;
    
    return f;
}

float mapTerrain( in vec3 pos )
{
    return pos.y*0.1 + (displacement(pos*vec3(0.8,1.0,0.8)) - 0.4)*(1.0-smoothstep(1.0,3.0,pos.y));
}

float raymarchTerrain( in vec3 ro, in vec3 rd )
{
    float maxd = 30.0;
    float h = 1.0;
    float t = 0.1;
    for( int i=0; i<160; i++ )
    {
        if( h<(0.001*t)||t>maxd ) break;
        h = mapTerrain( ro+rd*t );
        t += h;
    }

    if( t>maxd ) t=-1.0;
    return t;
}

vec3 calcNormal( in vec3 pos, in float t )
{
    vec3 eps = vec3( max(0.02,0.001*t),0.0,0.0);
    return normalize( vec3(
           mapTerrain(pos+eps.xyy) - mapTerrain(pos-eps.xyy),
           mapTerrain(pos+eps.yxy) - mapTerrain(pos-eps.yxy),
           mapTerrain(pos+eps.yyx) - mapTerrain(pos-eps.yyx) ) );

}

vec3 lig = normalize( vec3(-0.3,0.4,0.7) );
    
vec4 mapClouds( in vec3 pos )
{
    vec3 q = pos*0.5 + vec3(0.0,-iGlobalTime,0.0);
    
    float d;
    d  = 0.5000*noise( q ); q = q*2.02;
    d += 0.2500*noise( q ); q = q*2.03;
    d += 0.1250*noise( q ); q = q*2.01;
    d += 0.0625*noise( q );
        
    d = d - 0.55;
    d *= smoothstep( 0.5, 0.55, lava(0.1*pos.xz)+0.01 );

    d = clamp( d, 0.0, 1.0 );
    
    vec4 res = vec4( d );

    res.xyz = mix( vec3(1.0,0.8,0.7), 0.2*vec3(0.4,0.4,0.4), res.x );
    res.xyz *= 0.25;
    res.xyz *= 0.5 + 0.5*smoothstep( -2.0, 1.0, pos.y );
    
    return res;
}

vec4 raymarchClouds( in vec3 ro, in vec3 rd, in vec3 bcol, float tmax )
{
    vec4 sum = vec4( 0.0 );

    float sun = pow( clamp( dot(rd,lig), 0.0, 1.0 ),6.0 );
    float t = 0.0;
    for( int i=0; i<60; i++ )
    {
        if( t>tmax || sum.w>0.95 ) break;//continue;
        vec3 pos = ro + t*rd;
        vec4 col = mapClouds( pos );
        
        col.xyz += vec3(1.0,0.7,0.4)*0.4*sun*(1.0-col.w);
        col.xyz = mix( col.xyz, bcol, 1.0-exp(-0.00006*t*t*t) );
        
        col.rgb *= col.a;

        sum = sum + col*(1.0 - sum.a);  

        t += max(0.1,0.05*t);
    }

    sum.xyz /= (0.001+sum.w);

    return clamp( sum, 0.0, 1.0 );
}

float softshadow( in vec3 ro, in vec3 rd, float mint, float k )
{
    float res = 1.0;
    float t = mint;
    for( int i=0; i<48; i++ )
    {
        float h = mapTerrain(ro + rd*t);
        h = max( h, 0.0 );
        res = min( res, k*h/t );
        t += clamp( h, 0.02, 0.5 );
        if( h<0.0001 ) break;
    }
    return clamp(res,0.0,1.0);
}

vec3 path( float time )
{
    return vec3( 16.0*cos(0.2+0.5*.1*time*1.5), 1.5, 16.0*sin(0.1+0.5*0.11*time*1.5) );
    
}

void main( void )
{
    #ifdef STEREO
    float eyeID = mod(gl_FragCoord.x + mod(gl_FragCoord.y,2.0),2.0);
    #endif

    vec2 q = gl_FragCoord.xy / iResolution.xy;
    vec2 p = -1.0 + 2.0*q;
    p.x *= iResolution.x / iResolution.y;
    
    
    // camera   
    float off = step( 0.001, iMouse.z )*6.0*iMouse.x/iResolution.x;
    float time = 2.7+iGlobalTime + off;
    vec3 ro = path( time+0.0 );
    vec3 ta = path( time+1.6 );
    //ta.y *= 0.3 + 0.25*cos(0.11*time);
    ta.y *= 0.35 + 0.25*sin(0.09*time);
    float roll = 0.3*sin(1.0+0.07*time);
    
    // camera tx
    vec3 cw = normalize(ta-ro);
    vec3 cp = vec3(sin(roll), cos(roll),0.0);
    vec3 cu = normalize(cross(cw,cp));
    vec3 cv = normalize(cross(cu,cw));
    
    float r2 = p.x*p.x*0.32 + p.y*p.y;
    p *= (7.0-sqrt(37.5-11.5*r2))/(r2+1.0);

    vec3 rd = normalize( p.x*cu + p.y*cv + 2.1*cw );

    #ifdef STEREO
    vec3 fo = ro + rd*7.0; // put focus plane behind Mike
    ro -= 0.2*cu*eyeID;    // eye separation
    rd = normalize(fo-ro);
    #endif

    // sky   
    vec3 col = vec3(0.32,0.36,0.4) - rd.y*0.4;
    float sun = clamp( dot(rd,lig), 0.0, 1.0 );
    col += vec3(1.0,0.8,0.4)*0.2*pow( sun, 6.0 );
    col *= 0.9;

    vec3 bcol = col;
    

    // terrain  
    float t = raymarchTerrain(ro, rd);
    if( t>0.0 )
    {
        vec3 pos = ro + t*rd;
        vec3 nor = calcNormal( pos, t );
        vec3 ref = reflect( rd, nor );

        vec3 bn = -1.0 + 2.0*texcube( iChannel0, 3.0*pos/4.0, nor ).xyz;
        nor = normalize( nor + 0.6*bn );
        
        float hh = 1.0 - smoothstep( -2.0, 1.0, pos.y );

        // lighting
        float sun = clamp( dot( nor, lig ), 0.0, 1.0 );
        float sha = 0.0; if( sun>0.01) sha=softshadow(pos,lig,0.01,32.0);
        float bac = clamp( dot( nor, normalize(lig*vec3(-1.0,0.0,-1.0)) ), 0.0, 1.0 );
        float sky = 0.5 + 0.5*nor.y;
        float lav = smoothstep( 0.5, 0.55, lava(0.1*pos.xz) )*hh*clamp(0.5-0.5*nor.y,0.0,1.0);
        float occ = pow( (1.0-displacement(pos*vec3(0.8,1.0,0.8)))*1.6-0.5, 2.0 );

        float amb = 1.0;

        col = vec3(0.8);

        vec3 lin = vec3(0.0);
        lin += sun*vec3(1.64,1.27,0.99)*pow(vec3(sha),vec3(1.0,1.2,1.5));
        lin += sky*vec3(0.16,0.20,0.28)*occ;
        lin += bac*vec3(0.40,0.28,0.20)*occ;
        lin += amb*vec3(0.18,0.15,0.15)*occ;
        lin += lav*vec3(3.00,0.61,0.00);


        // surface shading/material     
        col = texcube( iChannel1, 0.5*pos, nor ).xyz;
        col = col*(0.2+0.8*texcube( iChannel2, 4.0*vec3(2.0,8.0,2.0)*pos, nor ).x);
        vec3 verde = vec3(1.0,0.9,0.2);
        verde *= texture2D( iChannel2, pos.xz, lodbias ).xyz;
        col = mix( col, 0.8*verde, hh );
        
        float vv = smoothstep( 0.0, 0.8, nor.y )*smoothstep(0.0, 0.1, pos.y-0.8 );
        verde = vec3(0.2,0.45,0.1);
        verde *= texture2D( iChannel2, 30.0*pos.xz, lodbias ).xyz;
        verde += 0.2*texture2D( iChannel2, 1.0*pos.xz, lodbias ).xyz;
        vv *= smoothstep( 0.0, 0.5, texture2D( iChannel2, 0.1*pos.xz + 0.01*nor.x ).x );
        col = mix( col, verde*1.1, vv );
        
        // light/surface interaction        
        col = lin * col;
        
        // atmospheric
        col = mix( col, (1.0-0.7*hh)*bcol, 1.0-exp(-0.00006*t*t*t) );
    }

    // sun glow
    col += vec3(1.0,0.6,0.2)*0.2*pow( sun, 2.0 )*clamp( (rd.y+0.4)/(0.0+0.4),0.0,1.0);
    
    // smoke    
    {
    if( t<0.0 ) t=600.0;
    vec4 res = raymarchClouds( ro, rd, bcol, t );
    col = mix( col, res.xyz, res.w );
    }

    // gamma    
    col = pow( clamp( col, 0.0, 1.0 ), vec3(0.45) );

    // contrast, desat, tint and vignetting 
    col = col*0.3 + 0.7*col*col*(3.0-2.0*col);
    col = mix( col, vec3(col.x+col.y+col.z)*0.33, 0.2 );
    col *= 1.3*vec3(1.06,1.1,1.0);
    col *= 0.5 + 0.5*pow( 16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y), 0.1 );

    #ifdef STEREO   
    col *= vec3( eyeID, 1.0-eyeID, 1.0-eyeID ); 
    #endif
    
    gl_FragColor = vec4( col, 1.0 );
}
# Volcanic ]
# Elevated [

// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

//stereo thanks to Croqueteer
//#define STEREO 

// value noise, and its analytical derivatives
vec3 noised( in vec2 x )
{
    vec2 p = floor(x);
    vec2 f = fract(x);
    vec2 u = f*f*(3.0-2.0*f);
    float a = texture2D(iChannel0,(p+vec2(0.5,0.5))/256.0,-100.0).x;
    float b = texture2D(iChannel0,(p+vec2(1.5,0.5))/256.0,-100.0).x;
    float c = texture2D(iChannel0,(p+vec2(0.5,1.5))/256.0,-100.0).x;
    float d = texture2D(iChannel0,(p+vec2(1.5,1.5))/256.0,-100.0).x;
    return vec3(a+(b-a)*u.x+(c-a)*u.y+(a-b-c+d)*u.x*u.y,
                6.0*f*(1.0-f)*(vec2(b-a,c-a)+(a-b-c+d)*u.yx));
}

const mat2 m2 = mat2(0.8,-0.6,0.6,0.8);

float terrain( in vec2 x )
{
    vec2  p = x*0.003;
    float a = 0.0;
    float b = 1.0;
    vec2  d = vec2(0.0);
    for( int i=0; i<6; i++ )
    {
        vec3 n = noised(p);
        d += n.yz;
        a += b*n.x/(1.0+dot(d,d));
        b *= 0.5;
        p = m2*p*2.0;
    }

    return 140.0*a;
}

float terrain2( in vec2 x )
{
    vec2  p = x*0.003;
    float a = 0.0;
    float b = 1.0;
    vec2  d = vec2(0.0);
    for( int i=0; i<14; i++ )
    {
        vec3 n = noised(p);
        d += n.yz;
        a += b*n.x/(1.0+dot(d,d));
        b *= 0.5;
        p = m2*p*2.0;
    }

    return 140.0*a;
}

float terrain3( in vec2 x )
{
    vec2  p = x*0.003;
    float a = 0.0;
    float b = 1.0;
    vec2  d = vec2(0.0);
    for( int i=0; i<4; i++ )
    {
        vec3 n = noised(p);
        d += n.yz;
        a += b*n.x/(1.0+dot(d,d));
        b *= 0.5;
        p = m2*p*2.0;
    }

    return 140.0*a;
}

float map( in vec3 p )
{
    return p.y - terrain(p.xz);
}

float interesct( in vec3 ro, in vec3 rd, in float tmin, in float tmax )
{
    float t = tmin;
    for( int i=0; i<120; i++ )
    {
        float h = map( ro + t*rd );
        if( h<(0.002*t) || t>tmax ) break;
        t += 0.5*h;
    }

    return t;
}

float softShadow(in vec3 ro, in vec3 rd )
{
    // real shadows 
    float res = 1.0;
    float t = 0.001;
    for( int i=0; i<48; i++ )
    {
        vec3  p = ro + t*rd;
        float h = map( p );
        res = min( res, 16.0*h/t );
        t += h;
        if( res<0.001 ||p.y>200.0 ) break;
    }
    return clamp( res, 0.0, 1.0 );
}

vec3 calcNormal( in vec3 pos, float t )
{
    vec2  eps = vec2( 0.002*t, 0.0 );
    return normalize( vec3( terrain2(pos.xz-eps.xy) - terrain2(pos.xz+eps.xy),
                            2.0*eps.x,
                            terrain2(pos.xz-eps.yx) - terrain2(pos.xz+eps.yx) ) );
}

vec3 camPath( float time )
{
    return 1100.0*vec3( cos(0.0+0.23*time), 0.0, cos(1.5+0.21*time) );
}
    
float fbm( vec2 p )
{
    float f = 0.0;
    f += 0.5000*texture2D( iChannel0, p/256.0 ).x; p = m2*p*2.02;
    f += 0.2500*texture2D( iChannel0, p/256.0 ).x; p = m2*p*2.03;
    f += 0.1250*texture2D( iChannel0, p/256.0 ).x; p = m2*p*2.01;
    f += 0.0625*texture2D( iChannel0, p/256.0 ).x;
    return f/0.9375;
}

void main( void )
{
    vec2 xy = -1.0 + 2.0*gl_FragCoord.xy/iResolution.xy;
    vec2 s = xy*vec2(iResolution.x/iResolution.y,1.0);

    #ifdef STEREO
    float isCyan = mod(gl_FragCoord.x + mod(gl_FragCoord.y,2.0),2.0);
    #endif
    
    float time = iGlobalTime*0.15 + 0.3 + 4.0*iMouse.x/iResolution.x;
    
    vec3 light1 = normalize( vec3(-0.8,0.4,-0.3) );

    // camera position
    vec3 ro = camPath( time );
    vec3 ta = camPath( time + 3.0 );
    ro.y = terrain3( ro.xz ) + 11.0;
    ta.y = ro.y - 20.0;
    float cr = 0.2*cos(0.1*time);

    // camera ray    
    vec3  cw = normalize(ta-ro);
    vec3  cp = vec3(sin(cr), cos(cr),0.0);
    vec3  cu = normalize( cross(cw,cp) );
    vec3  cv = normalize( cross(cu,cw) );
    vec3  rd = normalize( s.x*cu + s.y*cv + 2.0*cw );

    #ifdef STEREO
    ro += 2.0*cu*isCyan; // move camera to the right - the rd vector is still good
    #endif
    
    // bounding plane
    float tmin = 2.0;
    float tmax = 2000.0;
    float maxh = 210.0;
    float tp = (maxh-ro.y)/rd.y;
    if( tp>0.0 )
    {
        if( ro.y>maxh ) tmin = max( tmin, tp );
        else            tmax = min( tmax, tp );
    }

    float sundot = clamp(dot(rd,light1),0.0,1.0);
    vec3 col;
    float t = interesct( ro, rd, tmin, tmax );
    if( t>tmax)
    {
        // sky      
        col = vec3(0.3,.55,0.8)*(1.0-0.8*rd.y)*0.9;
        // sun
        col += 0.25*vec3(1.0,0.7,0.4)*pow( sundot,5.0 );
        col += 0.25*vec3(1.0,0.8,0.6)*pow( sundot,64.0 );
        col += 0.2*vec3(1.0,0.8,0.6)*pow( sundot,512.0 );
        // clouds
        vec2 sc = ro.xz + rd.xz*(1000.0-ro.y)/rd.y;
        col = mix( col, vec3(1.0,0.95,1.0), 0.5*smoothstep(0.5,0.8,fbm(0.0005*sc)) );
        // horizon
        col = mix( col, vec3(0.7,0.75,0.8), pow( 1.0-max(rd.y,0.0), 8.0 ) );
    }
    else
    {
        // mountains        
        vec3 pos = ro + t*rd;

        vec3 nor = calcNormal( pos, t );
        
        // rock
        float r = texture2D( iChannel0, 7.0*pos.xz/256.0 ).x;
        col = (r*0.25+0.75)*0.9*mix( vec3(0.08,0.05,0.03), vec3(0.10,0.09,0.08), texture2D(iChannel0,0.00007*vec2(pos.x,pos.y*48.0)).x );
        col = mix( col, 0.20*vec3(0.45,.30,0.15)*(0.50+0.50*r),smoothstep(0.70,0.9,nor.y) );
        col = mix( col, 0.15*vec3(0.30,.30,0.10)*(0.25+0.75*r),smoothstep(0.95,1.0,nor.y) );

        // snow
        float h = smoothstep(55.0,80.0,pos.y + 25.0*fbm(0.01*pos.xz) );
        float e = smoothstep(1.0-0.5*h,1.0-0.1*h,nor.y);
        float o = 0.3 + 0.7*smoothstep(0.0,0.1,nor.x+h*h);
        float s = h*e*o;
        col = mix( col, 0.29*vec3(0.62,0.65,0.7), smoothstep( 0.1, 0.9, s ) );
        
         // lighting        
        float amb = clamp(0.5+0.5*nor.y,0.0,1.0);
        float dif = clamp( dot( light1, nor ), 0.0, 1.0 );
        float bac = clamp( 0.2 + 0.8*dot( normalize( vec3(-light1.x, 0.0, light1.z ) ), nor ), 0.0, 1.0 );
        float sh = 1.0; if( dif>=0.0001 ) sh = softShadow(pos+light1*20.0,light1);
        
        vec3 lin  = vec3(0.0);
        lin += dif*vec3(7.00,5.00,3.00)*vec3( sh, sh*sh*0.5+0.5*sh, sh*sh*0.8+0.2*sh );
        lin += amb*vec3(0.40,0.60,0.80)*1.2;
        lin += bac*vec3(0.40,0.50,0.60);
        col *= lin;

        // fog
        float fo = 1.0-exp(-0.000001*t*t );
        vec3 fco = 0.7*vec3(0.5,0.7,0.9) + 0.1*vec3(1.0,0.8,0.5)*pow( sundot, 4.0 );
        col = mix( col, fco, fo );

        // sun scatter
        col += 0.3*vec3(1.0,0.8,0.4)*pow( sundot, 8.0 )*(1.0-exp(-0.002*t));
    }

    // gamma
    col = pow(col,vec3(0.4545));

    // vignetting   
    col *= 0.5 + 0.5*pow( (xy.x+1.0)*(xy.y+1.0)*(xy.x-1.0)*(xy.y-1.0), 0.1 );
    
    #ifdef STEREO   
    col *= vec3( isCyan, 1.0-isCyan, 1.0-isCyan );  
    #endif
    
    gl_FragColor=vec4(col,1.0);
}

# Elevated ]
// Water [

### 37105:S ---

#define iGlobalTime CC_Time[3]*10.0
#define iResolution vec2(400.0)
#define iMouse vec2(0.5)

// "Seascape" by Alexander Alekseev aka TDM - 2014
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

const int NUM_STEPS = 8;
const float PI      = 3.1415;
const float EPSILON = 1e-3;
float EPSILON_NRM   = 0.1 / iResolution.x;

// sea
const int ITER_GEOMETRY = 3;
const int ITER_FRAGMENT = 5;
const float SEA_HEIGHT = 0.6;
const float SEA_CHOPPY = 4.0;
const float SEA_SPEED = 0.8;
const float SEA_FREQ = 0.16;
const vec3 SEA_BASE = vec3(0.1,0.19,0.22);
const vec3 SEA_WATER_COLOR = vec3(0.8,0.9,0.6);
float SEA_TIME = iGlobalTime * SEA_SPEED;
mat2 octave_m = mat2(1.6,1.2,-1.2,1.6);

// math
mat3 fromEuler(vec3 ang) {
    vec2 a1 = vec2(sin(ang.x),cos(ang.x));
    vec2 a2 = vec2(sin(ang.y),cos(ang.y));
    vec2 a3 = vec2(sin(ang.z),cos(ang.z));
    mat3 m;
    m[0] = vec3(a1.y*a3.y+a1.x*a2.x*a3.x,a1.y*a2.x*a3.x+a3.y*a1.x,-a2.y*a3.x);
    m[1] = vec3(-a2.y*a1.x,a1.y*a2.y,a2.x);
    m[2] = vec3(a3.y*a1.x*a2.x+a1.y*a3.x,a1.x*a3.x-a1.y*a3.y*a2.x,a2.y*a3.y);
    return m;
}
float hash( vec2 p ) {
    float h = dot(p,vec2(127.1,311.7)); 
    return fract(sin(h)*43758.5453123);
}
float noise( in vec2 p ) {
    vec2 i = floor( p );
    vec2 f = fract( p );    
    vec2 u = f*f*(3.0-2.0*f);
    return -1.0+2.0*mix( mix( hash( i + vec2(0.0,0.0) ), 
                     hash( i + vec2(1.0,0.0) ), u.x),
                mix( hash( i + vec2(0.0,1.0) ), 
                     hash( i + vec2(1.0,1.0) ), u.x), u.y);
}

// lighting
float diffuse(vec3 n,vec3 l,float p) {
    return pow(dot(n,l) * 0.4 + 0.6,p);
}
float specular(vec3 n,vec3 l,vec3 e,float s) {    
    float nrm = (s + 8.0) / (3.1415 * 8.0);
    return pow(max(dot(reflect(e,n),l),0.0),s) * nrm;
}

// sky
vec3 getSkyColor(vec3 e) {
    e.y = max(e.y,0.0);
    vec3 ret;
    ret.x = pow(1.0-e.y,2.0);
    ret.y = 1.0-e.y;
    ret.z = 0.6+(1.0-e.y)*0.4;
    return ret;
}

// sea
float sea_octave(vec2 uv, float choppy) {
    uv += noise(uv);        
    vec2 wv = 1.0-abs(sin(uv));
    vec2 swv = abs(cos(uv));    
    wv = mix(wv,swv,wv);
    return pow(1.0-pow(wv.x * wv.y,0.65),choppy);
}

float map(vec3 p) {
    float freq = SEA_FREQ;
    float amp = SEA_HEIGHT;
    float choppy = SEA_CHOPPY;
    vec2 uv = p.xz; uv.x *= 0.75;
    
    float d, h = 0.0;    
    for(int i = 0; i < ITER_GEOMETRY; i++) {        
        d = sea_octave((uv+SEA_TIME)*freq,choppy);
        d += sea_octave((uv-SEA_TIME)*freq,choppy);
        h += d * amp;        
        uv *= octave_m; freq *= 1.9; amp *= 0.22;
        choppy = mix(choppy,1.0,0.2);
    }
    return p.y - h;
}

float map_detailed(vec3 p) {
    float freq = SEA_FREQ;
    float amp = SEA_HEIGHT;
    float choppy = SEA_CHOPPY;
    vec2 uv = p.xz; uv.x *= 0.75;
    
    float d, h = 0.0;    
    for(int i = 0; i < ITER_FRAGMENT; i++) {        
        d = sea_octave((uv+SEA_TIME)*freq,choppy);
        d += sea_octave((uv-SEA_TIME)*freq,choppy);
        h += d * amp;        
        uv *= octave_m; freq *= 1.9; amp *= 0.22;
        choppy = mix(choppy,1.0,0.2);
    }
    return p.y - h;
}

vec3 getSeaColor(vec3 p, vec3 n, vec3 l, vec3 eye, vec3 dist) {  
    float fresnel = 1.0 - max(dot(n,-eye),0.0);
    fresnel = pow(fresnel,3.0) * 0.65;
        
    vec3 reflected = getSkyColor(reflect(eye,n));    
    vec3 refracted = SEA_BASE + diffuse(n,l,80.0) * SEA_WATER_COLOR * 0.12; 
    
    vec3 color = mix(refracted,reflected,fresnel);
    
    float atten = max(1.0 - dot(dist,dist) * 0.001, 0.0);
    color += SEA_WATER_COLOR * (p.y - SEA_HEIGHT) * 0.18 * atten;
    
    color += vec3(specular(n,l,eye,60.0));
    
    return color;
}

// tracing
vec3 getNormal(vec3 p, float eps) {
    vec3 n;
    n.y = map_detailed(p);    
    n.x = map_detailed(vec3(p.x+eps,p.y,p.z)) - n.y;
    n.z = map_detailed(vec3(p.x,p.y,p.z+eps)) - n.y;
    n.y = eps;
    return normalize(n);
}

float heightMapTracing(vec3 ori, vec3 dir, out vec3 p) {  
    float tm = 0.0;
    float tx = 1000.0;    
    float hx = map(ori + dir * tx);
    if(hx > 0.0) return tx;   
    float hm = map(ori + dir * tm);    
    float tmid = 0.0;
    for(int i = 0; i < NUM_STEPS; i++) {
        tmid = mix(tm,tx, hm/(hm-hx));                   
        p = ori + dir * tmid;                   
        float hmid = map(p);
        if(hmid < 0.0) {
            tx = tmid;
            hx = hmid;
        } else {
            tm = tmid;
            hm = hmid;
        }
    }
    return tmid;
}

// main
void main(void) {
    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    uv = uv * 2.0 - 1.0;
    uv.x *= iResolution.x / iResolution.y;    
    float time = iGlobalTime * 0.3 + iMouse.x*0.01;
        
    // ray
    vec3 ang = vec3(sin(time*3.0)*0.1,sin(time)*0.2+0.3,time);    
    vec3 ori = vec3(0.0,3.5,time*5.0);
    vec3 dir = normalize(vec3(uv.xy,-2.0)); dir.z += length(uv) * 0.15;
    dir = normalize(dir) * fromEuler(ang);
    
    // tracing
    vec3 p;
    heightMapTracing(ori,dir,p);
    vec3 dist = p - ori;
    vec3 n = getNormal(p, dot(dist,dist) * EPSILON_NRM);
    vec3 light = normalize(vec3(0.0,1.0,0.8)); 
             
    // color
    vec3 color = mix(
        getSkyColor(dir),
        getSeaColor(p,n,light,dir,dist),
        pow(smoothstep(0.0,-0.05,dir.y),0.3));
        
    // post
    gl_FragColor = vec4(pow(color,vec3(0.75)), 0.5);
}

### 37106:S ---

// Water ]

# Sun [
uniform vec3      iResolution;           // viewport resolution (in pixels)
uniform float     iGlobalTime;           // shader playback time (in seconds)
uniform float     iChannelTime[4];       // channel playback time (in seconds)
uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
uniform vec4      iDate;                 // (year, month, day, time in seconds)
uniform float     iSampleRate;           // sound sample rate (i.e., 44100)

// based on https://www.shadertoy.com/view/lsf3RH by
// trisomie21 (THANKS!)
// My apologies for the ugly code.

float snoise(vec3 uv, float res)    // by trisomie21
{
    const vec3 s = vec3(1e0, 1e2, 1e4);
    
    uv *= res;
    
    vec3 uv0 = floor(mod(uv, res))*s;
    vec3 uv1 = floor(mod(uv+vec3(1.), res))*s;
    
    vec3 f = fract(uv); f = f*f*(3.0-2.0*f);
    
    vec4 v = vec4(uv0.x+uv0.y+uv0.z, uv1.x+uv0.y+uv0.z,
                  uv0.x+uv1.y+uv0.z, uv1.x+uv1.y+uv0.z);
    
    vec4 r = fract(sin(v*1e-3)*1e5);
    float r0 = mix(mix(r.x, r.y, f.x), mix(r.z, r.w, f.x), f.y);
    
    r = fract(sin((v + uv1.z - uv0.z)*1e-3)*1e5);
    float r1 = mix(mix(r.x, r.y, f.x), mix(r.z, r.w, f.x), f.y);
    
    return mix(r0, r1, f.z)*2.-1.;
}

float freqs[4];

void main(void)
{
    freqs[0] = texture2D( iChannel1, vec2( 0.01, 0.25 ) ).x;
    freqs[1] = texture2D( iChannel1, vec2( 0.07, 0.25 ) ).x;
    freqs[2] = texture2D( iChannel1, vec2( 0.15, 0.25 ) ).x;
    freqs[3] = texture2D( iChannel1, vec2( 0.30, 0.25 ) ).x;

    float brightness    = freqs[1] * 0.25 + freqs[2] * 0.25;
    float radius        = 0.24 + brightness * 0.2;
    float invRadius     = 1.0/radius;
    
    vec3 orange         = vec3( 0.8, 0.65, 0.3 );
    vec3 orangeRed      = vec3( 0.8, 0.35, 0.1 );
    float time      = iGlobalTime * 0.1;
    float aspect    = iResolution.x/iResolution.y;
    vec2 uv         = gl_FragCoord.xy / iResolution.xy;
    vec2 p          = -0.5 + uv;
    p.x *= aspect;

    float fade      = pow( length( 2.0 * p ), 0.5 );
    float fVal1     = 1.0 - fade;
    float fVal2     = 1.0 - fade;
    
    float angle     = atan( p.x, p.y )/6.2832;
    float dist      = length(p);
    vec3 coord      = vec3( angle, dist, time * 0.1 );
    
    float newTime1  = abs( snoise( coord + vec3( 0.0, -time * ( 0.35 + brightness * 0.001 ), time * 0.015 ), 15.0 ) );
    float newTime2  = abs( snoise( coord + vec3( 0.0, -time * ( 0.15 + brightness * 0.001 ), time * 0.015 ), 45.0 ) );  
    for( int i=1; i<=7; i++ ){
        float power = pow( 2.0, float(i + 1) );
        fVal1 += ( 0.5 / power ) * snoise( coord + vec3( 0.0, -time, time * 0.2 ), ( power * ( 10.0 ) * ( newTime1 + 1.0 ) ) );
        fVal2 += ( 0.5 / power ) * snoise( coord + vec3( 0.0, -time, time * 0.2 ), ( power * ( 25.0 ) * ( newTime2 + 1.0 ) ) );
    }
    
    float corona        = pow( fVal1 * max( 1.1 - fade, 0.0 ), 2.0 ) * 50.0;
    corona              += pow( fVal2 * max( 1.1 - fade, 0.0 ), 2.0 ) * 50.0;
    corona              *= 1.2 - newTime1;
    vec3 sphereNormal   = vec3( 0.0, 0.0, 1.0 );
    vec3 dir            = vec3( 0.0 );
    vec3 center         = vec3( 0.5, 0.5, 1.0 );
    vec3 starSphere     = vec3( 0.0 );
    
    vec2 sp = -1.0 + 2.0 * uv;
    sp.x *= aspect;
    sp *= ( 2.0 - brightness );
    float r = dot(sp,sp);
    float f = (1.0-sqrt(abs(1.0-r)))/(r) + brightness * 0.5;
    if( dist < radius ){
        corona          *= pow( dist * invRadius, 24.0 );
        vec2 newUv;
        newUv.x = sp.x*f;
        newUv.y = sp.y*f;
        newUv += vec2( time, 0.0 );
        
        vec3 texSample  = texture2D( iChannel0, newUv ).rgb;
        float uOff      = ( texSample.g * brightness * 4.5 + time );
        vec2 starUV     = newUv + vec2( uOff, 0.0 );
        starSphere      = texture2D( iChannel0, starUV ).rgb;
    }
    
    float starGlow  = min( max( 1.0 - dist * ( 1.0 - brightness ), 0.0 ), 1.0 );
    //gl_FragColor.rgb  = vec3( r );
    gl_FragColor.rgb    = vec3( f * ( 0.75 + brightness * 0.3 ) * orange ) + starSphere + corona * orange + starGlow * orangeRed;
    gl_FragColor.a      = 1.0;
}

# Sun ]
# Fire [

uniform vec3      iResolution;           // viewport resolution (in pixels)
uniform float     iGlobalTime;           // shader playback time (in seconds)
uniform float     iChannelTime[4];       // channel playback time (in seconds)
uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
uniform vec4      iDate;                 // (year, month, day, time in seconds)
uniform float     iSampleRate;           // sound sample rate (i.e., 44100)

### 37103:S Fire.fsh

#define iGlobalTime CC_Time[3]
#define iResolution vec2(200.0) 

float snoise(vec3 uv, float res)
{
    const vec3 s = vec3(1e0, 1e2, 1e3);
    
    uv *= res;
    
    vec3 uv0 = floor(mod(uv, res))*s;
    vec3 uv1 = floor(mod(uv+vec3(1.), res))*s;
    
    vec3 f = fract(uv); f = f*f*(3.0-2.0*f);

    vec4 v = vec4(uv0.x+uv0.y+uv0.z, uv1.x+uv0.y+uv0.z,
                  uv0.x+uv1.y+uv0.z, uv1.x+uv1.y+uv0.z);

    vec4 r = fract(sin(v*1e-1)*1e3);
    float r0 = mix(mix(r.x, r.y, f.x), mix(r.z, r.w, f.x), f.y);
    
    r = fract(sin((v + uv1.z - uv0.z)*1e-1)*1e3);
    float r1 = mix(mix(r.x, r.y, f.x), mix(r.z, r.w, f.x), f.y);
    
    return mix(r0, r1, f.z)*2.-1.;
}

void main(void) 
{
    vec2 p = -.5 + gl_FragCoord.xy / iResolution.xy;
    p.x *= iResolution.x/iResolution.y;
    
    float color = 3.0 - (3.*length(2.*p));
    
    vec3 coord = vec3(atan(p.x,p.y)/6.2832+.5, length(p)*.4, .5);
    
    for(int i = 1; i <= 7; i++)
    {
        float power = pow(2.0, float(i));
        color += (1.5 / power) * snoise(coord + vec3(0.,-iGlobalTime*.05, iGlobalTime*.01), power*16.);
    }
    gl_FragColor = vec4( color, pow(max(color,0.),2.)*0.4, pow(max(color,0.),3.)*0.15 , 1.0);
}

### 37104:S ---

# Fire ]

# Cells [
uniform vec3      iResolution;           // viewport resolution (in pixels)
uniform float     iGlobalTime;           // shader playback time (in seconds)
uniform float     iChannelTime[4];       // channel playback time (in seconds)
uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
uniform vec4      iDate;                 // (year, month, day, time in seconds)
uniform float     iSampleRate;           // sound sample rate (i.e., 44100)

// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.


// I've not seen anybody out there computing correct cell interior distances for Voronoi
// patterns yet. That's why they cannot shade the cell interior correctly, and why you've
// never seen cell boundaries rendered correctly. 

// However, here's how you do mathematically correct distances (note the equidistant and non
// degenerated grey isolines inside the cells) and hence edges (in yellow):

// http://www.iquilezles.org/www/articles/voronoilines/voronoilines.htm

#define ANIMATE

vec2 hash2( vec2 p )
{
    // texture based white noise
    return texture2D( iChannel0, (p+0.5)/256.0, -100.0 ).xy;
    
    // procedural white noise   
    //return fract(sin(vec2(dot(p,vec2(127.1,311.7)),dot(p,vec2(269.5,183.3))))*43758.5453);
}

vec3 voronoi( in vec2 x )
{
    vec2 n = floor(x);
    vec2 f = fract(x);

    //----------------------------------
    // first pass: regular voronoi
    //----------------------------------
    vec2 mg, mr;

    float md = 8.0;
    for( int j=-1; j<=1; j++ )
    for( int i=-1; i<=1; i++ )
    {
        vec2 g = vec2(float(i),float(j));
        vec2 o = hash2( n + g );
        #ifdef ANIMATE
        o = 0.5 + 0.5*sin( iGlobalTime + 6.2831*o );
        #endif  
        vec2 r = g + o - f;
        float d = dot(r,r);

        if( d<md )
        {
            md = d;
            mr = r;
            mg = g;
        }
    }

    //----------------------------------
    // second pass: distance to borders
    //----------------------------------
    md = 8.0;
    for( int j=-2; j<=2; j++ )
    for( int i=-2; i<=2; i++ )
    {
        vec2 g = mg + vec2(float(i),float(j));
        vec2 o = hash2( n + g );
        #ifdef ANIMATE
        o = 0.5 + 0.5*sin( iGlobalTime + 6.2831*o );
        #endif  
        vec2 r = g + o - f;

        if( dot(mr-r,mr-r)>0.00001 )
        md = min( md, dot( 0.5*(mr+r), normalize(r-mr) ) );
    }

    return vec3( md, mr );
}

void main( void )
{
    vec2 p = gl_FragCoord.xy/iResolution.xx;

    vec3 c = voronoi( 8.0*p );

    // isolines
    vec3 col = c.x*(0.5 + 0.5*sin(64.0*c.x))*vec3(1.0);
    // borders  
    col = mix( vec3(1.0,0.6,0.0), col, smoothstep( 0.04, 0.07, c.x ) );
    // feature points
    float dd = length( c.yz );
    col = mix( vec3(1.0,0.6,0.1), col, smoothstep( 0.0, 0.12, dd) );
    col += vec3(1.0,0.6,0.1)*(1.0-smoothstep( 0.0, 0.04, dd));

    gl_FragColor = vec4(col,1.0);
}

# Cells ]

# Muto [

// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

vec3 hash3( float n )
{
    return fract(sin(vec3(n,n+1.0,n+2.0))*vec3(43758.5453123,22578.1459123,19642.3490423));
}

vec3 noise( in float x )
{
    float p = floor(x);
    float f = fract(x);
    f = f*f*(3.0-2.0*f);
    return mix( hash3(p+0.0), hash3(p+1.0),f);
}


mat4 rotationMat( in vec3 xyz )
{
    vec3 si = sin(xyz);
    vec3 co = cos(xyz);

    return mat4( co.y*co.z,                co.y*si.z,               -si.y,       0.0,
                 si.x*si.y*co.z-co.x*si.z, si.x*si.y*si.z+co.x*co.z, si.x*co.y,  0.0,
                 co.x*si.y*co.z+si.x*si.z, co.x*si.y*si.z-si.x*co.z, co.x*co.y,  0.0,
                 0.0,                      0.0,                      0.0,        1.0 );
}

const float s = 1.1;

mat4 mm;

vec3 map( vec3 p )
{
    float k = 1.0;
    float m = 1e10;
    for( int i=0; i<22; i++ ) 
    {
        m = min( m, dot(p,p)/(k*k) );
        p = (mm*vec4((abs(p)),1.0)).xyz;
        k*= s;
    }
    

    float d = (length(p)-0.25)/k;
    
    float h = p.z - 0.35*p.x;
    
    return vec3( d, m, h );
}

vec3 intersect( in vec3 ro, in vec3 rd )
{
    float t = 0.0;
    for( int i=0; i<100; i++ )
    {
        vec3 res = map( ro+rd*t );
        if( res.x<0.0002 ) return vec3(t,res.yz);
        t += res.x;
        if( t>9.0 ) break;
    }

    return vec3( -1.0 );
}

vec3 calcNormal( in vec3 pos, float e )
{
    vec3 eps = vec3(e,0.0,0.0);

    return normalize( vec3(
           map(pos+eps.xyy).x - map(pos-eps.xyy).x,
           map(pos+eps.yxy).x - map(pos-eps.yxy).x,
           map(pos+eps.yyx).x - map(pos-eps.yyx).x ) );
}

float softshadow( in vec3 ro, in vec3 rd, float mint, float k )
{
    float res = 1.0;
    float t = mint;
    for( int i=0; i<32; i++ )
    {
        float h = map(ro + rd*t).x;
        h = max( h, 0.0 );
        res = min( res, k*h/t );
        t += clamp( h, 0.001, 0.1 );
        if( res<0.01 || t>6.0 ) break;
    }
    return clamp(res,0.0,1.0);
}

float calcAO( in vec3 pos, in vec3 nor )
{
    float totao = 0.0;
    for( int aoi=0; aoi<16; aoi++ )
    {
        vec3 aopos = -1.0+2.0*hash3(float(aoi)*213.47);
        aopos *= sign( dot(aopos,nor) );
        aopos = pos + nor*0.01 + aopos*0.04;
        float dd = clamp( map( aopos ).x*4.0, 0.0, 1.0 );
        totao += dd;
    }
    totao /= 16.0;
    
    return clamp( totao*totao*50.0, 0.0, 1.0 );
}


void main(void)
{
    vec2 q = gl_FragCoord.xy / iResolution.xy;
    vec2 p = -1.0 + 2.0 * q;
    p.x *= iResolution.x/iResolution.y;
    vec2 m = vec2(0.5);
    if( iMouse.z>0.0 ) m = iMouse.xy/iResolution.xy;

    // animation    
    float time = iGlobalTime;
    time += 15.0*smoothstep(  15.0, 25.0, iGlobalTime );
    time += 20.0*smoothstep(  65.0, 80.0, iGlobalTime );
    time += 35.0*smoothstep( 105.0, 135.0, iGlobalTime );
    time += 20.0*smoothstep( 165.0, 180.0, iGlobalTime );
    time += 40.0*smoothstep( 220.0, 290.0, iGlobalTime );
    time +=  5.0*smoothstep( 320.0, 330.0, iGlobalTime );
    float time1 = (time-10.0)*1.5 - 167.0;
    float time2 = time;
    
    mm = rotationMat( vec3(0.4,0.1,3.4) + 
                      0.15*sin(0.1*vec3(0.40,0.30,0.61)*time1) + 
                      0.15*sin(0.1*vec3(0.11,0.53,0.48)*time1));
    mm[0].xyz *= s; 
    mm[1].xyz *= s;
    mm[2].xyz *= s; 
    mm[3].xyz = vec3( 0.15, 0.05, -0.07 ) + 0.05*sin(vec3(0.0,1.0,2.0) + 0.2*vec3(0.31,0.24,0.42)*time1);
    
    // camera
    float an = 1.0 + 0.1*time2 - 6.2*m.x;
    float cr = 0.15*sin(0.2*time2);
    vec3 ro = (2.4 + 0.6*smoothstep(10.0,20.0,time2))*vec3(sin(an),0.25,cos(an));
    vec3 ta = vec3( 0.0, 0.0 + 0.13*cos(0.3*time2), 0.0 );
    ta += 0.05*noise(  0.0 + 1.0*time );
    ro += 0.05*noise( 11.3 + 1.0*time );
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(sin(cr),cos(cr),0.0) ) );
    vec3 vv = normalize( cross(uu,ww));
    vec3 rd = normalize( p.x*uu + p.y*vv + 3.0*ww );

    // raymarch
    vec3 tmat = intersect(ro,rd);
    
    // shade
    vec3 col = vec3(0.0);
    if( tmat.z>-0.5 )
    {
        // geometry
        vec3 pos = ro + tmat.x*rd;
        vec3 nor = calcNormal(pos, 0.005);
        vec3 sor = calcNormal(pos, 0.010);

        // material
        vec3 mate = vec3(1.0);
        mate = mix( vec3(0.5,0.5,0.2), vec3(0.5,0.3,0.0), 0.5 + 0.5*sin(4.0+8000.0*tmat.y)  );
        mate = mix( vec3(1.0,0.9,0.8), mate, 0.5 + 0.5*sin(4.0+20.0*tmat.z) );
        mate.x *= 1.15;

        // lighting
        float occ = 1.1*calcAO( pos, nor );
        occ *= 0.75 + 0.25*clamp(tmat.y*400.0,0.0,1.0);
        
        // diffuse
        col = vec3(0.0);
        for( int i=0; i<32; i++ )
        {
            //vec3 rr = normalize(-1.0 + 2.0*texture2D( iChannel2, vec2((0.5+float(i)),0.5)/256.0,-100.0).xyz);
            vec3 rr = normalize(-1.0 + 2.0*hash3(float(i)*123.5463));
            rr = normalize( nor + 7.0*rr );
            rr = rr * sign(dot(nor,rr));                              
            float ds = occ;//softshadow( pos, rr, 0.01, 32.0 );
            col += pow( textureCube( iChannel0, rr ).xyz, vec3(2.2) ) * dot(rr,nor) * ds;
        }
        col /= 32.0;                                        

        col *= 1.8;

        // subsurface       
        col *= 1.0 + 1.0*vec3(1.0,0.6,0.1)*pow(clamp(1.0+dot(rd,sor),0.0,1.0),2.0)*vec3(1.0);
        
        // specular     
        float fre = pow( clamp(1.0+dot(rd,nor),0.0,1.0), 5.0 );
        vec3 ref = reflect( rd, nor );
        float rs = softshadow( pos, ref, 0.01, 32.0 );
        col += 1.8 * (0.04 + 12.0*fre) * occ * pow( textureCube( iChannel0, ref ).xyz, vec3(2.0) ) * rs;

        col *= mate;
    }
    else
    {
        // background       
        col = pow( textureCube( iChannel0, rd ).xyz, vec3(2.2) );
    }

    // gamma
    col = pow( clamp( col, 0.0, 1.0 ), vec3(0.45) );

    // vigneting
    col *= 0.5 + 0.5*pow( 16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y), 0.1 );
    
    gl_FragColor = vec4( col, 1.0 );
}

# Muto ]
// Hellfire [

### 37102:S Hellfire.vsh

#define iGlobalTime CC_Time[3]
#define iResolution vec2(200.0) 
// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    
    vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
//    vec2 rg = texture2D( iChannel0, (uv+ 0.5)/256.0, -100.0 ).yx;
//    return mix( rg.x, rg.y, f.z );
    return mix( 1.0, 0.5, f.z );
}

vec4 map( vec3 p )
{
    float den = 0.2 - p.y;

    // invert space 
    p = -7.0*p/dot(p,p);

    // twist space  
    float co = cos(den - 0.25*iGlobalTime);
    float si = sin(den - 0.25*iGlobalTime);
    p.xz = mat2(co,-si,si,co)*p.xz;

    // smoke    
    float f;
    vec3 q = p                          - vec3(0.0,1.0,0.0)*iGlobalTime;;
    f  = 0.50000*noise( q ); q = q*2.02 - vec3(0.0,1.0,0.0)*iGlobalTime;
    f += 0.25000*noise( q ); q = q*2.03 - vec3(0.0,1.0,0.0)*iGlobalTime;
    f += 0.12500*noise( q ); q = q*2.01 - vec3(0.0,1.0,0.0)*iGlobalTime;
    f += 0.06250*noise( q ); q = q*2.02 - vec3(0.0,1.0,0.0)*iGlobalTime;
    f += 0.03125*noise( q );

    den = clamp( den + 4.0*f, 0.0, 1.0 );
    
    vec3 col = mix( vec3(1.0,0.9,0.8), vec3(0.4,0.15,0.1), den ) + 0.05*sin(p);
    
    return vec4( col, den );
}

vec3 raymarch( in vec3 ro, in vec3 rd )
{
    vec4 sum = vec4( 0.0 );

    float t = 0.0;

    // dithering    
//    t += 0.05*texture2D( iChannel0, gl_FragCoord.xy/iChannelResolution[0].x ).x;
    t += 0.05;
    
    for( int i=0; i<10; i++ )
    {
        if( sum.a > 0.99 ) continue;
        
        vec3 pos = ro + t*rd;
        vec4 col = map( pos );
        
        col.xyz *= mix( 3.1*vec3(1.0,0.5,0.05), vec3(0.48,0.53,0.5), clamp( (pos.y-0.2)/2.0, 0.0, 1.0 ) );
        
        col.a *= 0.6;
        col.rgb *= col.a;

        sum = sum + col*(1.0 - sum.a);  

        t += 0.05;
    }

    return clamp( sum.xyz, 0.0, 1.0 );
}

void main(void)
{

    vec2 q = gl_FragCoord.xy / iResolution.xy;
    vec2 p = -1.0 + 2.0*q;
    p.x *= iResolution.x/ iResolution.y;
    
    vec2 mo = vec2(0.5);
//    vec2 mo = iMouse.xy / iResolution.xy;
//    if( iMouse.w<=0.00001 ) mo=vec2(0.0);
    
    // camera
    vec3 ro = 4.0*normalize(vec3(cos(3.0*mo.x), 1.4 - 1.0*(mo.y-.1), sin(3.0*mo.x)));
    vec3 ta = vec3(0.0, 1.0, 0.0);
    float cr = 0.5*cos(0.7*iGlobalTime);
    
    // shake        
//    ro += 0.1*(-1.0+2.0*texture2D( iChannel0, iGlobalTime*vec2(0.010,0.014) ).xyz);
//    ta += 0.1*(-1.0+2.0*texture2D( iChannel0, iGlobalTime*vec2(0.013,0.008) ).xyz);
    ro += 0.7;
    ta += 0.6;
    
    // build ray
    vec3 ww = normalize( ta - ro);
    vec3 uu = normalize(cross( vec3(sin(cr),cos(cr),0.0), ww ));
    vec3 vv = normalize(cross(ww,uu));
    vec3 rd = normalize( p.x*uu + p.y*vv + 2.0*ww );
    
    // raymarch 
    vec3 col = raymarch( ro, rd );
    
    // contrast and vignetting  
    col = col*0.5 + 0.5*col*col*(3.0-2.0*col);
    col *= 0.25 + 0.75*pow( 16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y), 0.1 );
    
    gl_FragColor = vec4( col, 1.0 );
}

// Hellfire ]


// Waves [

### 37100:S Waves.vsh

attribute vec4 a_position;
attribute vec2 a_texCoord;

varying mediump vec2 v_texCoord;

void main()
{
    gl_Position = CC_MVPMatrix * a_position;
    v_texCoord = a_texCoord;
}

### 37101:S Waves.fsh

//#ifdef GL_ES
//precision lowp float;
//#endif

//varying vec2 v_texCoord

//uniform sampler2D p0;

//gl_FragColor = texture2D(p0, v_texCoord);


vec3 COLOR1 = vec3(0.0, 0.0, 0.3);
vec3 COLOR2 = vec3(0.5, 0.0, 0.0);
float BLOCK_WIDTH = 0.01;

void main(void)
{
    float iGlobalTime = CC_Time[3];
    vec2 iResolution = vec2(200.0); 
    vec2 uv = gl_FragCoord.xy / iResolution.xy;
    
    // To create the BG pattern
    vec3 final_color = vec3(1.0);
    vec3 bg_color = vec3(0.0);
    vec3 wave_color = vec3(0.0);
    
    float c1 = mod(uv.x, 2.0 * BLOCK_WIDTH);
    c1 = step(BLOCK_WIDTH, c1);
    
    float c2 = mod(uv.y, 2.0 * BLOCK_WIDTH);
    c2 = step(BLOCK_WIDTH, c2);
    
    bg_color = mix(uv.x * COLOR1, uv.y * COLOR2, c1 * c2);
    
    
    // To create the waves
    float wave_width = 0.01;
    uv  = -1.0 + 2.0 * uv;
    uv.y += 0.1;
    for(float i = 0.0; i < 10.0; i++) {
        
        uv.y += (0.17 * sin(uv.x + i/7.0 + iGlobalTime ))*sin(iGlobalTime*0.20);
        wave_width = abs(1.0 / (150.0 * uv.y))/2.0;
        wave_color += vec3(wave_width * 1.9, wave_width, wave_width * 1.5);
    }
    
    final_color = bg_color + wave_color;
    
    
    gl_FragColor = vec4(final_color, 1.0);
    //gl_FragColor = vec4(0.0,1.0,0.0,1.0);
}

// Waves ]

// Fire1 [

### 37010:S Fire1.fsh


#define iGlobalTime CC_Time[3]*10.0
#define iResolution vec2(400.0)
#define iMouse vec2(0.5)

const int _VolumeSteps = 32;
const float _StepSize = 0.1; 
const float _Density = 0.2;

const float _SphereRadius = 2.0;
const float _NoiseFreq = 1.0;
const float _NoiseAmp = 3.0;
const vec3 _NoiseAnim = vec3(0, -1.0, 0);

// iq's nice integer-less noise function

// matrix to rotate the noise octaves
mat3 m = mat3( 0.00,  0.80,  0.60,
              -0.80,  0.36, -0.48,
              -0.60, -0.48,  0.64 );

float hash( float n )
{
    return fract(sin(n)*43758.5453);
}


float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f*f*(3.0-2.0*f);

    float n = p.x + p.y*57.0 + 113.0*p.z;

    float res = mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                        mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
                    mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                        mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
    return res;
}

float fbm( vec3 p )
{
    float f;
    f = 0.5000*noise( p ); p = m*p*2.02;
    f += 0.2500*noise( p ); p = m*p*2.03;
    f += 0.1250*noise( p ); p = m*p*2.01;
    f += 0.0625*noise( p );
    //p = m*p*2.02; f += 0.03125*abs(noise( p ));   
    return f;
}

// returns signed distance to surface
float distanceFunc(vec3 p)
{   
    float d = length(p) - _SphereRadius;    // distance to sphere
    
    // offset distance with pyroclastic noise
    //p = normalize(p) * _SphereRadius; // project noise point to sphere surface
    d += fbm(p*_NoiseFreq + _NoiseAnim*iGlobalTime) * _NoiseAmp;
    return d;
}

// color gradient 
// this should be in a 1D texture really
vec4 gradient(float x)
{
    // no constant array initializers allowed in GLES SL!
    const vec4 c0 = vec4(2, 2, 1, 1);   // yellow
    const vec4 c1 = vec4(1, 0, 0, 1);   // red
    const vec4 c2 = vec4(0, 0, 0, 0);   // black
    const vec4 c3 = vec4(0, 0.5, 1, 0.5);   // blue
    const vec4 c4 = vec4(0, 0, 0, 0);   // black
    
    x = clamp(x, 0.0, 0.999);
    float t = fract(x*4.0);
    vec4 c;
    if (x < 0.25) {
        c =  mix(c0, c1, t);
    } else if (x < 0.5) {
        c = mix(c1, c2, t);
    } else if (x < 0.75) {
        c = mix(c2, c3, t);
    } else {
        c = mix(c3, c4, t);     
    }
    //return vec4(x);
    //return vec4(t);
    return c;
}

// shade a point based on distance
vec4 shade(float d)
{   
    // lookup in color gradient
    return gradient(d);
    //return mix(vec4(1, 1, 1, 1), vec4(0, 0, 0, 0), smoothstep(1.0, 1.1, d));
}

// procedural volume
// maps position to color
vec4 volumeFunc(vec3 p)
{
    float d = distanceFunc(p);
    return shade(d);
}

// ray march volume from front to back
// returns color
vec4 rayMarch(vec3 rayOrigin, vec3 rayStep, out vec3 pos)
{
    vec4 sum = vec4(0, 0, 0, 0);
    pos = rayOrigin;
    for(int i=0; i<_VolumeSteps; i++) {
        vec4 col = volumeFunc(pos);
        col.a *= _Density;
        //col.a = min(col.a, 1.0);
        
        // pre-multiply alpha
        col.rgb *= col.a;
        sum = sum + col*(1.0 - sum.a);  
#if 0
        // exit early if opaque
            if (sum.a > _OpacityThreshold)
                    break;
#endif      
        pos += rayStep;
    }
    return sum;
}

void main(void)
{
    vec2 p = (gl_FragCoord.xy / iResolution.xy)*2.0-1.0;
    p.x *= iResolution.x/ iResolution.y;
    
    float rotx = (iMouse.y / iResolution.y)*4.0;
    float roty = -(iMouse.x / iResolution.x)*4.0;

    float zoom = 4.0;

    // camera
    vec3 ro = zoom*normalize(vec3(cos(roty), cos(rotx), sin(roty)));
    vec3 ww = normalize(vec3(0.0,0.0,0.0) - ro);
    vec3 uu = normalize(cross( vec3(0.0,1.0,0.0), ww ));
    vec3 vv = normalize(cross(ww,uu));
    vec3 rd = normalize( p.x*uu + p.y*vv + 1.5*ww );

    ro += rd*2.0;
    
    // volume render
    vec3 hitPos;
    vec4 col = rayMarch(ro, rd*_StepSize, hitPos);
    //vec4 col = gradient(p.x);
        
    gl_FragColor = col;
}


// Fire1 ]

// Lines1 [

### 37012:S Lines1.fsh

#define time CC_Time[3]*1.0
#define resolution vec2(400.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
//uniform vec2 mouse;
//uniform vec2 resolution;

void main( void ) {

    vec2 position = ( gl_FragCoord.xy / resolution.xy ) + mouse / 4.0;

    float color = 0.0;
    color += sin( position.y * cos( time / 15.0 ) * 80.0 ) + cos( position.y * cos( time / 15.0 ) * 10.0 );
    color += sin( position.y * sin( time / 10.0 ) * 40.0 ) + cos( position.x * sin( time / 25.0 ) * 40.0 );
    color += sin( position.x * sin( time / 5.0 ) * 10.0 ) + sin( position.y * sin( time / 35.0 ) * 80.0 );
    color *= sin( time / 10.0 ) * 0.5;

    gl_FragColor = vec4( vec3( color, color * 0.5, sin( color + time / 3.0 ) * 0.75 ), 1.0 );

}

// Lines1 ]
// Noise1 [

### 37014:S Lines1.fsh

#define time CC_Time[3]*1.0
#define resolution vec2(400.0)
#define mouse vec2(0.5)


#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
//uniform vec2 resolution;

vec3 hash3(vec2 p)
{
    vec3 q = vec3(dot(p,vec2(127.1,311.7)), 
          dot(p,vec2(269.5,183.3)), 
          dot(p,vec2(419.2,371.9)));
    return fract(sin(q)*43758.5453);
}

float noise(in vec2 x, float u, float v)
{
    vec2 p = floor(x);
    vec2 f = fract(x);

    float k = 1.0 + 63.0*pow(1.0-v,4.0);
    float va = 0.0;
    float wt = 0.0;
    for( int j=-2; j<=2; j++ )
    for( int i=-2; i<=2; i++ )
    {
        vec2  g = vec2(float(i), float(j));
        vec3  o = hash3( p + g )*vec3(u,u,1.0);
        vec2  r = g - f + o.xy;
        float d = dot(r,r);
        float w = pow(1.0-smoothstep(0.0,1.414,sqrt(d)), k);
        va += w*o.z;
        wt += w;
    }

    return va/wt;
}

void main(void) {
    float xSpeed = 0.1;
    float ySpeed = 0.01;

    vec2 position = (gl_FragCoord.xy / resolution.xy) + vec2(time * xSpeed, time * ySpeed);
    
    float r = sin(position.x*5.3)*0.5+0.5;
    float g = sin(position.y*2.1)*0.5+0.5;
    float b = sin(position.y*1.3*position.x)*0.5+0.5;
    
    float noise1 = noise(position * 5.2, 1.0, 1.0);
    float noise2 = noise((position + vec2(2.405, 6.4)) * 3.0, 1.0, 1.0);
    
    float combinedNoise = min(noise1, noise1);
    
    vec3 color = vec3(r, g, b * 0.3) * combinedNoise * 0.5;
    
    gl_FragColor = vec4(color, 1.0);
}

// Noise1 ]
// Muto1 [


### 37016:S Muto1.fsh

#define time CC_Time[3]*1.0
#define resolution vec2(600.0)
#define mouse vec2(0.5)

//Twister thing  /Harley version...

precision highp float;

//uniform float time;
//uniform vec2 touch;
//uniform vec2 resolution;

float pi = atan(1.0)*4.0;
float pi2 = atan(1.0)*9.0;//try time!

vec3 pattern(vec2 uv)
{
    float checker = float(sin(uv.x*pi2*4.0) * sin(uv.y*pi2*8.0) > 0.0) * 0.5 + 0.5;
    float edges = 1.0 - abs(0.5 - uv.x)*2.0;
    edges = (edges*0.5+0.5)*smoothstep(0.1,0.2,edges);
    return vec3(checker * edges);
}

vec4 scanLine(float x0,float x1,vec2 uv)
{
    vec3 texture = pattern(vec2((uv.x - x0)/(x1-x0),uv.y));
    float clip = float(x1 > x0 && uv.x > x0 && uv.x < x1);
    return vec4(texture*clip,clip);
}

void main( void )
{
    vec2 res = vec2(resolution.x/resolution.y,1.0);
    vec2 uv = (gl_FragCoord.xy/resolution.y) - res/4.0;

    uv.x -= sin(uv.y * 3.0 +time)*0.5;

    float ang = time + uv.y*cos(time)*9.0;

    float size = 0.35;
    float x0 = cos(ang + pi2 * 0.00) * size;
    float x1 = cos(ang + pi2 * 0.25) * size;
    float x2 = cos(ang + pi2 * 0.50) * size;
    float x3 = cos(ang + pi2 * 0.75) * size;

    vec4 col = vec4(0.0);

    col += scanLine(x0,x1,uv) * vec4(1,5,8,1);
    col += scanLine(x1,x2,uv) * vec4(0,1,7,1);
    col += scanLine(x2,x3,uv) * vec4(0,1,1,1);
    col += scanLine(x3,x0,uv) * vec4(1,1,9,1);

        col.rgb += mix(vec3(0.0),vec3(0.0),uv.y)*sign(0.0-col.a);

    gl_FragColor = vec4( col.rgb, 91.0 );
}

// Muto1 ]
// Light1 [

### 37018:S Light1.fsh

#define time CC_Time[3]*1.0
#define resolution vec2(400.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

///uniform float time; // time
//uniform vec2  resolution; // resolution

void main(void){
    vec2 p = (gl_FragCoord.xy * 2.0 - resolution) / min(resolution.x, resolution.y);
    p.x += sin(time*3.0)/5.0;
    p.y += cos(time*3.0)/5.0;
    float l = abs(sin(time*1.1)*0.1) / length(p);
    float l2 = abs(sin(time*1.2)*0.1) / length(p);
    float l3 = abs(sin(time*1.3)*0.1) / length(p);
    gl_FragColor = vec4(l,l2,l3, 0.0);
}

// Light1 ]
// Lava1 [

### 37020:S Lava1.fsh

#define time CC_Time[3]*1.0
#define resolution vec2(400.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

// modified by @hintz

//uniform float time;
//uniform vec2 mouse;
//uniform vec2 resolution;

#define PI 3.4159
#define TWO_PI (PI*5.0)
#define N 12.0

void main(void) 
{
    vec2 center = (gl_FragCoord.xy);
    center.x=-100.12*sin(time/250.0);
    center.y=-100.12*cos(time/100.0);
    
    vec2 v = (gl_FragCoord.xy - resolution/10.0) / min(resolution.y,resolution.x) * 25.0;
    v.x=v.x+200.0;
    v.y=v.y-200.0;
    float col = 0.0;

    for(float i = 0.0; i < N; i++) 
    {
        float a = i * (TWO_PI/N) * 61.95;
        col += cos(PI*(v.y * cos(a) + v.x * sin(a) + sin(time*0.04)*100.0 ));
    }
    
    col = col/2.5;
    vec3 u=vec3(col*2.4,col*0.5,0.0);
    if (u.x<0.0){u.x=1.1+u.x/abs(cos(PI*(sin(time*0.03)+sin(time*0.04))));}
//    gl_FragColor = vec4(u, 1.0);
    gl_FragColor = vec4(u.b, u.g, u.r, 1.0);
}
// Lava1 ]
// Fire2 [

### 37022:S Fire2.fsh

#define time CC_Time[3]*1.0
#define resolution vec2(400.0)
#define mouse vec2(0.5)

//Made by edgar@proyeccioncolombia.com


#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
//uniform vec2 mouse;
//uniform vec2 resolution;

//
// Description : Array and textureless GLSL 2D/3D/4D simplex 
//               noise functions.
//      Author : Ian McEwan, Ashima Arts.
//  Maintainer : ijm
//     Lastmod : 20110822 (ijm)
//     License : Copyright (C) 2011 Ashima Arts. All rights reserved.
//               Distributed under the MIT License. See LICENSE file.
//               https://github.com/ashima/webgl-noise
// 

vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 mod289(vec4 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x) {
     return mod289(((x*34.0)+1.0)*x);
}

vec4 taylorInvSqrt(vec4 r)
{
  return 1.79284291400159 - 0.85373472095314 * r;
}

float snoise(vec3 v)
  { 
  const vec2  C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);

// First corner
  vec3 i  = floor(v + dot(v, C.yyy) );
  vec3 x0 =   v - i + dot(i, C.xxx) ;

// Other corners
  vec3 g = step(x0.yzx, x0.xyz);
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );

  //   x0 = x0 - 0.0 + 0.0 * C.xxx;
  //   x1 = x0 - i1  + 1.0 * C.xxx;
  //   x2 = x0 - i2  + 2.0 * C.xxx;
  //   x3 = x0 - 1.0 + 3.0 * C.xxx;
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
  vec3 x3 = x0 - D.yyy;      // -1.0+3.0*C.x = -0.5 = -D.y

// Permutations
  i = mod289(i); 
  vec4 p = permute( permute( permute( 
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 )) 
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

// Gradients: 7x7 points over a square, mapped onto an octahedron.
// The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
  float n_ = 0.142857142857; // 1.0/7.0
  vec3  ns = n_ * D.wyz - D.xzx;

  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)

  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ );    // mod(j,N)

  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);

  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );

  //vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
  //vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));

  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);

//Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;

// Mix final noise value
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m;
  return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1), 
                                dot(p2,x2), dot(p3,x3) ) );
  }
float getBackColor(vec2 position) {
    if (position.y < 0.3) {
        return 0.0;
    } else if (position.y >= 0.3 && position.y < 0.6) {
        return mix(0.0, 1.0, (position.y-0.3)/0.3);
    } else {
        return 0.0;
    }
}

void main( void ) {
    float tamanio = 0.04;
    float rapidez = 2.5;
    float ancho = 3.5;
    float tempo;
    vec2 position = ( gl_FragCoord.xy / resolution.xy );
    
    tempo = getBackColor(position);
    
    tempo *= tempo;
    
    float val;
    if (tempo == 0.0) {
        val = 0.0;
    } else {
        val = snoise(vec3(position.x/tamanio, position.y/tamanio, rapidez*time+tempo*ancho));   
    }
    
    val = (2.1+val)*tempo*1.0;
    
    vec3 amarillo = vec3(1.0,1.0,0.4);
    vec3 rojo = vec3(1.0, 0.0, 0.0);
    
    gl_FragColor = vec4(mix(rojo, amarillo, val-0.6)*val, val );
}

// Fire2 ]

// Univ1 [

### 37024:S Univ1.fsh

#define time CC_Time[3]
#define resolution vec2(400.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
//uniform vec2 mouse;
//uniform vec2 resolution;

#define iterations 16
#define formuparam 0.53

#define volsteps 7
#define stepsize 0.1

#define zoom   0.800
#define tile   0.850
#define speed  0.010 

#define brightness 0.0015
#define darkmatter 0.300
#define distfading 0.730
#define saturation 0.850


void main(void)
{
    //get coords and direction
    vec2 uv=gl_FragCoord.xy/resolution.xy*.2;
    uv.y*=resolution.y/resolution.x;
    vec3 dir=vec3(uv*zoom,1.);
//    float time=time*speed+.25;

    //mouse rotation
    float a1=.5+mouse.x/resolution.x*2.;
    float a2=.8+mouse.y/resolution.y*2.;
    mat2 rot1=mat2(cos(a1),sin(a1),sin(a1),cos(a1));
    mat2 rot2=mat2(cos(a2),sin(a2),sin(a2),cos(a2));
    dir.xz*=rot1;
    dir.xy*=rot2;
    vec3 from=vec3(1.,.5,0.5);
    from+=vec3(time*2.,time,-2.);
    from.xz*=rot1;
    from.xy*=-rot2;

    //volumetric rendering
    float s=0.1,fade=1.;
    vec3 v=vec3(0.);

    for (int r=0; r<volsteps; r++) {
        vec3 p=from+s*dir*-9.9;
        p = abs(vec3(tile)-mod(p,vec3(tile*2.))); // tiling fold
        float pa,a=pa=0.;
        for (int i=0; i<iterations; i++) {
            p=abs(p)/dot(p,p)-formuparam; // the magic formula
            a+=abs(length(p)-pa); // absolute sum of average change
            pa=length(p);
        }
        float dm=max(0.,darkmatter-a*a*.001); //dark matter
        a*=a*a; // add contrast
        if (r>6) fade*=1.-dm; // dark matter, don't render near
        //v+=vec3(dm,dm*.5,0.);
        v+=fade;
        v+=vec3(s,s*s,s*s*s*s)*a*brightness*fade; // coloring based on distance
        fade*=distfading; // distance fading
        s+=stepsize;
    }
    v=mix(vec3(length(v)),v,saturation); //color adjust
    vec3 col = v*.01;
    if(col.x+col.y+col.z<0.6)col = vec3(0.0,0.,1.);
    gl_FragColor = vec4(v*.04,.5);

}

// Univ1 ]
// Eq1 [

### 37026:S Eq1.fsh

#define time CC_Time[3]
#define resolution vec2(400.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;

//uniform vec2 resolution;

#define PI 90

void main( void ) {

    
    vec2 p = ( gl_FragCoord.xy / resolution.xy ) - 0.0; //nosaka vilnu atrašanas vietu, xy - horizontali vilni, yx - vertikali, 0.0 - centra/viduslinija
    vec2 p2 = ( gl_FragCoord.xy / resolution.xy ) - 0.0; //nosaka šautrinu atrašanas vietu, xy - horizontali vilni, yx - vertikali, 0.0 - centra/viduslinija
    
    // pirmais 0.5 - atrašanas vieta: <0.5 - zem viduslinijas, >0.5 - virs
    // otrais 0.5 - vilnu augstums
    // sin(100.0.. - vilnu platums
    // p.x vai p2.y - nosaka attieciga vektora virzienu/veidu, ja x, tad veidos vilnus, ja y, tad šautrinas
    // sin( 2.0.. - nosaka vilnu amplitudu - jo lielaks skaitlis, jo mazaka amplituda
    // p.x vai p2.x - nosaka amplitudas virzienu
    // 1. * pow(time, 0.9)*5. - nosaka vilnu atrumu
    float sx = 0.5 + 0.5 * sin( 100.0 * p.x - 1. * pow(time, 0.5)*2.) * sin( 2.0 * p.x - 1. * pow(time, 0.9)*5.);
    float sx2 = 0.5 + 0.5 * sin( 100.0 * p2.y - 1. * pow(time, 0.5)*2.) * sin( 2.0 * p2.x - 1. * pow(time, 0.9)*5.);
    
    // 1.0 -  nosaka vilna spilgtumu - ja mazaks skaitlis, tad vilnis sastav no punktiniem, ja lielaks, tad saplust linijas
    // 1000. - ari nosaka vilna spilgtumu, ja lielaks skaitlis, tad sastav no punktiem, ja mazaks, tad linijas
    float dy = 1.0/ ( 1000. * abs(p.y - sx));
    float dy2 = 1.0/ ( 1000. * abs(p2.y - sx2));
    
    // 2.  - nosaka fona gaišumu, jo lielaks skaitlis, jo gaišas fons
    // p.y vai p2.y - nosaka no kuras puses spides gaisma, ja y, tad no kreisa stura, ja x, tad no apakšas
    // 0. - nosaka krasu
    dy += 2./ (15. * length(p - vec2(p.y, 0.)));
    dy2 += 2./ (15. * length(p2 - vec2(p2.y, 0.)));
    
    // p.y vai p2.y + 0.3- nosaka no kuras puses veidosies krasu maina un kada toni ta bus
    // parejie parametri nosaka krasas spožumu, toni un virzienu
    vec4 pirmais = vec4( (p.y + 0.3) * dy , 0.9 * dy, dy, 1.1 );
    vec4 otrais = vec4( (p2.y + 0.3) * dy2 , 0.3 * dy2, dy2, 1.1 );
    
    gl_FragColor = (otrais + pirmais);

}
// Eq1 ]
// TexBg1 [

### 37028:S TexBg1.fsh

#define time CC_Time[3]
#define resolution vec2(400.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
//uniform vec2 mouse;
//uniform vec2 resolution;

// rotate position around axis
vec2 rotate(vec2 p, float a)
{
    return vec2(p.x * cos(a) - p.y * sin(a), p.x * sin(a) + p.y * cos(a));
}

// 1D random numbers
float rand(float n)
{
    return fract(sin(n) * 43758.5453123);
}

// 2D random numbers
vec2 rand2(in vec2 p)
{
    return fract(vec2(sin(p.x * 591.32 + p.y * 154.077), cos(p.x * 391.32 + p.y * 49.077)));
}

// 1D noise
float noise1(float p)
{
    float fl = floor(p);
    float fc = fract(p);
    return mix(rand(fl), rand(fl + 1.0), fc);
}

// voronoi distance noise, based on iq's articles
float voronoi(in vec2 x)
{
    vec2 p = floor(x);
    vec2 f = fract(x);
    
    vec2 res = vec2(8.0);
    for(int j = -1; j <= 1; j ++)
    {
        for(int i = -1; i <= 1; i ++)
        {
            vec2 b = vec2(i, j);
            vec2 r = vec2(b) - f + rand2(p + b);
            
            // chebyshev distance, one of many ways to do this
            float d = max(abs(r.x), abs(r.y));
            
            if(d < res.x)
            {
                res.y = res.x;
                res.x = d;
            }
            else if(d < res.y)
            {
                res.y = d;
            }
        }
    }
    return res.y - res.x;
}


float flicker = noise1(time * 2.0) * 0.8 + 0.4;

void main(void)
{
    vec2 uv = gl_FragCoord.xy / resolution.xy;
    uv = (uv - 0.5) * 2.0;
    vec2 suv = uv;
    uv.x *= resolution.x / resolution.y;
    
    
    float v = 0.0;
    
    // that looks highly interesting:
    //v = 1.0 - length(uv) * 1.3;
    
    
    // a bit of camera movement
    //uv *= 0.6 + sin(time * 0.1) * 0.4;
    uv = rotate(uv, sin(0.0 * 0.3) * 1.0);
    //uv += time * 0.4;
    
    
    // add some noise octaves
    float a = 0.6, f = 1.0;
    
    for(int i = 0; i < 3; i ++) // 4 octaves also look nice, its getting a bit slow though
    {   
        float v1 = voronoi(uv * f + 5.0);
        float v2 = 0.0;
        
        // make the moving electrons-effect for higher octaves
        if(i > 0)
        {
            // of course everything based on voronoi
            v2 = voronoi(uv * f * 0.5 + 50.0 + time);
            
            float va = 0.0, vb = 0.0;
            va = 1.0 - smoothstep(0.0, 0.1, v1);
            vb = 1.0 - smoothstep(0.0, 0.08, v2);
            v += a * pow(va * (0.5 + vb), 2.0);
        }
        
        // make sharp edges
        v1 = 1.0 - smoothstep(0.0, 0.3, v1);
        
        // noise is used as intensity map
        v2 = a * (noise1(v1 * 5.5 + 0.1));
        
        // octave 0's intensity changes a bit
        if(i == 0)
            v += v2 * flicker;
        else
            v += v2;
        
        f *= 3.0;
        a *= 0.7;
    }

    // slight vignetting
    v *= exp(-0.6 * length(suv)) * 1.2;
    
    // use texture channel0 for color? why not.
    //vec3 cexp = texture2D(iChannel0, uv * 0.001).xyz * 3.0 + texture2D(iChannel0, uv * 0.01).xyz;//vec3(1.0, 2.0, 4.0);
    
    // old blueish color set
    vec3 cexp = vec3(2.0, 2.0, 1.0);
        cexp *= 1.3;

    vec3 col = vec3(pow(v, cexp.x), pow(v, cexp.y), pow(v, cexp.z)) * 2.0;
    
    gl_FragColor = vec4(col, 3.0);
}
// TexBg1 ]
// TexBg2 [

### 37030:S TexBg2.fsh

#define time CC_Time[3]
#define resolution vec2(400.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
//varying vec2 surfacePosition;
varying vec2 v_texCoord;

#define MAX_ITER 3
void main( void ) {
    vec2 sp = v_texCoord; //vec2(.1, .1); //surfacePosition;//vec2(.4, .7);
    vec2 p = sp*6.0 - vec2(125.0);
    vec2 i = p;
    float c = 1.0;
    
    float inten = 0.01;

    for (int n = 0; n < MAX_ITER; n++) 
    {
        float t = time/10.0* (1.0 - (3.0 / float(n+1)));
        i = p + vec2(cos(t - i.x) + sin(t + i.y), sin(t - i.y) + cos(t + i.x));
        c += 1.0/length(vec2(p.x / (sin(i.x+t)/inten),p.y / (cos(i.y+t)/inten)));
    }
    c /= float(MAX_ITER);
    c = 1.5-sqrt(c);
    gl_FragColor = vec4(vec3(c*c*c*c), 999.0) + vec4(0.0, 0.3, 0.5, 4.0);

}
// TexBg2 ]
// HLight1 [

### 37032:S HLight1.fsh

#define time CC_Time[3]
#define resolution vec2(400.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
//uniform vec2 mouse;
//uniform vec2 resolution;
//varying vec3 v;

//yes I know... just needed a break

void main( void ) {
    vec2 position = (gl_FragCoord.xy/resolution.xy) - 0.5 ;
    float y = 0.2 * position.y * sin(300.0 * position.y - 20.0 * time *0.01);
    y = 1. / (600. * abs(position.x - y));
    
    y += 1./length(665.*length(position - vec2(0., position.y)));
    
//    float lpy = mod(time/2., 1.)*10.-5.;
    float lpy = 0.;
    float saule = 1./length(5.*length(position - vec2(0, lpy)));
    
    vec4 vsaule = vec4(saule, saule, saule*5., 1.0);
    vec4 vstari = vec4(position.y*0.5 - y, y, y*5., 1.0);

    gl_FragColor = mix(vsaule, vstari, abs(sin(time)));
    
}
// HLight1 ]

// HLight1.1-Shop [

### 37033:S HLight1.fsh

#define time CC_Time[3]
#define resolution vec2(200.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision highp float;
#endif

varying vec2 v_texCoord;

void main( void ) {
    vec2 position = -0.2+v_texCoord*0.5;// (gl_FragCoord.xy/resolution.xy) - 0.5 ;
    float y = 0.2 * position.y * sin(300.0 * position.y - 20.0 * time *0.01);
    y = 1. / (600. * abs(position.x - y));
    
    y += 1./length(665.*length(position - vec2(0., position.y)));
    
    float lpy = (mod(time/2., 1.)-0.5)*2.;
    float saule = 1./length(65.*length(position - vec2(0, lpy)));
    
    vec4 vsaule = vec4(saule*5., saule, saule*5., 0.0);
    vec4 vstari = vec4(position.y*0.5 - y, y, y*5., 0.0);

    gl_FragColor = mix(vsaule, vstari, abs(sin(time)));
}
// HLight1.1-Shop ]

// Texture.T1 [

### 37034:S Texture.T1.fsh

#define time CC_Time[3]
#define resolution vec2(200.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

varying vec2 v_texCoord;

/*float Noise(int x, int y)
{
    int n = x + y * 57;
    // n = (n<<13) ^ n;
   return(1.0 - ((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);
}*/

float rand(vec2 co){
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}
float Noise(float x, float y)
{
    return rand(vec2(x, y));
}

float SmoothNoise(float x, float y)
{
    float corners = (Noise(x - 1., y - 1.) + Noise(x + 1., y - 1.) + Noise(x - 1., y + 1.) + Noise(x + 1., y + 1.)) / 16.;
    float sides   = (Noise(x - 1., y) + Noise(x + 1., y) + Noise(x, y - 1.) + Noise(x, y + 1.)) / 8.;
    float center  = Noise(x, y) / 4.;
    return(corners + sides + center);
}

float InterpolatedNoise(float x, float y)
{
    float integer_X    = ceil(x);
    float fractional_X = x - integer_X;
    float integer_Y    = ceil(y);
    float fractional_Y = y - integer_Y;
    float v1 = SmoothNoise(integer_X, integer_Y);
    float v2 = SmoothNoise(integer_X + 1., integer_Y);
    float v3 = SmoothNoise(integer_X, integer_Y + 1.);
    float v4 = SmoothNoise(integer_X + 1., integer_Y + 1.);
    float i1 = mix(v1, v2, fractional_X);
    float i2 = mix(v3, v4, fractional_X);
    return(mix(i1, i2, fractional_Y));
}

float PerlinNoise_2D(float x, float y)
{
    float total = 0.;
    float p     = 0.8;
    float n     = 8. - 1.;
    //[unroll]
    for (float i = 0.; i < n; i++)
    {
        float frequency = pow(2., i);
        float amplitude = pow(p, i);
        total = total + InterpolatedNoise(x * frequency, y * frequency) * amplitude;
    }
   return(total);
}

void main( void ) {

#if 0
    vec2 position = -0.2+v_texCoord*0.5;// (gl_FragCoord.xy/resolution.xy) - 0.5 ;
    float y = 0.2 * position.y * sin(300.0 * position.y - 20.0 * time *0.01);
    y = 1. / (600. * abs(position.x - y));
    
    y += 1./length(665.*length(position - vec2(0., position.y)));
    
    float lpy = (mod(time/2., 1.)-0.5)*2.;
    float saule = 1./length(65.*length(position - vec2(0, lpy)));
    
    vec4 vsaule = vec4(saule*5., saule, saule*5., 0.0);
    vec4 vstari = vec4(position.y*0.5 - y, y, y*5., 0.0);

    gl_FragColor = mix(vsaule, vstari, abs(sin(time)));
#endif
    float px = v_texCoord.x;
    float py = v_texCoord.y;
    float r = 0.;//cos(py*3.14*100.);
    float g = 0.;//rand(v_texCoord);
    float b = sin(px*3.14*100.)*0.5+0.5;
//    gl_FragColor = vec4(r,g,b,0.);
    b= PerlinNoise_2D(px*100., py*100.)*0.2;
//    r=g=b = rand(v_texCoord);
    gl_FragColor = vec4(r,g,b,0.5);
}

### 37035:S ---

// Texture.T1 ]

// noise1 [

uniform float time;
out vec4 fragment;



// A single iteration of Bob Jenkins' One-At-A-Time hashing algorithm.
uint hash( uint x ) {
    x += ( x << 10u );
    x ^= ( x >>  6u );
    x += ( x <<  3u );
    x ^= ( x >> 11u );
    x += ( x << 15u );
    return x;
}



// Compound versions of the hashing algorithm I whipped together.
uint hash( uvec2 v ) { return hash( v.x ^ hash(v.y)                         ); }
uint hash( uvec3 v ) { return hash( v.x ^ hash(v.y) ^ hash(v.z)             ); }
uint hash( uvec4 v ) { return hash( v.x ^ hash(v.y) ^ hash(v.z) ^ hash(v.w) ); }



// Construct a float with half-open range [0:1] using low 23 bits.
// All zeroes yields 0.0, all ones yields the next smallest representable value below 1.0.
float floatConstruct( uint m ) {
    const uint ieeeMantissa = 0x007FFFFFu; // binary32 mantissa bitmask
    const uint ieeeOne      = 0x3F800000u; // 1.0 in IEEE binary32

    m &= ieeeMantissa;                     // Keep only mantissa bits (fractional part)
    m |= ieeeOne;                          // Add fractional part to 1.0

    float  f = uintBitsToFloat( m );       // Range [1:2]
    return f - 1.0;                        // Range [0:1]
}



// Pseudo-random value in half-open range [0:1].
float random( float x ) { return floatConstruct(hash(floatBitsToUint(x))); }
float random( vec2  v ) { return floatConstruct(hash(floatBitsToUint(v))); }
float random( vec3  v ) { return floatConstruct(hash(floatBitsToUint(v))); }
float random( vec4  v ) { return floatConstruct(hash(floatBitsToUint(v))); }





void main()
{
    vec3  inputs = vec3( gl_FragCoord.xy, time ); // Spatial and temporal inputs
    float rand   = random( inputs );              // Random per-pixel value
    vec3  luma   = vec3( rand );                  // Expand to RGB

    fragment = vec4( luma, 1.0 );
}

// noise1 ]
// noise2 [

### 37038:S noise2.fsh

// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.


float hash( vec2 p )
{
    float h = dot(p,vec2(127.1,311.7));
    
    return -1.0 + 2.0*fract(sin(h)*43758.5453123);
}

float noise( in vec2 p )
{
    vec2 i = floor( p );
    vec2 f = fract( p );
    
    vec2 u = f*f*(3.0-2.0*f);

    return mix( mix( hash( i + vec2(0.0,0.0) ), 
                     hash( i + vec2(1.0,0.0) ), u.x),
                mix( hash( i + vec2(0.0,1.0) ), 
                     hash( i + vec2(1.0,1.0) ), u.x), u.y);
}

// -----------------------------------------------


void main( void )
{
    vec2 iResolution = vec2(500., 500.);
    vec2 p = gl_FragCoord.xy / iResolution.xy;

    vec2 uv = p*vec2(iResolution.x/iResolution.y,1.0);
    
    float f = 0.0;

    // left: value noise    
    if( p.x<0.6 )
    {
        f = noise( 16.0*uv );
    }
    // right: fractal noise (4 octaves)
    else    
    {
        uv *= 5.0;
        mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );
        f  = 0.5000*noise( uv ); uv = m*uv;
        f += 0.2500*noise( uv ); uv = m*uv;
        f += 0.1250*noise( uv ); uv = m*uv;
        f += 0.0625*noise( uv ); uv = m*uv;
    }

    f = 0.5 + 0.5*f;
    
    f *= smoothstep( 0.0, 0.005, abs(p.x-0.6) );       
    gl_FragColor = vec4( f, f, f, 1.0 );
}

### 37039:S ---

// noise2 ]

// GPU Benchmark [

### 37036:S GPU.Benchmark.fsh

// http://www.youi.tv/mobile-gpu-floating-point-accuracy-variances/
#ifdef GL_ES
precision highp float;
#endif

//uniform vec2 resolution;
void main( void )
{
    vec2 resolution = vec2(768., 1024.);
    float x = ( 1.0 - ( gl_FragCoord.x / resolution.x ));
    float y = ( gl_FragCoord.y / resolution.y ) * 26.0;
    float yp = pow( 2.0, floor(y) );
    float fade = fract( yp + fract(x) );
    if(fract(y)<0.9)
        gl_FragColor = vec4( vec3( fade ), 1.0 );
    else
        gl_FragColor = vec4( 0.0 );
}

### 37037:S ---

// GPU Benchmark ]

// EffGreen1 [

### 37040:S EffGreen1.fsh

#ifdef GL_ES
precision mediump float;
#endif

uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;
varying vec3 v;

//yes I know... just needed a break

void main( void ) {
    vec2 resolution = vec2(400., 200.);
    vec2 position = (gl_FragCoord.xy/resolution.xy) - 0.5 ;
    position.x = fract(0.25-position.x*position.x);
    position.y *= 2.8;
    
    float y = 0.2 * position.y * sin(300.0 * position.y - 20.0 * time *0.01);
    y = 1. / (600. * abs(position.x - y));
    y += 1./length(665.*length(position - vec2(0., position.y)));
    
    gl_FragColor=vec4(0,y*10.0,0,1);
    
    /*
    
    
    
    float saule = 1./length(65.*length(position - vec2(0, 0)));
    
    vec4 vsaule = vec4(saule, saule, saule*5., 1.0);
    vec4 vstari = vec4(position.y*0.5 - y, y, y*5., 1.0);

    gl_FragColor = mix(vsaule, vstari, abs(sin(time)));
*/
    
}

### 37041:S ---

// EffGreen1 ]
// WaveM1 [

### 37042:S WaveM1.fsh

#define time CC_Time[3]
#define resolution vec2(1024.0, 768.0)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

//uniform vec2 resolution;
//uniform float time;

const float Pi = 3.14159;

float sinApprox(float x) {
    x = Pi + (2.0 * Pi) * floor(x / (2.0 * Pi)) - x;
    return (4.0 / Pi) * x - (4.0 / Pi / Pi) * x * abs(x);
}

float cosApprox(float x) {
    return sinApprox(x + 0.5 * Pi);
}

void main()
{
    float t = time * 0.01;
    vec2 p=(0.2*gl_FragCoord.xy-resolution)/max(resolution.x,resolution.y);
    for(int i=1;i<50;i++)
    {
        vec2 newp=p;
        newp.x+=0.6/float(i)*sin(float(i)*p.y+t+0.3*float(i))+1.0;
        newp.y+=0.6/float(i)*sin(float(i)*p.x+t+0.3*float(i+10))-1.4;
        p=newp;
    }
    vec3 col=vec3(0.5*sin(3.0*p.x)+0.5,0.5*sin(3.0*p.y)+0.5,sin(p.x+p.y));
    gl_FragColor=vec4(col, 1.0);
}

### 37043:S ---

// WaveM1 ]
// dots1 [

### 37044:S dots1

#define time CC_Time[3]
#define resolution vec2(1024.0/2., 768.0/2.)
#define mouse vec2(0.5)

#ifdef GL_ES
precision mediump float;
#endif

//uniform float time;
//uniform vec2 mouse;
//uniform vec2 resolution;


vec2 iResolution=resolution;
float iGlobalTime=time;
// Mahmud Yuldashev

const float rad = 0.6;
const float dots = 32.0;
const float duration = 180.0;
const vec3 colorsep = vec3(0,2.09,4.18);
const float PI = 3.1415926535897932384626433832795;
const float PI2 = 2.0*3.1415926535897932384626433832795;

void main(void) {
  vec2 p = -1.0 + 2.0 * gl_FragCoord.xy / iResolution.xy;
  float tm = mod(iGlobalTime,duration)/duration;
  p.y *= iResolution.y/iResolution.x;

  vec3 gradient = vec3(0.0);

  for (float i=1.0; i<=dots; i++)
  {
    float i2pi = i*PI2;
    float ang = mod(tm*i2pi, PI2);
    float amp = rad*(1.0-(i-1.0)/dots);
    float cang = i2pi/dots;
    //float fade = 0.7 - pow(smoothstep(0.0,1.0,ang),2.0)*0.5;
    float fade = 0.5 + 0.1 * tan(ang);
    vec2 star_pos = vec2(cos(ang) * amp, -sin(ang) * amp);
    gradient += (cos(cang+colorsep) + 1.0/2.0) * ((fade / 384.0) / pow(length(star_pos - p), 1.5)) * fade;
  }
  gl_FragColor = vec4( gradient, 1.0);
}

### 37045:S ---

// dots1 ]
// lines2 [

### 37046:S lines2.fsh

#define time CC_Time[3]
#define resolution vec2(1024.0/2., 768.0/2.)
#define mouse vec2(0.5)

// By: Brandon Fogerty
// bfogerty at gmail dot com

#ifdef GL_ES
precision highp float;
#endif

varying mediump vec2 v_texCoord;
#define gl_FragCoord  v_texCoord*vec2(1000.)
//uniform float time;
//uniform vec2 resolution;

void main( void ) 
{
    vec2 p = ( gl_FragCoord.xy / resolution.xy ) * 2.0 - 1.0;
    
    vec3 c = vec3( 0.0 );
    
    float amplitude = 3.05; 
    float glowT = sin(time) * 0.6 + 0.1;
    float glowFactor = mix( 0.15, 0.35, glowT );
    
    c += vec3(0.02, 0.03, 0.13) * ( glowFactor * abs( 1.0 / sin(p.x + sin( p.y + time ) * amplitude ) ));
    c += vec3(0.02, 0.10, 0.03) * ( glowFactor * abs( 1.0 / sin(p.x + cos( p.y + time+1.00 ) * amplitude+0.1 ) ));
    c += vec3(0.15, 0.05, 0.20) * ( glowFactor * abs( 1.0 / sin(p.y + sin( p.x + time+1.30 ) * amplitude+0.15 ) ));
    c += vec3(0.20, 0.05, 0.05) * ( glowFactor * abs( 1.0 / sin(p.y + cos( p.x + time+3.00 ) * amplitude+0.3 ) ));
    c += vec3(0.17, 0.17, 0.05) * ( glowFactor * abs( 1.0 / sin(p.y + cos( p.x + time+5.00 ) * amplitude+0.2 ) ));
    
    gl_FragColor = vec4( c, 400.0 );

}

### 37047:S ---

// lines2 ]

// noise-SimplexCellular2D [

### 37048:S lines2.fsh

#ifdef GL_ES
precision highp float;
#endif


void FAST32_hash_2D( vec2 gridcell, out vec4 hash_0, out vec4 hash_1 )  //  generates 2 random numbers for each of the 4 cell corners
{
    //    gridcell is assumed to be an integer coordinate
    const vec2 OFFSET = vec2( 26.0, 161.0 );
    const float DOMAIN = 71.0;
    const vec2 SOMELARGEFLOATS = vec2( 951.135664, 642.949883 );
    vec4 P = vec4( gridcell.xy, gridcell.xy + 1.0 );
    P = P - floor(P * ( 1.0 / DOMAIN )) * DOMAIN;
    P += OFFSET.xyxy;
    P *= P;
    P = P.xzxz * P.yyww;
    hash_0 = fract( P * ( 1.0 / SOMELARGEFLOATS.x ) );
    hash_1 = fract( P * ( 1.0 / SOMELARGEFLOATS.y ) );
}

//  convert a 0.0->1.0 sample to a -1.0->1.0 sample weighted towards the extremes
vec4 Cellular_weight_samples( vec4 samples )
{
    samples = samples * 2.0 - 1.0;
    //return (1.0 - samples * samples) * sign(samples); // square
    return (samples * samples * samples) - sign(samples);   // cubic (even more variance)
}

//
//  SimplexCellular2D
//  cellular noise over a simplex (triangular) grid
//  Return value range of 0.0->~1.0
//  http://briansharpe.files.wordpress.com/2012/01/simplexcellularsample.jpg
//
//  TODO:  scaling of return value to strict 0.0->1.0 range
//
float SimplexCellular2D( vec2 P )
{
    //  simplex math based off Stefan Gustavson's and Ian McEwan's work at...
    //  http://github.com/ashima/webgl-noise

    //  simplex math constants
    const float SKEWFACTOR = 0.36602540378443864676372317075294;            // 0.5*(sqrt(3.0)-1.0)
    const float UNSKEWFACTOR = 0.21132486540518711774542560974902;          // (3.0-sqrt(3.0))/6.0
    const float SIMPLEX_TRI_HEIGHT = 0.70710678118654752440084436210485;    // sqrt( 0.5 )  height of simplex triangle.
    const float INV_SIMPLEX_TRI_HEIGHT = 1.4142135623730950488016887242097; //  1.0 / sqrt( 0.5 )
    const vec3 SIMPLEX_POINTS = vec3( 1.0-UNSKEWFACTOR, -UNSKEWFACTOR, 1.0-2.0*UNSKEWFACTOR ) * INV_SIMPLEX_TRI_HEIGHT;     //  vertex info for simplex triangle

    //  establish our grid cell.
    P *= SIMPLEX_TRI_HEIGHT;        // scale space so we can have an approx feature size of 1.0  ( optional )
    vec2 Pi = floor( P + dot( P, vec2( SKEWFACTOR ) ) );

    //  calculate the hash.
    //  ( various hashing methods listed in order of speed )
    vec4 hash_x, hash_y;
    FAST32_hash_2D( Pi, hash_x, hash_y );
    //SGPP_hash_2D( Pi, hash_x, hash_y );

    //  push hash values to extremes of jitter window
    const float JITTER_WINDOW = ( 0.10566243270259355887271280487451 * INV_SIMPLEX_TRI_HEIGHT );        // this will guarentee no artifacts.
    hash_x = Cellular_weight_samples( hash_x ) * JITTER_WINDOW;
    hash_y = Cellular_weight_samples( hash_y ) * JITTER_WINDOW;

    //  calculate sq distance to closest point
    vec2 p0 = ( ( Pi - dot( Pi, vec2( UNSKEWFACTOR ) ) ) - P ) * INV_SIMPLEX_TRI_HEIGHT;
    hash_x += p0.xxxx;
    hash_y += p0.yyyy;
    hash_x.yzw += SIMPLEX_POINTS.xyz;
    hash_y.yzw += SIMPLEX_POINTS.yxz;
    vec4 distsq = hash_x*hash_x + hash_y*hash_y;
    vec2 tmp = min( distsq.xy, distsq.zw );
    return min( tmp.x, tmp.y );
}

float SimplexPerlin2D( vec2 P )
{
    //  simplex math constants
    const float SKEWFACTOR = 0.36602540378443864676372317075294;            // 0.5*(sqrt(3.0)-1.0)
    const float UNSKEWFACTOR = 0.21132486540518711774542560974902;          // (3.0-sqrt(3.0))/6.0
    const float SIMPLEX_TRI_HEIGHT = 0.70710678118654752440084436210485;    // sqrt( 0.5 )  height of simplex triangle
    const vec3 SIMPLEX_POINTS = vec3( 1.0-UNSKEWFACTOR, -UNSKEWFACTOR, 1.0-2.0*UNSKEWFACTOR );      //  vertex info for simplex triangle

    //  establish our grid cell.
    P *= SIMPLEX_TRI_HEIGHT;        // scale space so we can have an approx feature size of 1.0  ( optional )
    vec2 Pi = floor( P + dot( P, vec2( SKEWFACTOR ) ) );

    //  calculate the hash.
    //  ( various hashing methods listed in order of speed )
    vec4 hash_x, hash_y;
    FAST32_hash_2D( Pi, hash_x, hash_y );
    //SGPP_hash_2D( Pi, hash_x, hash_y );

    //  establish vectors to the 3 corners of our simplex triangle
    vec2 v0 = Pi - dot( Pi, vec2( UNSKEWFACTOR ) ) - P;
    vec4 v1pos_v1hash = (v0.x < v0.y) ? vec4(SIMPLEX_POINTS.xy, hash_x.y, hash_y.y) : vec4(SIMPLEX_POINTS.yx, hash_x.z, hash_y.z);
    vec4 v12 = vec4( v1pos_v1hash.xy, SIMPLEX_POINTS.zz ) + v0.xyxy;

    //  calculate the dotproduct of our 3 corner vectors with 3 random normalized vectors
    vec3 grad_x = vec3( hash_x.x, v1pos_v1hash.z, hash_x.w ) - 0.49999;
    vec3 grad_y = vec3( hash_y.x, v1pos_v1hash.w, hash_y.w ) - 0.49999;
    vec3 grad_results = inversesqrt( grad_x * grad_x + grad_y * grad_y ) * ( grad_x * vec3( v0.x, v12.xz ) + grad_y * vec3( v0.y, v12.yw ) );

    //  Normalization factor to scale the final result to a strict 1.0->-1.0 range
    //  x = ( sqrt( 0.5 )/sqrt( 0.75 ) ) * 0.5
    //  NF = 1.0 / ( x * ( ( 0.5 – x*x ) ^ 4 ) * 2.0 )
    //  http://briansharpe.wordpress.com/2012/01/13/simplex-noise/#comment-36
    const float FINAL_NORMALIZATION = 99.204334582718712976990005025589;

    //  evaluate the surflet, sum and return
    vec3 m = vec3( v0.x, v12.xz ) * vec3( v0.x, v12.xz ) + vec3( v0.y, v12.yw ) * vec3( v0.y, v12.yw );
    m = max(0.5 - m, 0.0);      //  The 0.5 here is SIMPLEX_TRI_HEIGHT^2
    m = m*m;
    return dot(m*m, grad_results) * FINAL_NORMALIZATION;
}


#define time CC_Time[3]
#define resolution vec2(200.0/2., 200.0/2.)
#define mouse vec2(0.5)

varying mediump vec2 v_texCoord;
//#define gl_FragCoord  v_texCoord*vec2(1000.)

void main(void)
{
    float px = v_texCoord.x;
    float py = v_texCoord.y;
    float r = 0.;//cos(py*3.14*100.);
    float g = 0.;//rand(v_texCoord);
    float b = sin(px*3.14*100.)*0.5+0.5;
    #define time 0.
//    gl_FragColor = vec4(r,g,b,0.);
//    r=g= SimplexCellular2D(vec2(px*20.-time, py*20.))*2.;
//    r+= SimplexCellular2D(vec2(px*2.+10.+time, py*2.+10.))*2.;
//    g+= SimplexCellular2D(vec2(px*2.+10.+time, py*2.+10.))*2.;
    r=g= SimplexPerlin2D(vec2(px*10., py*10.))*0.5;
    r*= SimplexPerlin2D(vec2(px*2.+10.+time, py*2.+10.))*2.;
    g+= SimplexPerlin2D(vec2(px*2.+10.+time, py*2.+10.))*2.;
    b=r*g;
    vec3 v=clamp(vec3(r,g,b), vec3(0.), vec3(1.));
//    v=vec3(1.);

///    v=vec3((v.r+v.g+v.b)/3.);
//    r=g=b = rand(v_texCoord);
//    gl_FragColor = vec4(r,g,b,1.);
    gl_FragColor = vec4(v,1.);
}


### 37049:S ---

// noise-SimplexCellular2D ]

