# Fire [
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

# Fire ]

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
# Water [
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
    gl_FragColor = vec4(pow(color,vec3(0.75)), 1.0);
}
# Water ]
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
# Hellfire [
// Created by inigo quilez - iq/2013
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.

float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
    
    vec2 uv = (p.xy+vec2(37.0,17.0)*p.z) + f.xy;
    vec2 rg = texture2D( iChannel0, (uv+ 0.5)/256.0, -100.0 ).yx;
    return mix( rg.x, rg.y, f.z );
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
    t += 0.05*texture2D( iChannel0, gl_FragCoord.xy/iChannelResolution[0].x ).x;
    
    for( int i=0; i<100; i++ )
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
    
    vec2 mo = iMouse.xy / iResolution.xy;
    if( iMouse.w<=0.00001 ) mo=vec2(0.0);
    
    // camera
    vec3 ro = 4.0*normalize(vec3(cos(3.0*mo.x), 1.4 - 1.0*(mo.y-.1), sin(3.0*mo.x)));
    vec3 ta = vec3(0.0, 1.0, 0.0);
    float cr = 0.5*cos(0.7*iGlobalTime);
    
    // shake        
    ro += 0.1*(-1.0+2.0*texture2D( iChannel0, iGlobalTime*vec2(0.010,0.014) ).xyz);
    ta += 0.1*(-1.0+2.0*texture2D( iChannel0, iGlobalTime*vec2(0.013,0.008) ).xyz);
    
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

# Hellfire ]

### 37000 shader-test

r0:=k==0;

draw(r0, 37001);
x:=-100;
y:=200;
k:=1;

### 37001 shader-test.draw

color(u0, 0.1,0.1,0.1,0.1);
rect(u0, 6, 0,0, 200,200);

_ @1010;
_ @sys.node;
_ @sys.display;     

node.ex.shader.use(u0, 37100, 37101);
node.ex.gl.rect(u0, 0,0, 200,200);

### 37002 ---

// Waves [

### 37100:S Waves.vsh

attribute vec4 a_position;
//attribute vec2 a_texCoord;

//varying mediump vec2 v_texCoord;

void main()
{
    gl_Position = CC_MVPMatrix * a_position;
//    v_texCoord = a_texCoord;
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
    // vec2 uv = gl_FragCoord.xy / iResolution.xy;
    
    // // To create the BG pattern
    // vec3 final_color = vec3(1.0);
    // vec3 bg_color = vec3(0.0);
    // vec3 wave_color = vec3(0.0);
    
    // float c1 = mod(uv.x, 2.0 * BLOCK_WIDTH);
    // c1 = step(BLOCK_WIDTH, c1);
    
    // float c2 = mod(uv.y, 2.0 * BLOCK_WIDTH);
    // c2 = step(BLOCK_WIDTH, c2);
    
    // bg_color = mix(uv.x * COLOR1, uv.y * COLOR2, c1 * c2);
    
    
    // // To create the waves
    // float wave_width = 0.01;
    // uv  = -1.0 + 2.0 * uv;
    // uv.y += 0.1;
    // for(float i = 0.0; i < 10.0; i++) {
        
    //     uv.y += (0.07 * sin(uv.x + i/7.0 + iGlobalTime ));
    //     wave_width = abs(1.0 / (150.0 * uv.y));
    //     wave_color += vec3(wave_width * 1.9, wave_width, wave_width * 1.5);
    // }
    
    // final_color = bg_color + wave_color;
    
    
    // gl_FragColor = vec4(final_color, 1.0);
    gl_FragColor = vec4(0.0,1.0,0.0,1.0);
}

// Waves ]
