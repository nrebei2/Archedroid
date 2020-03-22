#include "./Math.cginc"


float2 fractalTexMap(float3 p, float3 txt, float fct1, float fct2, float fct3, float4 TexParam)
{
    //x = patternScale_
    //y = colmapYscale_
    //z = colmapPattern_
    //w = colmapFractal_
    
    #ifdef FCT_ORBIT
         return float2(txt.r *0.1 + fct1 * TexParam.w + p.y * TexParam.y,pow(fct2,0.001)* 20 + txt.g * TexParam.z);
    #elif defined (FCT_BBSK)
        //return vec2(pow(fct1,0.05)*2.,TexParam.w+p.y * TexParam.y + txt.r * TexParam.z);
        return float2(txt.g*0.1+pow(fct1,0.01)*10.-fct3 * 3.,fct2*TexParam.w+p.y * TexParam.y + txt.r * TexParam.z);
    #elif defined (FCT_DLBT)
         return float2(0.5+0.5*sin(txt.g * 1.2 + fct2 * TexParam.w + fct3),TexParam.y*(pow(fct1,0.01)/0.2+p.y / 22.5 + txt.r * TexParam.z) - fct3 * 0.3);
    #elif defined (FCT_MZGN)
        return float2(pow(fct1,0.05)*2.+txt.g*0.1 + fct3,fct2 * TexParam.w + p.z * TexParam.y + txt.r * TexParam.z + fct3*0.5);
    #elif defined (FCT_PIPES)
         return float2(txt.r *0.1 + fct1 * TexParam.w + p.y * TexParam.y ,pow(fct2,0.001)* 20 + txt.g * TexParam.z  + fct3 * txt.g);
    #elif defined (FCT_HTVT)
        return float2(txt.g*0.1+pow(fct1,0.005)*20. +fct3,fct2*TexParam.w+p.y * TexParam.y + txt.r * TexParam.z - fct3*0.3);
    #else
        return float2(txt.g + fct2*TexParam.z + fct3,fct1*TexParam.w+p.y * TexParam.y + txt.r * TexParam.z + fct3 * 0.3);
        //return vec2(1.,fct3);
    #endif
}

// Sphere
// s: radius
float sdSphere(float3 p, float s) {
	return length(p) - s;
}

// Box
inline float Box(float3 pos, float3 size)
{
    float3 d = abs(pos) - size;
    return length(max(abs(pos) - size, 0.0))
        + min(max(d.x, max(d.y, d.z)), 0.0);
}

// Rounded Box
float sdRoundBox(in float3 p, in float3 b, in float r) {
	float3 q = abs(p) - b;
	return min(max(q.x, max(q.y, q.z)), 0.0) + length(max(q, 0.0)) - r;
}

inline float Cylinder(float3 pos, float2 r)
{
    float2 d = abs(float2(length(pos.xy), pos.z)) - r;
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0)) - 0.1;
}

inline float HexagonalPrismX(float3 pos, float2 h)
{
    float3 p = abs(pos);
    return max(
        p.x - h.y, 
        max(
            (p.z * 0.866025 + p.y * 0.5),
            p.y
        ) - h.x
    );
}

inline float HexagonalPrismY(float3 pos, float2 h)
{
    float3 p = abs(pos);
    return max(
        p.y - h.y, 
        max(
            (p.z * 0.866025 + p.x * 0.5),
            p.x
        ) - h.x
    );
}

inline float HexagonalPrismZ(float3 pos, float2 h)
{
    float3 p = abs(pos);
    return max(
        p.z - h.y, 
        max(
            (p.x * 0.866025 + p.y * 0.5),
            p.y
        ) - h.x
    );
}

// Torus
float sdTorus(float3 p, float2 t) {
	float2 q = float2(length(p.xz)-t.x, p.y);
	return length(q) - t.y;
}


float4 printBulb(float4 d1, float k) {
	// Lerp between color and distance between objects
	float3 color = d1.rgb;
	float dist =  d1.w;
	return float4(color, dist);
}

// Subtraction
float opSS(float d1, float d2, float k) {
	float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
	return lerp(d2, -d1, h) + k * h * (1.0 - h);
}

// Intersection
float opIS(float d1, float d2, float k) {
	float h = clamp(0.5 - 0.5 * (d2 - d1) / k, 0.0, 1.0);
	return lerp(d2, d1, h) + k * h * (1.0 - h);
}

// trinoise
float tri(float x)
{
    return abs(frac(x) - .5);
}

float3 tri3(float3 p)
{
    return float3(
        tri(p.z + tri(p.y * 1.)),
        tri(p.z + tri(p.x * 1.)),
        tri(p.y + tri(p.x * 1.)) );
}

float trinoise(float3 p, float spd, float time)
{
    float z = 1.4;
    float rz = 0.;
    float3  bp = p;
    for (float i = 0.; i <= 3.; i++) {
        float3 dg = tri3(bp * 2.);
        p += (dg + time * .1 * spd);
        bp *= 1.8;
        z *= 1.5;
        p *= 1.2;
        float t = tri(p.z + tri(p.x + tri(p.y)));
        rz += t / z;
        bp += 0.14;
    }
    return rz;
}


// original code: https://www.shadertoy.com/view/ldfGWn
inline float truchetarc(float3 pos)
{
    float p = 4.0 + 2.0 * _SinTime.w;
    float r = length(pos.xy);
    float t = 0.12 + 0.02 * _SinTime.w;
    return pow(
        pow(abs(r - 0.5), p) + pow(abs(pos.z - 0.5), p), 
        rcp(p)
    ) - t;
}
inline float truchetcell(float3 pos)
{
    return min(min(
        truchetarc(pos),
        truchetarc(float3(pos.z, 1.0 - pos.x, pos.y))), 
        truchetarc(float3(1.0 - pos.y, 1.0 - pos.z, pos.x)));
}


float4 sdfmap(float3 pos, float power, float3 Params){
    float dist = 1000.;
    float hdist = -1000000.;
    float tex = 0;
    float tex2 = 0;
    float tex3 = 0;
    float3 cFcParams = Params;
  
    #ifdef StaticMandelbulb
        float3 w = pos;
        float m = dot(w, w);
    
        float dz = 1.0;
            
        for(int i = 0; i < 2; i++)
        {
            dz = 8 * pow(sqrt(m), 7.0)*dz + 1.0;
            float r = length(w);
            float b = 8 * acos(w.y / r);
            float a = 8 * atan2(w.x, w.z);
            w = pos + pow(r, 8) * float3(sin(b) * sin(a), cos(b), sin(b) * cos(a));
    
            m = dot(w, w);
            if(m > 256.0)
                break;
        }
        dist = 0.25*log(m)*sqrt(m)/dz;
    #elif defined (DynamicMandelbulb)
        float sinTime = sin(_Time.y / 1);
        float dpower = remap(sinTime, -1, 1, 4, 9);
        float3 z = pos;
        float r = 0;
        float dr = 1;
        for(int i = 0; i < 5; i++) 
        {
            r = length(z);
            if(r > 100) break;
            
            float theta = acos(z.z / r);
            float phi = atan2(z.y, z.x);
            
            dr = dpower * pow(r, dpower-1)*dr+1;
            
            r = pow(r, dpower);
            theta *= dpower;
            phi *= dpower;
            
            z = r * float3(sin(theta) * cos(phi), 
                    sin(theta) * sin(phi), 
                    cos(theta));
    
            z += pos;
        }
        dist = 0.5 * log(r) * r / dr;
    #elif defined (Julia)
        float4 c = 0.45* cos( float4(0.5,3.9,1.4,1.1) + _Time.y * float4(1.2,1.7,1.3,2.5) ) - float4(0.3,0.0,0.0,0.0);
        float4 z = float4(pos, 0);
        float md2 = 1;
        float mz2 = dot(z, z);
    
        [loop]
        for(int i = 0; i < 11; i++)
        {
            md2 *= 4.0 * mz2; // dz -> 2·z·dz, meaning |dz| -> 2·|z|·|dz| (can take the 4 out of the loop and do an exp2() afterwards)
            z = qsqr(z) + c; // z  -> z^2 + c
    
            mz2 = dot(z,z);
    
            if(mz2 > 4.0) break;
        }
        
        dist = 0.25 * sqrt(mz2/md2) * log(mz2);
    #elif defined (Juliabulb)
        float4 c = 0.45* cos( float4(0.5,3.9,1.4,1.1) + _Time.y * float4(1.2,1.7,1.3,2.5) ) - float4(0.3,0.0,0.0,0.0);
        float3 orbit = pos;
        float dz = 1;
        
        for (int i = 0; i < 4; i++) 
        {
            float r = length(orbit);
            float o = acos(orbit.z/r);
            float p = atan(orbit.y/orbit.x);
            
            dz = 8*r*r*r*r*r*r*r*dz;
            
            r = r*r*r*r*r*r*r*r;
            o = 8*o;
            p = 8*p;
            
            orbit = float3(r*sin(o) * cos(p), 
                    r*sin(o) * sin(p), 
                    r*cos(o)) + c;
            
            if (dot(orbit, orbit) > 4.0) break;
        }
        float z = length(orbit);
        dist = 0.5*z*log(z)/dz;
    #elif defined (Sierpinski)
        const float3 va = float3(  0.0,  0.575735,  0.0 );
        const float3 vb = float3(  0.0, -1.0,  1.15470 );
        const float3 vc = float3(  1.0, -1.0, -0.57735 );
        const float3 vd = float3( -1.0, -1.0, -0.57735 );
    
        float a = 0;
        float s = 1;
        float r = 1;
        float dm;
        float3 v;
        [loop]
        for(int i = 0; i < 15; i++)
        {
            float d, t;
            d = dot(pos - va, pos - va);
    
            v = va; 
            dm = d; 
            t = 0;
            
            d = dot(pos - vb, pos - vb); 
            if(d < dm) 
            { 
                v = vb; 
                dm=d; 
                t = 1.0; 
            }
            
            d = dot(pos-vc, pos-vc); 
    
            if(d < dm) { v = vc; dm = d; t = 2.0; }
            d = dot(pos-vd,pos-vd); 
            if(d < dm) { v = vd; dm = d; t = 3.0; }
    
            pos = v + 2*(pos - v); 
            r*= 2;
            a = t + 4*a; 
            s*= 4;
        }
        
        dist = float2((sqrt(dm)-1.0)/r, a/s);
    #elif defined (Mandelbox)
        float SCALE = 2.75;
        float fixedRadius = 1.0;
        float FR2 = fixedRadius * fixedRadius;
        float minRadius = 0.5;
        float MR2 = minRadius * minRadius;
        float4 scalefloat = float4(SCALE, SCALE, SCALE, abs(SCALE)) / MR2;
        float C1 = abs(SCALE-1.0);
        float C2 = pow(abs(SCALE), float(1-5));
        float4 p = float4(pos.xyz, 1.0); 
        float4 p0 = float4(pos.xyz, 1.0);  // p.w is knighty's DEfactor
        for (int i=0; i<5; i++) {
            p.xyz = clamp(p.xyz *0.5+0.5, 0.0, 1.0) *4.0-2.0 - p.xyz; // box fold: min3, max3, mad3
            float r2 = dot(p.xyz, p.xyz);  // dp3
            p.xyzw *= clamp(max(MR2/r2, MR2), 0.0, 1.0);  // sphere fold: div1, max1.sat, mul4
            p.xyzw = p*scalefloat + p0;  // mad4
        }
      return (length(p.xyz) - C1) / p.w - C2;
    #elif defined (KaleidoscopicIFS)
        int FRACT_ITER      = 20;
        float FRACT_SCALE   = 1.8;
        float FRACT_OFFSET  = 1.0;
    
        float c = 2.0;
        pos.y = modc(pos.y, c)-c/2.0;
        pos = RotateZ(pos, (PI/2.0)/0.0174532925);
        float r;
        int n1 = 0;
        for (int n = 0; n < FRACT_ITER; n++) {
            float rotate = PI*0.5;
            pos = RotateX(pos, rotate/0.0174532925);
            pos = RotateY(pos, rotate/0.0174532925);
            pos = RotateZ(pos, rotate/0.0174532925);
    
            pos.xy = abs(pos.xy);
            if (pos.x+pos.y<0.0) pos.xy = -pos.yx; // fold 1
            if (pos.x+pos.z<0.0) pos.xz = -pos.zx; // fold 2
            if (pos.y+pos.z<0.0) pos.zy = -pos.yz; // fold 3
            pos = pos*FRACT_SCALE - FRACT_OFFSET*(FRACT_SCALE-1.0);
        }
        dist = (length(pos) ) * pow(FRACT_SCALE, -float(FRACT_ITER));
    #elif defined (Tglad)
        pos = modc(pos, 2.0);
    
        float mr=0.25, mxr=1.0;
        float4 scale=float4(-3.12,-3.12,-3.12,3.12), p0=float4(0.0,1.59,-1.0,0.0);
        float4 z = float4(pos,1.0);
        for (int n = 0; n < 3; n++) {
            z.xyz=clamp(z.xyz, -0.94, 0.94)*2.0-z.xyz;
            z*=scale/clamp(dot(z.xyz,z.xyz),mr,mxr);
            z+=p0;
        }
        dist =(length(max(abs(z.xyz)-float3(1.2,49.0,1.4),0.0))-0.06)/z.w;
    #elif defined (Hartverdrahtet)
        // distance function from Hartverdrahtet
        // ( http://www.pouet.net/prod.php?which=59086 )
        float3 cs=float3(.808,.808,1.167);
        float fs=1.;
        float3 fc=0;
        float fu=10.;
        float fd=.763;
        
        // scene selection
        {
            float time = _Time.y;
            int i = int(modc(time/2.0, 9.0));
            if(i==0) cs.y=.58;
            if(i==1) cs.xy=.5;
            if(i==2) cs.xy=.5;
            if(i==3) fu=1.01,cs.x=.9;
            if(i==4) fu=1.01,cs.x=.9;
            if(i==6) cs=float3(.5,.5,1.04);
            if(i==5) fu=.9;
            if(i==7) fd=.7,fs=1.34,cs.xy=.5;
            if(i==8) fc.z=-.38;
        }
        
        //cs += sin(time)*0.2;
    
        float v=1.;
        for(int i=0; i<12; i++){
            pos=2.*clamp(pos,-cs,cs)-pos;
            float c=max(fs/dot(pos,pos),1.);
            pos*=c;
            v*=c;
            pos+=fc;
        }
        float z=length(pos.xy)-fu;
        return fd*max(z,abs(length(pos.xy)*pos.z)/sqrt(dot(pos,pos)))/abs(v);
    #elif defined (PseudoKleinian)
        float3 CSize = float3(0.92436,0.90756,0.92436);
        float Size = 1.0;
        float3 C = float3(0.0,0.0,0.0);
        float DEfactor=1.;
        float3 Offset = float3(0.0,0.0,0.0);
        float3 ap=pos+1.;
        for(int i=0;i<10 ;i++){
            ap=pos;
            pos=2.*clamp(pos, -CSize, CSize)-pos;
            float r2 = dot(pos,pos);
            float k = max(Size/r2,1.);
            pos *= k;
            DEfactor *= k + 0.05;
            pos += C;
        }
        dist = abs(0.5*abs(pos.z-Offset.z)/DEfactor);
    #elif defined (PseudoKnightyan)
        float3 CSize = float3(0.63248,0.78632,0.875);
        float DEfactor=1.;
        for(int i=0;i<6;i++){
            pos = 2.*clamp(pos, -CSize, CSize)-pos;
            float k = max(0.70968/dot(pos,pos),1.);
            pos *= k;
            DEfactor *= k + 0.05;
        }
        float rxy=length(pos.xy);
        dist =  max(rxy-0.92784, abs(rxy*pos.z) / length(pos))/DEfactor;
    #elif defined (Mandelbulb2)
        float3 z = pos;
        float dr = 1.0;
        float r = 0.0;
        int iterations = 0;
    
        for (int i = 0; i < 100 ; i++) {
            iterations = i;
            r = length(z);
    
            if (r>3) {
                break;
            }
            
            // convert to polar coordinates
            float theta = asin( z.z/r );
            float phi = atan2( z.y,z.x );
            dr =  pow( r, power-1.0)*power*dr + 1.0;
    
            // scale and rotate the point
            float zr = pow( r,power);
            theta = theta*power;
            phi = phi*power;
            
            // convert back to cartesian coordinates
            z = zr*float3( cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta) );
            z+=pos;
        }
        dist = 0.5*log(r)*r/dr;
    #elif defined (MengerSponge)
        // MengerSponge distance estimation:
        // http://www.iquilezles.org/www/articles/menger/menger.htm
        float d = Box(pos, 50.0); // size of box
    
        float s = 0.02; // size of squares in box
        int iterations = 0;
       
        for( int m=0; m<5; m++ )
        {
          //p = p + float3(power-1.0,power-1.0,power-1.0);
          iterations = m;
          float3 a = (pos*s - 2.0 * floor(pos*s/2.0))-1.0;
          s *= 3.0;
          float3 r = abs(1.0 - 3.0*abs(a));
    
          float da = max(r.x,r.y);
          float db = max(r.y,r.z);
          float dc = max(r.z,r.x);
          float c = (min(da,min(db,dc))-1.0)/s;
    
          d = max(d,c);
        }
    
        dist = d;
    #elif defined (apo)
          float seed = .0274;
          float3 CSize = float3(1., 1., 1.3);
          float3 C = float3(0., 0., 0.);
          //float3 CSize = float3(1., 1., 1.3);
          float3 p = pos.xzy;
          float scale = 1.0;
         // p *= cFcRot;
          float r2 = 0.;
          float k = 0.;
          //float uggg = 0.;
          for( int i=0; i < 12;i++ )
          {
              p = 2.0*clamp(p, -CSize, CSize) - p;
              r2 = dot(p,p);
              //r2 = dot(p,p+sin(p.z*.3)); //Alternate fractal
              k = max((2.0)/(r2), seed); //.378888 //.13345 max((2.6)/(r2), .03211); //max((1.8)/(r2), .0018);
              p     *= k;
              scale *= k;
              //uggg += r2;
              p+=C;
               //p.xyz = float3(-1.0*p.z,1.0*p.x,1.0*p.y);
          }
          float l = length(p.xy);
          float rxy = l - 4.0;
          float n = 1.0 * p.z;
          rxy = max(rxy, -(n) / 4.);
          dist = (rxy) / abs(scale);
    #elif defined (plane)
        // (Infinite) Plane
        // n.xyz: normal of the plane (normalized)
        // n.w: offset
        float4 n = float4(0,1,0,0);
        dist = dot(pos, n.xyz) + n.w;
    #elif defined (FCT_BBSK)
        // Notable Params: (2.18 -0.18 0)

        float3 CSize = float3(1.4,0.87, 1.1);
        float3 p = pos.xzy * 2.0;
        float scale = 1.0;
        
        for( int i=0; i < 4;i++ )
        {
            p = 2.0*clamp(p, -CSize, CSize) - p;
            //float r2 = dot(p,p);
            float r2 = dot(p,p+sin(p.z*.5)); //Alternate fractal
            float k = max((2.)/(r2), .17);
            p *= k;
            //p *=rot;
            //p= p.yzx;
            p+=float3(0.2,0.2,-0.5);
            scale *= k;
        }
        tex2 = p.x;
        
        p = 2.0*clamp(p, -CSize * 4., CSize * 4.) - p;
       
        for(int i=0; i < 8; i++ )
        {
            p = 2.0*clamp(p, -CSize, CSize) - p;
            float r2 = dot(p,p);
            //float r2 = dot(p,p+sin(p.z*.3)); //Alternate fractal
            float k = max((cFcParams.x)/(r2),  0.027);
            p     *= k;
            scale *= k;
            p.y += cFcParams.y;
        }
        
        float l = length(p.xy);
        //l = mix(l,l2,0.5);
        float rxy = l - 4.0;
        float n = 1.0 * p.z;
        rxy = max(rxy, -(n) / 4.);
        dist = (rxy) / abs(scale);
        dist *=.75;
        tex = min(scale,150.);
    #elif defined (trinoise)
        dist = trinoise(pos, 1.0, 1.0);
    #elif defined (RecursiveTetrahedron)
        pos = Repeat(pos / 2, 3.0);

        const float3 a1 = float3( 1.0,  1.0,  1.0);
        const float3 a2 = float3(-1.0, -1.0,  1.0);
        const float3 a3 = float3( 1.0, -1.0, -1.0);
        const float3 a4 = float3(-1.0,  1.0, -1.0);
    
        const float scale = 2.0;
        float d;
        for (int n = 0; n < 3; ++n) {
            float3 c = a1; 
            float minDist = length(pos - a1);
            d = length(pos - a2); if (d < minDist) { c = a2; minDist = d; }
            d = length(pos - a3); if (d < minDist) { c = a3; minDist = d; }
            d = length(pos - a4); if (d < minDist) { c = a4; minDist = d; }
            pos = scale * pos - c * (scale - 1.0);
        }
     
        dist = length(pos) * pow(scale, float(-n));
    #elif defined (TruchetTentacles)
        float3 c = frac(pos);
        float  r = rand(floor(pos));
    
        if      (r < 0.125) dist = truchetcell(float3(c.x, c.y, c.z));
        else if (r < 0.250) dist = truchetcell(float3(c.x, 1.0 - c.y, c.z));
        else if (r < 0.375) dist = truchetcell(float3(1.0 - c.x, c.y, c.z));
        else if (r < 0.500) dist = truchetcell(float3(1.0 - c.x, 1.0 - c.y, c.z));
        else if (r < 0.625) dist = truchetcell(float3(c.y, c.x, c.z));
        else if (r < 0.750) dist = truchetcell(float3(c.y, 1.0 - c.x, c.z));
        else if (r < 0.875) dist = truchetcell(float3(1.0 - c.y, c.x, c.z));
        else                dist = truchetcell(float3(1.0 - c.y, 1.0 - c.x, c.z));
    #elif defined (FCT_PROTEIN)
        // Notable Params: (-12.82 -0.63 -16.18)
        float3 p = pos * 0.002;
        float4 q = float4(p - 1., 1.);
        for(int i = 0; i < 7; i++) {
          q.xyz = abs(q.xyz + 1.3) - 1.0;
          q /= clamp(dot(q.xyz, q.xyz), 0.0, 0.8);
          //q.xyz *= cFcRot;
          q *= 1.567+cFcParams.x;// + p.y*0.8;
          q = q.zxyw;
          q+= float4(0.2119,0.8132,0.077213,0.);
        }
        for(int i = 0; i < 4; i++) {
          q.xyz = abs(q.xyz + 1.3) - 1.0;
          q /= clamp(dot(q.xyz, q.xyz), 0.0, 0.8);
          q *= 1.867;// + p.y*0.8;
        }
        dist = (length(q.xyz) - max(-240.-p.y*700.,2.5))/q.w * 300.;
    #elif defined (FCT_ORBIT)
        float3 p = pos * 0.001;
        float4 q = float4(p - 1.0, 1);
        for(int i = 0; i < 11; i++) {
          //q.xyz = mod(q.xyz,2.0)-0.5*q.xyz;
           q.xyz = abs(q.xyz + 1.0) - 1.0;
          q.xyz = 2.0*clamp(q.xyz, -73.0174, 73.0174) - q.xyz;
          q /= min(dot(q.xyz, q.xyz), 0.9);
          q *= 1.817;// + p.y*0.8;
          //q += float4(0.2,.02,-.2,0.2);
          //q.xyz *= rot;
        }
        dist = (length(q.xz) - 2.1)/q.w * 800.;
        tex = q.y;
        tex2 = q.w;
    #elif defined (FCT_MNMT)
        float3 CSize = float3(1., 1., 1.3); // <-- CSize Constant
        float3 p = pos.yzx; // <-- 3D Position
        float scale = 1.0; // <-- Scale
        
        for (int i = 0; i < 8; i++) { // <-- Primary iteration
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            float r2 = dot(p, p);
            float k = max((2.) / (r2), .1274);
            p *= k;
            //p *= cFcRot;                                    // \ Lines present only
            p.xyz = float3(1.0 * p.z, 1.0 * p.x, -1.0 * p.y); // / in this loop
            scale *= k;
        }
        
        CSize = float3(1.2, 0.4, 1.4); // <-- CSize Constant
        tex = p.y; //<-- Texture Primary
        
        for (int i = 0; i < 4; i++) { // <-- Secondary Iteration
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            float r2 = dot(p, p);
            float k = max((1.6) / (r2), 0.0274);
            p *= k;
            scale *= k;
        }
        
        float l = length(p.xyz);
        float rxy = l - 1.4;
        float n = 1.0 * p.z;
        rxy = max(rxy, -(n) / 4.);
        dist = (rxy) / abs(scale);
        
        dist *= 1.5;
    #elif defined (FCT_CRAB)
        float3 p = pos * 0.002; // <-- 3D Position
        float4 q = float4(p - 1.0, 1); // <-- 4D Map onto Position
        
        for (int i = 0; i < 7; i++) { // <-- Primary iteration
            q.xyz = abs(q.xyz + 1.3) - 1.0;
            q /= clamp(dot(q.xyz, q.xyz), 0.0, 0.8);
            q *= 1.867; // + p.y*0.8;
            q = q.zxyw;                              // \ Lines present only
            q += float4(0.2119, 0.8132, 0.077213, 0.); // / in this loop
        }
        
        for (int i = 0; i < 4; i++) { // <-- Secondary iteration
            q.xyz = abs(q.xyz + 1.3) - 1.0;
            q /= clamp(dot(q.xyz, q.xyz), 0.0, 0.8);
            q *= 1.867; // + p.y*0.8;
        }
        
        dist = (length(q.yz) - 1.2) / q.w * 300.; // <-- Distance
    #elif defined (FCT_HUB)
        // Notable Params: (0.514 -2.28 0.83), (0.644 -2.28 0.83)
        
        float3 p = pos.yzx * 0.05; // <-- 3D Position
        float4 q = float4(p - 1.0, 1); // <-- 4D Map onto Position
        
        for (int i = 0; i < 8; i++) { // <-- Primary iteration (All lines shared)
            q.xyz = abs(q.xyz + cFcParams.z) - 1.0;
            q /= clamp(dot(q.xyz, q.xyz), 0.12, 1.0);
            q *= 1.0 + cFcParams.x;
            //q.xyz *= cFcRot;
        }
        
        tex = q.x; // <-- Texture Primary
        
        for (int i = 0; i < 2; i++) { // <-- Secondary iteration
            q.xyz = abs(q.xyz + cFcParams.z) - 1.0;
            q /= clamp(dot(q.xyz, q.xyz), 0.12, 1.0);
            q *= 1.0 + cFcParams.x;
            //q.xyz *= cFcRot;
        }
        
        dist = (length(q.xyz) - 1. + cFcParams.y) / q.w * 19.; // <-- Distance
    #elif defined (FCT_HYPERAPO)
        // Notable Params: (0.644 -2.28 0.83), (0.067 1.05 -5.58), (0.76 0.91 4.22)
        float3 CSize = float3(1., 1., 1.1); // <-- CSize Constant
        float3 p = pos.yzx; // <-- 3D Position
        float scale = 1.0; // <-- Scale
        
        for (int i = 0; i < 4; i++) { // <-- Primary iteration
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            float r2 = dot(p, p);
            float k = max((2.) / (r2), 0.067);
            p *= k;
            p.xyz = float3(1.0 * p.z, 1.0 * p.x, -1.0 * p.y); // Line present only in this loop
            scale *= k;
        }
        
        p = p.zxy;
        //p *= cFcRot;
        
        CSize = float3(1.2, cFcParams.x, 1.2); // <-- CSize Constant
        tex = p.y; // <-- Texture Primary
        
        for (int i = 0; i < 7; i++) { // <-- Secondary iteration
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            float r2 = dot(p, p);
            float k = max((1.6) / (r2), cFcParams.y);
            p *= k;
            scale *= k;
        }
        
        float l = length(p.xyz);
        float rxy = l - 1.4;
        float n = 1.0 * p.z;
        rxy = max(rxy, -(n) / 4.);
        
        dist = (rxy) / abs(scale); // <-- Distance
    #elif defined (FCT_DLBT)
        // Notable Params: (-24.63 16.16 101.07)
        float3 npos = pos;
        float noise = noise3d(npos * 0.07) * 10.;
        float3 apopos = pos.xzy;
        apopos.z += cFcParams.x;
        domainRep3(apopos, float3(250., 250., 0));
        float r2 = 0.;
        float k = 0.;
        
        float3 CSize = float3(1., 1., 1.3); // <-- CSize Constant
        float3 p = apopos; // <-- 3D Position
        float scale = 1.0; // <-- Scale
        
        for (int i = 0; i < 5; i++) { // <-- Primary iteration (All lines shared)
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            r2 = dot(p, p);
            k = max((2.0) / (r2), .0274);
            p *= k;
            scale *= k;
        
        }
        
        tex = scale; // <-- Texture Primary
        
        for (int i = 0; i < 7; i++) { // <-- Secondary iteration
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            r2 = dot(p, p);
            k = max((2.0) / (r2), .0274);
            p *= k;
            scale *= k;
        }
        
        float l = length(p.xy);
        float rxy = l - 4.0;
        float n = 1.0 * p.z;
        rxy = max(rxy, -(n) / 4.);
        float apodist = (rxy) / abs(scale);
        apodist = 2. - apodist * 2.;
        dist = noise - apodist + 0.01 * (pos.y - 67.);
        dist *= 0.7;
        dist = smin(dist, pos.y + cFcParams.y, 15.);
        
        tex2 = p.z / scale; // <-- Texture Secondary
        dist = smin(dist, pos.y - cFcParams.z, -15.); // <-- Distance
    #elif defined (FCT_MZGN)
        // Notable Params: A lot
        float3 p = pos * 0.01; // <-- 3D Position
        float4 q = float4(p, 1); // <-- 4D Map onto Position
        
        float4 qd;
        
        for (int i = 0; i < 12; i++) { // <-- Primary iteration
            q.xyz = abs(q.xyz + 1.0 + cFcParams.y) - 1.0;
            qd = q;
            qd.w = qd.w * 0.002;
            q /= clamp(dot(qd, qd), 0.0, 0.8);
           // q.xyz *= cFcRot;
            q *= 1.567 + cFcParams.x;
        }
        
        tex = q.w; // <-- Texture Primary
        tex2 = q.x; // <-- Texture Secondary
        dist = (length(q.xyz) - 3.5) / q.w * 100.; // <-- Distance
    #elif defined (FCT_PIPES)
        float3 p = pos * 0.002; // <-- 3D Position
        float4 q = float4(p - 1.0, 1); // <-- 4D Map onto Position
        
        for (int i = 0; i < 11; i++) { // <-- Primary iteration
            q.xyz = abs(q.xyz + 1.0) - 1.0;
            q /= clamp(dot(q.xyz, q.xyz), 0.12, 1.0);
            q *= 1.837;
        }
        
        tex = q.y; // <-- Texture Primary
        tex2 = q.w; // <-- Texture Secondary
        dist = (length(q.xz) - 1.2) / q.w * 500.; // <-- Distance
    #elif defined (FCT_APOP)
        // Notable Params: (1.11 4.48 0.07)
        float3 p = pos * 0.01; // <-- 3D Position
        float scale = 1.0; // <-- Scale
        
        for (int i = 0; i < 10; i++) { // <-- Primary iteration
            p = -1.0 + 2.0 * fract(0.5 * p + 0.5);
            float r2 = dot(p, p);
            float k = cFcParams.x / r2;
            p *= k;
            scale *= k;
            //p *= cFcRot;
        }
        
        dist = length(p) / scale;
        
        dist *= 25.; // <-- Distance
    #elif defined (FCT_APO)
        float3 CSize = float3(1.7, 1.7, 1.3); // <-- CSize Constant
        float3 p = pos.xzy; // <-- 3D Position
        float scale = 1.; // <-- Scale
        
        for (int i = 0; i < 12; i++) { // <-- Primary iteration
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            float r2 = dot(p, p);
            float k = max((2.) / (r2), .027);
            p *= k;
            scale *= k;
        }
        
        float l = length(p.xyz);
        float rxy = l - 4.0;
        float n = l * p.z;
        rxy = max(rxy, -(n) / 4.);
        dist = (rxy) / abs(scale);
        
        dist = max(dist, pos.x); // <-- Distance
    #elif defined (FCT_HTVT)
        // Notable Params: (0.067 0.75 -0.02), (0.09 2.72 -0.51)
        float3 CSize = float3(1., 1., 1.3); // <-- CSize Constant
        float3 p = pos.yzx; // <-- 3D Position
        float scale = 1.0; // <-- Scale
        
        for (int i = 0; i < 4; i++) { // <-- Primary iteration
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            float r2 = dot(p, p);
            float k = max((2.) / (r2), cFcParams.x);
            p *= k;
            p.xyz = float3(1.0 * p.z, 1.0 * p.x, -1.0 * p.y); // Line present only in this loop
            scale *= k;
        }
        
        p = p.zxy;
        //p *= cFcRot;
        
        CSize = float3(1.2, 0.4, 1.4); // <-- CSize Constant
        tex2 = p.x; // <-- Texture Secondary
        
        for (int i = 0; i < 8; i++) { // <-- Secondary iteration
            p = 2.0 * clamp(p, -CSize, CSize) - p;
            float r2 = dot(p, p);
            float k = max((1.6) / (r2), cFcParams.y);
            p *= k;
            scale *= k;
        }
        
        float l = length(p.xyz);
        float rxy = l - 1.4;
        float n = 1.0 * p.z;
        rxy = max(rxy, -(n) / 4.);
        dist = (rxy) / abs(scale) - 0.0005;
        
        tex = min(scale, 150.); // <-- Texture Primary
        dist *= 1.4; // <-- Distance
    #elif defined (FCT_KNKL)
        // Knighty's Pseudo Kleinian Fractal
        // Notable Params: (-0.05 -0.37 0.07), (-0.84 0.97 0.13)
        float3 CSize = float3(0.97478, 1.4202, 0.97478); // <-- CSize Constant
        float3 p = pos.xzy * 0.01; // <-- 3D Position
        float3 C = cFcParams;
        
        float DEfactor = 1.;
        float3 ap = p + float3(1.,1.,1.);
        
        for(int i=0;i<11 && ap!=p;i++){
            ap=p;
            p=2.*clamp(p, -CSize, CSize)-p;  
            float r2=dot(p,p);
            float k=max(1.0/r2,1.);
            p*=k;DEfactor*=k;  
            p+=C;
        }
        
        tex = p.y; // <-- Texture Primary
        tex2 = DEfactor; // <-- Texture Secondary
        dist = (0.5 * (p.z) / DEfactor) * 70; // <-- Distance
    #elif defined (FCT_KIFS)
        // No notable Params, just experiment!
        float3 p = pos * 0.01; // <-- 3D Position
        
        float l = 0.0;
        float3 Fold = float3(3.,3.,3.);
        
        for (int i = 0; i < 12; i++) { // <-- Primary iteration
            p.xyz = abs(p.xyz + Fold.xyz) - Fold.xyz;
            p = p * cFcParams.x;
            //p *= cFcRot;
        }
        
        l = length(p);
        
        dist = (l * pow(cFcParams.x, -12.) - 0.001) * 100.; // <-- Distance
    #elif defined (FCT_TEXT)
        float3 p = pos; // <-- 3D Position

        float3 cell = domainRep3Idx(p, float3(12., 8., 12.));
        float3 h = float3(hash(cell.y), hash(cell.z + cell.x), hash(cell.x));
        float rot = h.x + h.y + h.z * 6.3;
        float2 sc = float2(sin(rot), cos(rot));
        p.xz = float2(sc.y * p.x - sc.x * p.z, sc.y * p.z + sc.x * p.x);
        float3 cube = abs(p) - float3(3.,3.,3.);
        cube *= 0.9;
        dist = length(max(cube, 0.0)) - 1.5;
        p += float3(5.,5.,5.);
        dist = max(dist, pos.y - 54.);
        if (dist < 2.) {
            //float4 vol = texture3D(sVolumeMap, clamp(p.zxy * 0.1, 0., 1.)) * 0.9;
            float d = 2.; //* (vol.r + vol.g + vol.b);
            dist = smax(dist - 1.5, d * 1.0 - 0.1, 2.);
        } else {
            dist += 1.;
        }
        
        dist = min(dist, pos.y + 4.); // <-- Distance
    #elif defined (FCT_TEST)
        float3 p = pos;
        p+=float3(0.,0.,14.);
        float r2 = dot(p,p);
        float k = max(20. / max(5., r2), 1.0);
        p*=k;
        dist = length(p) - 35.;
        dist/=k;
        tex = dist;
        dist = min(dist,pos.y);
    #else
        dist = 12000.;
    #endif
    
    return float4(dist, tex, tex2, tex3);
}
