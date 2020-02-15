using System;
using System.Collections;
using System.Collections.Generic;
using System.Numerics;
using UnityEngine;
using static Unity.Mathematics.math; // Used to implement shader like functions
using Unity.Mathematics; // Used to implement shader like functions 
using UnityEngine.UIElements;
using Vector2 = UnityEngine.Vector2;
using Vector3 = UnityEngine.Vector3;
using Vector4 = UnityEngine.Vector4;

// Use regex to find all decimals and add f at the end ot turn it to a float.
// Find: (\d\.\d+)    Remember that the parantheses create a group
// Replace: $1f       Remember to access the (first) group use $1

[ExecuteInEditMode]
public class RaymarchCollider : MonoBehaviour
{
    private float _fractalNumber;
    private Vector3 Params;
    private Vector3 _fractalPosition;
    private Vector3 _fractaldegreeRotate;
    private Vector3 _fractalScale;
    private float _power;

    public float distance;
    
    public Vector3 p;
    SphereCollider myCollider;
    Rigidbody rb;
    public Vector3 sdfNormal;
    public bool okToJump;
    
    Vector3 _gravityVec = new Vector3(0,-0.7071f,0.7071f);
    
    // TODO: do definitions of shader code and turn it to readable C# code, making it easy to translate code from shadercode to c#
    // TODO: implement collision similar to Yeomada
    // TODO: implement the sdf FCT_BBSK from distancefuction to see if something is really wrong with what im doing (try this first)
    private void Start()
    {
        myCollider = GetComponent<SphereCollider>();
        rb = GetComponent<Rigidbody>();
    }

    // Update is called once per frame
    void Update()
    {
        _fractalNumber = RaymarchingCamera._fractalNumber;
        Params = RaymarchingCamera.Params;
        _fractalPosition = RaymarchingCamera._fractalPosition;
        _fractaldegreeRotate = RaymarchingCamera._fractaldegreeRotate;
        _fractalScale = RaymarchingCamera._fractalScale;
        _power = RaymarchingCamera._power;
        
        p = GetComponent<Transform>().position;

        distance = GetDistance();

        okToJump = false;
        Vector3 playerV = rb.velocity;
        //sdfNormal = getNormal(position);

        // First of all, the merger sponge seems to be working perfectly. The object now knows that it isn't just a solid cube, which is nice.
        // The return value for the sdf still works such that:
        // 1. if the center of the object is inside the surface, sdf returns a negative number.
        // 2. if the center of the object is outside the surface, sdf returns a positive number.
        // 3. if the center of the object is EXACTLY on the surface, sdf returns a 0.
        // 4. sdf return value increases as object moves away from surface, and decreases as object moves further inside surface.
        // but, the scales for the return values are all wrong. Look at the sierpenski, for example.
        
        // Which means, there is a possible fix...
        // instead of using the position at the center of the object, make multiple positions around the edges of the object.
        // if distance <= 0, we have a collision!
        // for example, in terms of a player, we only really need one for the feet, left and right sides of the body, and the top of the head (4 separate calculations).
        
        // this kinda sucks as we would have to use more calculations for the same result, but isn't this exactly how normal collisions are calculated? 
        
        // im pretty sure there still is a simple fix to this problem with the scale..., i mean why do the merger sponge and plane, and i assume others, work perfectly? And why do some don't?


        //Check if the distance estimate indicates a collision
        if (distance <= 0.5)
        {
            // Push object back to surface of fractal
            okToJump = true;
            rb.AddForce(playerV * -1 , ForceMode.Impulse);
            //rb.AddForce(FeetN*(0.3-FeetDist)*rb.mass * 16 * invAbsFdot);
        }
        else
        {
            okToJump = false;
        }
    }
    
    public float GetDistance()
    {
        //TODO: Rotation is currently broke. FIX
        
        p = RotateZ(RotateY(RotateX(p - _fractalPosition, _fractaldegreeRotate.x), _fractaldegreeRotate.y), _fractaldegreeRotate.z);
        p = ScaleX(p, _fractalScale.x);
        p = ScaleY(p, _fractalScale.y);
        p = ScaleZ(p, _fractalScale.z);
        
        /*if (_modBool.x == 1) {
            p.x = Repeat(p.x, _modInterval.x);
        }
        if (_modBool.y == 1) {
            p.y = Repeat(p.y, _modInterval.y);
        }
        if (_modBool.z == 1) {
            p.z = Repeat(p.z, _modInterval.z);
        }*/
        switch (_fractalNumber)
        {
            case 0:
                return sdMandelbulb(p);
                break;
            case 1:
                return sdDinamMandelbulb(p);
                break;
            case 2:
                return sdJulia(p);
                break;
            case 3:
                return sdJuliabulb(p);
                break;
            case 4:
                return sierpinski(p);
                break;
            case 5:
                return mandelbox(p);
                break;
            case 6:
                return kaleidoscopic_IFS(p);
                break;
            case 7:
                return tglad_formula(p);
                break;
            case 8:
                return hartverdrahtet(p);
                break;
            case 9:
                return pseudo_kleinian(p);
                break;
            case 10:
                return pseudo_knightyan(p);
                break;
            case 11:
                return mandelbulb2(p, _power);
                break;        
            case 12:
                return MengerSponge(p);
                break;      
            case 13:
                return apo(p, .0274f, float3(1, 1, 1.3f), float3(0, 0, 0));
                break; 
            case 14:
                return sdPlane(p, float4(0,1,0,0));
                break;    
            case 15:
                return FCT_BBSK(p, Params);     
                break;    
            case 16:
                return trinoise(p);     
                break;  
            /*case 17:
                return RecursiveTetrahedron(p, 3);     
                break; */ 
            case 18:
                return TruchetTentacles(p);     
                break;
            case 19:
                return FCT_PROTEIN(p, Params);     
                break;
            case 20:
                return FCT_ORBIT(p);     
                break;
            case 21:
                return FCT_MNMT(p);     
                break;
            case 22:
                return FCT_CRAB(p);     
                break;
            case 23:
                return FCT_HUB(p, Params);     
                break;
            case 24:
                return FCT_HYPERAPO(p, Params);     
                break;
            case 25:
                return FCT_DLBT(p, Params);     
                break;
            case 26:
                return FCT_MZGN(p, Params);     
                break;
            case 27:
                return FCT_PIPES(p);     
                break;
            case 28:
                return FCT_APOP(p, Params);     
                break;
            case 29:
                return FCT_APO(p);     
                break;
            case 30:
                return FCT_HTVT(p, Params);     
                break;
            case 31:
                return FCT_KNKL(p, Params);     
                break;   
            case 32:
                return FCT_KIFS(p, Params);     
                break;
            case 33:
                return FCT_TEXT(p);     
                break;
            case 34:
                return FCT_TEST(p);     
                break;                                                                
            default:
                return apo(p, .0274f, float3(1f, 1f, 1.3f), float3(0f, 0f, 0f));
                break;         
        }
    }


    float rand(float3 r)
{
    return frac(sin(dot(r.xy, r.yz)));
}

float Mod(float a, float b)
{
    return frac(abs(a / b)) * abs(b);
}

float2 Mod(float2 a, float2 b)
{
    return frac(abs(a / b)) * abs(b);
}

float3 Mod(float3 a, float3 b)
{
    return frac(abs(a / b)) * abs(b);
}

float SmoothMin(float d1, float d2, float k)
{
    float h = exp(-k * d1) + exp(-k * d2);
    return -log(h) / k;
}

float Repeat(float pos, float span)
{
    return Mod(pos, span) - span * 0.5f;
}

float2 Repeat(float2 pos, float2 span)
{
    return Mod(pos, span) - span * 0.5f;
}

float3 Repeat(float3 pos, float3 span)
{
    return Mod(pos, span) - span * 0.5f;
}

float3 Rotate(float3 p, float angle, float3 axis)
{
    float3 a = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float r = 1.0f - c;
    float3x3 m = float3x3(
        a.x * a.x * r + c,
        a.y * a.x * r + a.z * s,
        a.z * a.x * r - a.y * s,
        a.x * a.y * r - a.z * s,
        a.y * a.y * r + c,
        a.z * a.y * r + a.x * s,
        a.x * a.z * r + a.y * s,
        a.y * a.z * r - a.x * s,
        a.z * a.z * r + c
    );
    return mul(m, p);
}

float3 TwistY(float3 p, float power)
{
    float s = sin(power * p.y);
    float c = cos(power * p.y);
    float3x3 m = float3x3(
          c, 0.0f,  -s,
        0.0f, 1.0f, 0.0f,
          s, 0.0f,   c
    );
    return mul(m, p);
}

float3 TwistX(float3 p, float power)
{
    float s = sin(power * p.y);
    float c = cos(power * p.y);
    float3x3 m = float3x3(
        1.0f, 0.0f, 0.0f,
        0.0f,   c,   s,
        0.0f,  -s,   c
    );
    return mul(m, p);
}

float3 TwistZ(float3 p, float power)
{
    float s = sin(power * p.y);
    float c = cos(power * p.y);
    float3x3 m = float3x3(
          c,   s, 0.0f,
         -s,   c, 0.0f,
        0.0f, 0.0f, 1.0f
    );
    return mul(m, p);
}

// BOOLEAN OPERATORS //

// Union
float4 opU(float4 d1, float4 d2) {
	return (d1.w < d2.w) ? d1 : d2;
}

// Subtraction
float opS(float d1, float d2) {
	return max(-d1, d2);
}

// Intersection
float opI(float d1, float d2) {
	return max(d1, d2);
}

// SMOOTH BOOLEAN OPERATORS

    float4 qsqr( in float4 a ) // square a quaterion
{
    return float4( a.x*a.x - a.y*a.y - a.z*a.z - a.w*a.w,
                 2.0f*a.x*a.y,
                 2.0f*a.x*a.z,
                 2.0f*a.x*a.w );
}

void sphereFold( float3 z,  float dz)
{
	float r2 = dot(z,z);
	if (r2 < 0.5f)
    { 
		float temp = 2.0f;
		z *= temp;
		dz*= temp;
	}
    else if (r2 < 1.0f)
    { 
		float temp = 1.0f / r2;
		z *= temp;
		dz*= temp;
	}
}

void boxFold( float3 z,  float dz)
{
	z = clamp(z, -1.0f, 1.0f) * 2.0f - z;
}


float remap(float value, float low1, float high1, float low2, float high2)
{
    return low2 + (value - low1) * (high2 - low2) / (high1 - low1);
}

float  modc(float  a, float  b) { return a - b * floor(a/b); }
float2 modc(float2 a, float2 b) { return a - b * floor(a/b); }
float3 modc(float3 a, float3 b) { return a - b * floor(a/b); }
float4 modc(float4 a, float4 b) { return a - b * floor(a/b); }


float3 RotateX(float3 p, float angle)
{
    float rad = 0.0174532925f * angle;
    float s = sin(angle);
    float c = cos(angle);
    //sincos(rad, s, c);
    return float3(p.x, c*p.y + s*p.z, -s*p.y + c*p.z);
}
float3 RotateY(float3 v, float degree)
{
	float rad = 0.0174532925f * degree;
	float cosY = cos(rad);
	float sinY = sin(rad);
	return float3(cosY * v.x - sinY * v.z, v.y, sinY * v.x + cosY * v.z);
}
float3 RotateZ(float3 p, float angle)
{
    float rad = 0.0174532925f * angle;
    float s = sin(rad);
    float c = cos(rad);
    return float3(c*p.x + s*p.y, -s*p.x + c*p.y, p.z);
}

float3 ScaleX(float3 p, float scale)
{
    return float3(p.x / scale, p.y, p.z);
}

float3 ScaleY(float3 p, float scale)
{
    return float3(p.x, p.y / scale, p.z);
}

float3 ScaleZ(float3 p, float scale)
{
    return float3(p.x, p.y, p.z / scale);
}

float2x2 rot(float a) {
	return float2x2(cos(a),sin(a),-sin(a),cos(a));	
}

float AngleBetween(float3 a, float3 b)
{
    return acos(dot(a, b));
}
float AngleBetween(float3 pos1, float3 pos2, float3 center)
{
    return AngleBetween(
        normalize(pos1 - center),
        normalize(pos2 - center));
}

float fract(float x) {
 return x - floor(x);
}

float3 fract(float3 x) {
    return x - floor(x);
}

float mix(float x, float y, float z) {
 return x*(1f - z) + y*z;
}

void domainRep3(float3 p, float3 size) {
	p = modc(p + size,2f * size) - size;
}

float3 domainRep3Idx( float3 p, float3 size) {
	float3 idx = floor((p + size)/size/2f);
	p = modc(p + size, 2f * size) - size;
	return idx;
}

float smin( float a, float b, float k )
{
    float h = clamp( 0.5f+0.5f*(b-a)/k, 0.0f, 1.0f );
    return mix( b, a, h ) - k*h*(1.0f-h);
}

float smax( in float a, in float b, in float s ){
    float h = clamp( 0.5f + 0.5f*(a-b)/s, 0.0f, 1.0f );
    return mix(b, a, h) + h*(1.0f-h)*s;
}

float hash(float h) {
	return fract(sin(h) * 43758.5453123f);
}

float noise3d(float3 x) {
	float3 p = floor(x);
	float3 f = fract(x);
	f = f * f * (3 - 2 * f);

	float n = p.x + p.y * 157 + 113 * p.z;
	return mix(
			mix(mix(hash(n + 0.0f), hash(n + 1.0f), f.x),
					mix(hash(n + 157), hash(n + 158.0f), f.x), f.y),
			mix(mix(hash(n + 113.0f), hash(n + 114.0f), f.x),
					mix(hash(n + 270.0f), hash(n + 271.0f), f.x), f.y), f.z);
}

    

float sdSphere(float3 p, float s) {
	return length(p) - s;
}

// Box
float Box(float3 pos, float3 size)
{
    float3 d = abs(pos) - size;
    return length(max(abs(pos) - size, 0.0f))
        + min(max(d.x, max(d.y, d.z)), 0.0f);
}

// Rounded Box
float sdRoundBox(in float3 p, in float3 b, in float r) {
	float3 q = abs(p) - b;
	return min(max(q.x, max(q.y, q.z)), 0.0f) + length(max(q, 0.0f)) - r;
}

float Cylinder(float3 pos, float2 r)
{
    float2 d = abs(float2(length(pos.xy), pos.z)) - r;
    return min(max(d.x, d.y), 0.0f) + length(max(d, 0.0f)) - 0.1f;
}

float HexagonalPrismX(float3 pos, float2 h)
{
    float3 p = abs(pos);
    return max(
        p.x - h.y, 
        max(
            (p.z * 0.866025f + p.y * 0.5f),
            p.y
        ) - h.x
    );
}

float HexagonalPrismY(float3 pos, float2 h)
{
    float3 p = abs(pos);
    return max(
        p.y - h.y, 
        max(
            (p.z * 0.866025f + p.x * 0.5f),
            p.x
        ) - h.x
    );
}

float HexagonalPrismZ(float3 pos, float2 h)
{
    float3 p = abs(pos);
    return max(
        p.z - h.y, 
        max(
            (p.x * 0.866025f + p.y * 0.5f),
            p.y
        ) - h.x
    );
}

// Torus
float sdTorus(float3 p, float2 t) {
	float2 q = float2(length(p.xz)-t.x, p.y);
	return length(q) - t.y;
}

// (Infinite) Plane
// n.xyz: normal of the plane (normalized)
// n.w: offset
float sdPlane(float3 p, float4 n) {
	// n must be normalized
	return dot(p, n.xyz) + n.w;
}

// Subtraction
float opSS(float d1, float d2, float k) {
	float h = clamp(0.5f - 0.5f * (d2 + d1) / k, 0.0f, 1.0f);
	return lerp(d2, -d1, h) + k * h * (1.0f - h);
}

// Intersection
float opIS(float d1, float d2, float k) {
	float h = clamp(0.5f - 0.5f * (d2 - d1) / k, 0.0f, 1.0f);
	return lerp(d2, d1, h) + k * h * (1.0f - h);
}

//static mandelbulb
float sdMandelbulb(float3 p)
{
	float3 w = p;
    float m = dot(w, w);

	float dz = 1.0f;
        
	for(int i = 0; i < 2; i++)
    {
        dz = 8 * pow(sqrt(m), 7.0f)*dz + 1.0f;
        float r = length(w);
        float b = 8 * acos(w.y / r);
        float a = 8 * atan2(w.x, w.z);
        w = p + pow(r, 8) * float3(sin(b) * sin(a), cos(b), sin(b) * cos(a));

        m = dot(w, w);
		if(m > 256.0f)
            break;
    }
    return 0.25f*log(m)*sqrt(m)/dz;
}

float sdDinamMandelbulb(float3 pos)
{
    float sinTime = sin(Time.timeSinceLevelLoad / 1);
    float power = remap(sinTime, -1, 1, 4, 9);
    float3 z = pos;
    float r = 0;
    float dr = 1;
    for(int i = 0; i < 5; i++) 
    {
        r = length(z);
        if(r > 100) break;
        
        float theta = acos(z.z / r);
        float phi = atan2(z.y, z.x);
        
        dr = power * pow(r, power-1)*dr+1;
        
        r = pow(r, power);
        theta *= power;
        phi *= power;
        
        z = r * float3(sin(theta) * cos(phi), 
                sin(theta) * sin(phi), 
                cos(theta));

        z += pos;
    }
    return 0.5f * log(r) * r / dr;

}

float mandelbulb2(float3 pos, float power) {
    float3 z = pos;
	float dr = 1.0f;
	float r = 0.0f;
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
		dr =  pow( r, power-1.0f)*power*dr + 1.0f;

		// scale and rotate the point
		float zr = pow( r,power);
		theta = theta*power;
		phi = phi*power;
		
		// convert back to cartesian coordinates
		z = zr*float3( cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta) );
		z+=pos;
	}
    float dst = 0.5f*log(r)*r/dr;
	return dst;
}

// MengerSponge distance estimation:
// http://www.iquilezles.org/www/articles/menger/menger.htm
float MengerSponge( in float3 p )
{
   float d = Box(p, 50.0f); // size of box

   float s = 0.02f; // size of squares in box
   int iterations = 0;
   
   for( int m=0; m<5; m++ )
   {
      //p = p + float3(power-1.0,power-1.0,power-1.0f);
      iterations = m;
      float3 a = (p*s - 2.0f * floor(p*s/2.0f))-1.0f;
      s *= 3.0f;
      float3 r = abs(1.0f - 3.0f*abs(a));

      float da = max(r.x,r.y);
      float db = max(r.y,r.z);
      float dc = max(r.z,r.x);
      float c = (min(da,min(db,dc))-1.0f)/s;

      d = max(d,c);
   }

   return d;
}

float sdJulia(float3 pos)
{
    float4 c = 0.45f* cos( float4(0.5f,3.9f,1.4f,1.1f) + Time.timeSinceLevelLoad * float4(1.2f,1.7f,1.3f,2.5f) ) - float4(0.3f,0.0f,0.0f,0.0f);
	float4 z = float4(pos, 0);
    float md2 = 1;
    float mz2 = dot(z, z);

	//[loop]
    for(int i = 0; i < 11; i++)
    {
        md2 *= 4.0f * mz2; // dz -> 2·z·dz, meaning |dz| -> 2·|z|·|dz| (can take the 4 out of the loop and do an exp2() afterwards)
        z = qsqr(z) + c; // z  -> z^2 + c

        mz2 = dot(z,z);

        if(mz2 > 4.0f) break;
    }
    
    return 0.25f * sqrt(mz2/md2) * log(mz2);
}

float sdJuliabulb(float3 pos)
{
    float4 c = 0.45f* cos( float4(0.5f,3.9f,1.4f,1.1f) + Time.timeSinceLevelLoad * float4(1.2f,1.7f,1.3f,2.5f) ) - float4(0.3f,0.0f,0.0f,0.0f);
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
                    r*cos(o)) + c.xyz;
        
        if (dot(orbit, orbit) > 4.0f) break;
    }
    float z = length(orbit);
    return 0.5f*z*log(z)/dz;
}

float dm;
float sierpinski(float3 p)
{
     float3 va = float3(  0.0f,  0.575735f,  0.0f );
     float3 vb = float3(  0.0f, -1.0f,  1.15470f );
     float3 vc = float3(  1.0f, -1.0f, -0.57735f );
     float3 vd = float3( -1.0f, -1.0f, -0.57735f );

    float a = 0;
    float s = 1;
    float r = 1;

    
    float3 v;
    for(int i = 0; i < 15; i++)
	{
	    float d, t;
		d = dot(p - va, p - va);

        v = va; 
        dm = d; 
        t = 0;
        
        d = dot(p - vb, p - vb); 
        if(d < dm) 
        { 
            v = vb; 
            dm=d; 
            t = 1.0f;
        }
        
        d = dot(p-vc, p-vc); 

        if(d < dm) { v = vc; dm = d; t = 2.0f;}
        d = dot(p-vd,p-vd); 
        if(d < dm) { v = vd; dm = d; t = 3.0f;}

		p = v + 2*(p - v); 
        r*= 2;
		a = t + 4*a; 
        s*= 4;
	}
	
	//return float2((sqrt(dm)-1.0f)/r, a/s);
    return sqrt(dm)-1.0f;
}

float mandelbox(float3 position) {
    float SCALE = 2.75f;
    float fixedRadius = 1.0f;
    float FR2 = fixedRadius * fixedRadius;
    float minRadius = 0.5f;
    float MR2 = minRadius * minRadius;
    float4 scalefloat = float4(SCALE, SCALE, SCALE, abs(SCALE)) / MR2;
    float C1 = abs(SCALE-1.0f);
    float C2 = pow(abs(SCALE), -4);
    float4 p = float4(position.xyz, 1.0f); 
    float4 p0 = float4(position.xyz, 1.0f);  // p.w is knighty's DEfactor
    for (int i=0; i<5; i++) {
        p.xyz = clamp(p.xyz *0.5f+0.5f, 0.0f, 1.0f) *4.0f-2.0f - p.xyz; // box fold: min3, max3, mad3
        float r2 = dot(p.xyz, p.xyz);  // dp3
        p.xyzw *= clamp(max(MR2/r2, MR2), 0.0f, 1.0f);  // sphere fold: div1, max1.sat, mul4
        p.xyzw = p*scalefloat + p0;  // mad4
    }
  return (length(p.xyz) - C1) / p.w - C2;

}

float mandelbox2(float3 p)
{
    float scale = 2;
	float3 offset = p;
	float dr = 1.0f;
	for (int n = 0; n < 10; n++)
    {
		boxFold(p, dr);
		sphereFold(p, dr);
        p = scale * p + offset;
        dr = dr * abs(scale) + 1.0f;
	}
	float r = length(p);
	return r / abs(dr);
}

float kaleidoscopic_IFS(float3 z)
{
    float FRACT_ITER      = 20;
    float FRACT_SCALE   = 1.8f;
    float FRACT_OFFSET  = 1.0f;

    float c = 2.0f;
    z.y = modc(z.y, c)-c/2.0f;
    z = RotateZ(z, (PI/2.0f)/0.0174532925f);
    float r;
    int n1 = 0;
    for (int n = 0; n < FRACT_ITER; n++) {
        float rotate = PI*0.5f;
        z = RotateX(z, rotate/0.0174532925f);
        z = RotateY(z, rotate/0.0174532925f);
        z = RotateZ(z, rotate/0.0174532925f);

        z.xy = abs(z.xy);
        if (z.x+z.y<0.0f) z.xy = -z.yx; // fold 1
        if (z.x+z.z<0.0f) z.xz = -z.zx; // fold 2
        if (z.y+z.z<0.0f) z.zy = -z.yz; // fold 3
        z = z*FRACT_SCALE - FRACT_OFFSET*(FRACT_SCALE-1.0f);
    }
    return (length(z) ) * pow(FRACT_SCALE, - FRACT_ITER);
}

float tglad_formula(float3 z0)
{
    z0 = modc(z0, 2.0f);

    float mr=0.25f, mxr=1.0f;
    float4 scale=float4(-3.12f,-3.12f,-3.12f,3.12f), p0=float4(0.0f,1.59f,-1.0f,0.0f);
    float4 z = float4(z0,1.0f);
    for (int n = 0; n < 3; n++) {
        z.xyz=clamp(z.xyz, -0.94f, 0.94f)*2.0f-z.xyz;
        z*=scale/clamp(dot(z.xyz,z.xyz),mr,mxr);
        z+=p0;
    }
    float dS=(length(max(abs(z.xyz)-float3(1.2f,49.0f,1.4f),0.0f))-0.06f)/z.w;
    return dS;
}


// distance function from Hartverdrahtet
// ( http://www.pouet.net/prod.php?which=59086 )
float hartverdrahtet(float3 f)
{
    float3 cs=float3(.808f,.808f,1.167f);
    float fs=1;
    float3 fc=0;
    float fu=10;
    float fd=.763f;
    
    // scene selection
    {
        float time = Time.timeSinceLevelLoad;
        float i = modc(time/2.0f, 9.0f);
        if(i==0) cs.y=.58f;
        if(i==1) cs.xy=.5f;
        if(i==2) cs.xy=.5f;
        if (i == 3)
        {
            fu = 1.01f;
            cs.x=.9f;
        }
        if (i == 4)
        {
            fu = 1.01f;
            cs.x=.9f;
        }
        if(i==6) cs=float3(.5f,.5f,1.04f);
        if(i==5) fu=.9f;
        if (i == 7)
        {
            fd = .7f;
            fs = 1.34f;
            cs.xy=.5f;
        }
        if(i==8) fc.z=-.38f;
    }
    
    //cs += sin(time)*0.2;

    float v=1;
    for(int i=0; i<12; i++){
        f=2*clamp(f,-cs,cs)-f;
        float c=max(fs/dot(f,f),1);
        f*=c;
        v*=c;
        f+=fc;
    }
    float z=length(f.xy)-fu;
    return fd*max(z,abs(length(f.xy)*f.z)/sqrt(dot(f,f)))/abs(v);
}

float pseudo_kleinian(float3 p)
{
    float3 CSize = float3(0.92436f,0.90756f,0.92436f);
    float Size = 1.0f;
    float3 C = float3(0.0f,0.0f,0.0f);
    float DEfactor=1;
    float3 Offset = float3(0.0f,0.0f,0.0f);
    float3 ap=p+1;
    for(int i=0;i<10 ;i++){
        ap=p;
        p=2*clamp(p, -CSize, CSize)-p;
        float r2 = dot(p,p);
        float k = max(Size/r2,1);
        p *= k;
        DEfactor *= k + 0.05f;
        p += C;
    }
    float r = abs(0.5f*abs(p.z-Offset.z)/DEfactor);
    return r;
}

float pseudo_knightyan(float3 p)
{
    float3 CSize = float3(0.63248f,0.78632f,0.875f);
    float DEfactor=1;
    for(int i=0;i<6;i++){
        p = 2*clamp(p, -CSize, CSize)-p;
        float k = max(0.70968f/dot(p,p),1);
        p *= k;
        DEfactor *= k + 0.05f;
    }
    float rxy=length(p.xy);
    return max(rxy-0.92784f, abs(rxy*p.z) / length(p))/DEfactor;
}

// trinoise
float tri(float x)
{
    return abs(frac(x) - .5f);
}

float3 tri3(float3 p)
{
    return float3(
        tri(p.z + tri(p.y * 1)),
        tri(p.z + tri(p.x * 1)),
        tri(p.y + tri(p.x * 1)) );
}

float trinoise(float3 p, float spd, float time)
{
    float z = 1.4f;
    float rz = 0;
    float3  bp = p;
    for (float i = 0; i <= 3; i++) {
        float3 dg = tri3(bp * 2);
        p += (dg + time * .1f * spd);
        bp *= 1.8f;
        z *= 1.5f;
        p *= 1.2f;
        float t = tri(p.z + tri(p.x + tri(p.y)));
        rz += t / z;
        bp += 0.14f;
    }
    return rz;
}

float trinoise(float3 p)
{
    return trinoise(p, 1.0f, 1.0f);
}
/*
float RecursiveTetrahedron(float3 p, int loop)
{

    p = Repeat(p / 2, 3.0f);

    const float3 a1 = float3( 1.0,  1.0,  1.0f);
    const float3 a2 = float3(-1.0, -1.0,  1.0f);
    const float3 a3 = float3( 1.0, -1.0, -1.0f);
    const float3 a4 = float3(-1.0,  1.0, -1.0f);

    const float scale = 2.0f;
    float d;
    for (int n = 0; n < loop; ++n) {
        float3 c = a1; 
        float minDist = length(p - a1);
        d = length(p - a2); if (d < minDist) { c = a2; minDist = d; }
        d = length(p - a3); if (d < minDist) { c = a3; minDist = d; }
        d = length(p - a4); if (d < minDist) { c = a4; minDist = d; }
        p = scale * p - c * (scale - 1.0f);
    }
 
    return length(p) * pow(scale, float(-n));
    
    return 0;
}
*/
// original code: https://www.shadertoy.com/view/ldfGWn
float truchetarc(float3 pos)
{
    float p = 4.0f + 2.0f * sin(Time.timeSinceLevelLoad);
    float r = length(pos.xy);
    float t = 0.12f + 0.02f * sin(Time.timeSinceLevelLoad);
    return pow(
        pow(abs(r - 0.5f), p) + pow(abs(pos.z - 0.5f), p), 
        rcp(p)
    ) - t;
}
float truchetcell(float3 pos)
{
    return min(min(
        truchetarc(pos),
        truchetarc(float3(pos.z, 1.0f - pos.x, pos.y))), 
        truchetarc(float3(1.0f - pos.y, 1.0f - pos.z, pos.x)));
}
float TruchetTentacles(float3 pos)
{
    float3 c = frac(pos);
    float  r = rand(floor(pos));

    if      (r < 0.125) return truchetcell(float3(c.x, c.y, c.z));
    else if (r < 0.250) return truchetcell(float3(c.x, 1.0f - c.y, c.z));
    else if (r < 0.375) return truchetcell(float3(1.0f - c.x, c.y, c.z));
    else if (r < 0.500) return truchetcell(float3(1.0f - c.x, 1.0f - c.y, c.z));
    else if (r < 0.625) return truchetcell(float3(c.y, c.x, c.z));
    else if (r < 0.750) return truchetcell(float3(c.y, 1.0f - c.x, c.z));
    else if (r < 0.875) return truchetcell(float3(1.0f - c.y, c.x, c.z));
    else                return truchetcell(float3(1.0f - c.y, 1.0f - c.x, c.z));
}

//uniform float3 cFcParams = float3(0,0,0);
//uniform float3x3 cFcRot;
float dist = 10000f;
float tex,tex2,tex3 = 0;


float apo(float3 pos, float seed, float3 CSize, float3 C) 
{
  float dist;
  //float3 CSize = float3(1., 1., 1.3);
  float3 p = pos.xzy;
  float scale = 1.0f;
 // p *= cFcRot;
  float r2 = 0;
  float k = 0;
  //float uggg = 0.;
  for( int i=0; i < 12;i++ )
  {
      p = 2.0f*clamp(p, -CSize, CSize) - p;
      r2 = dot(p,p);
      //r2 = dot(p,p+sin(p.z*.3)); //Alternate fractal
      k = max((2.0f)/(r2), seed); //.378888 //.13345 max((2.6)/(r2), .03211); //max((1.8)/(r2), .0018);
      p     *= k;
      scale *= k;
      //uggg += r2;
      p+=C;
       //p.xyz = float3(-1.0f*p.z,1.0f*p.x,1.0f*p.y);
  }
  float l = length(p.xy);
  float rxy = l - 4.0f;
  float n = 1.0f * p.z;
  rxy = max(rxy, -(n) / 4);
  dist = (rxy) / abs(scale);
  return dist;
}

float FCT_PROTEIN(float3 pos, float3 Params) {
// Notable Params: (-12.82 -0.63 -16.18)
    float3 cFcParams = Params;
    float3 p = pos * 0.002f;
    float4 q = float4(p - 1, 1);
    for(int i = 0; i < 7; i++) {
      q.xyz = abs(q.xyz + 1.3f) - 1.0f;
      q /= clamp(dot(q.xyz, q.xyz), 0.0f, 0.8f);
      //q.xyz *= cFcRot;
      q *= 1.567f+cFcParams.x;// + p.y*0.8;
      q = q.zxyw;
      q+= float4(0.2119f,0.8132f,0.077213f,0);
    }
    for(int i = 0; i < 4; i++) {
      q.xyz = abs(q.xyz + 1.3f) - 1.0f;
      q /= clamp(dot(q.xyz, q.xyz), 0.0f, 0.8f);
      q *= 1.867f;// + p.y*0.8;
    }
    return (length(q.xyz) - max(-240-p.y*700,2.5f))/q.w * 300;
}

float FCT_ORBIT(float3 pos) {
    float3 p = pos * 0.001f;
    float4 q = float4(p - 1.0f, 1);
    for(int i = 0; i < 11; i++) {
      //q.xyz = mod(q.xyz,2.0f)-0.5*q.xyz;
       q.xyz = abs(q.xyz + 1.0f) - 1.0f;
      q.xyz = 2.0f*clamp(q.xyz, -73.0174f, 73.0174f) - q.xyz;
      q /= min(dot(q.xyz, q.xyz), 0.9f);
      q *= 1.817f;// + p.y*0.8;
      //q += float4(0.2,.02,-.2,0.2);
      //q.xyz *= rot;
    }
    return dist = (length(q.xz) - 2.1f)/q.w * 800;
}

float FCT_MNMT(float3 pos) {
    float3 CSize = float3(1, 1, 1.3f); // <-- CSize Constant
float3 p = pos.yzx; // <-- 3D Position
float scale = 1.0f;// <-- Scale

for (int i = 0; i < 8; i++) { // <-- Primary iteration
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    float r2 = dot(p, p);
    float k = max((2) / (r2), .1274f);
    p *= k;
    //p *= cFcRot;                                    // \ Lines present only
    p.xyz = float3(1.0f * p.z, 1.0f * p.x, -1.0f * p.y); // / in this loop
    scale *= k;
}

CSize = float3(1.2f, 0.4f, 1.4f); // <-- CSize Constant
tex = p.y; //<-- Texture Primary

for (int i = 0; i < 4; i++) { // <-- Secondary Iteration
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    float r2 = dot(p, p);
    float k = max((1.6f) / (r2), 0.0274f);
    p *= k;
    scale *= k;
}

float l = length(p.xyz);
float rxy = l - 1.4f;
float n = 1.0f * p.z;
rxy = max(rxy, -(n) / 4);
dist = (rxy) / abs(scale);

return dist * 1.5f;
}

float FCT_CRAB(float3 pos) {
float3 p = pos * 0.002f; // <-- 3D Position
float4 q = float4(p - 1.0f, 1); // <-- 4D Map onto Position

for (int i = 0; i < 7; i++) { // <-- Primary iteration
    q.xyz = abs(q.xyz + 1.3f) - 1.0f;
    q /= clamp(dot(q.xyz, q.xyz), 0.0f, 0.8f);
    q *= 1.867f; // + p.y*0.8;
    q = q.zxyw;                              // \ Lines present only
    q += float4(0.2119f, 0.8132f, 0.077213f, 0); // / in this loop
}

for (int i = 0; i < 4; i++) { // <-- Secondary iteration
    q.xyz = abs(q.xyz + 1.3f) - 1.0f;
    q /= clamp(dot(q.xyz, q.xyz), 0.0f, 0.8f);
    q *= 1.867f; // + p.y*0.8;
}

return (length(q.yz) - 1.2f) / q.w * 300; // <-- Distance
}

float FCT_BBSK(float3 pos, float3 Params) {
// Notable Params: (2.18 -0.18 0)
    //float3 cFcParams = float3(2.18, -0.18, 0);
    float3 cFcParams = Params;
    float3 CSize = float3(1.4f,0.87f, 1.1f);
    float3 p = pos.xzy * 2.0f;
    float scale = 1.0f;
    
    for( int i=0; i < 4;i++ )
    {
        p = 2.0f*clamp(p, -CSize, CSize) - p;
        //float r2 = dot(p,p);
        float r2 = dot(p,p+sin(p.z*.5f)); //Alternate fractal
        float k = max((2)/(r2), .17f);
        p *= k;
        //p *=rot;
        //p= p.yzx;
        p+=float3(0.2f,0.2f,-0.5f);
        scale *= k;
    }

    p = 2.0f*clamp(p, -CSize * 4, CSize * 4) - p;
   
    for(int i=0; i < 8; i++ )
    {
        p = 2.0f*clamp(p, -CSize, CSize) - p;
        float r2 = dot(p,p);
        //float r2 = dot(p,p+sin(p.z*.3)); //Alternate fractal
        float k = max((cFcParams.x)/(r2),  0.027f);
        p     *= k;
        scale *= k;
        p.y += cFcParams.y;
    }
    
    float l = length(p.xy);
    //l = mix(l,l2,0.5);
    float rxy = l - 4.0f;
    float n = 1.0f * p.z;
    rxy = max(rxy, -(n) / 4);
    float dist = (rxy) / abs(scale);
    dist *=.75f;

    return dist;

}

float FCT_HUB(float3 pos, float3 Params) {
// Notable Params: (0.514 -2.28 0.83), (0.644 -2.28 0.83)
float3 cFcParams = Params;

float3 p = pos.yzx * 0.05f; // <-- 3D Position
float4 q = float4(p - 1.0f, 1); // <-- 4D Map onto Position

for (int i = 0; i < 8; i++) { // <-- Primary iteration (All lines shared)
    q.xyz = abs(q.xyz + cFcParams.z) - 1.0f;
    q /= clamp(dot(q.xyz, q.xyz), 0.12f, 1.0f);
    q *= 1.0f + cFcParams.x;
    //q.xyz *= cFcRot;
}

tex = q.x; // <-- Texture Primary

for (int i = 0; i < 2; i++) { // <-- Secondary iteration
    q.xyz = abs(q.xyz + cFcParams.z) - 1.0f;
    q /= clamp(dot(q.xyz, q.xyz), 0.12f, 1.0f);
    q *= 1.0f + cFcParams.x;
    //q.xyz *= cFcRot;
}

return (length(q.xyz) - 1 + cFcParams.y) / q.w * 19; // <-- Distance
}

float FCT_HYPERAPO(float3 pos, float3 Params) {
// Notable Params: (0.644 -2.28 0.83), (0.067 1.05 -5.58), (0.76 0.91 4.22)
float3 cFcParams = Params;
float3 CSize = float3(1, 1, 1.1f); // <-- CSize Constant
float3 p = pos.yzx; // <-- 3D Position
float scale = 1.0f;// <-- Scale

for (int i = 0; i < 4; i++) { // <-- Primary iteration
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    float r2 = dot(p, p);
    float k = max((2) / (r2), 0.067f);
    p *= k;
    p.xyz = float3(1.0f * p.z, 1.0f * p.x, -1.0f * p.y); // Line present only in this loop
    scale *= k;
}

p = p.zxy;
//p *= cFcRot;

CSize = float3(1.2f, cFcParams.x, 1.2f); // <-- CSize Constant
tex = p.y; // <-- Texture Primary

for (int i = 0; i < 7; i++) { // <-- Secondary iteration
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    float r2 = dot(p, p);
    float k = max((1.6f) / (r2), cFcParams.y);
    p *= k;
    scale *= k;
}

float l = length(p.xyz);
float rxy = l - 1.4f;
float n = 1.0f * p.z;
rxy = max(rxy, -(n) / 4);

return (rxy) / abs(scale); // <-- Distance
}

float FCT_DLBT(float3 pos, float3 Params) {
// Notable Params: (-24.63 16.16 101.07)
float3 cFcParams = Params;
float3 npos = pos;
float noise = noise3d(npos * 0.07f) * 10;
float3 apopos = pos.xzy;
apopos.z += cFcParams.x;
domainRep3(apopos, float3(250, 250, 0));
float r2 = 0;
float k = 0;

float3 CSize = float3(1, 1, 1.3f); // <-- CSize Constant
float3 p = apopos; // <-- 3D Position
float scale = 1.0f;// <-- Scale

for (int i = 0; i < 5; i++) { // <-- Primary iteration (All lines shared)
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    r2 = dot(p, p);
    k = max((2.0f) / (r2), .0274f);
    p *= k;
    scale *= k;

}

tex = scale; // <-- Texture Primary

for (int i = 0; i < 7; i++) { // <-- Secondary iteration
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    r2 = dot(p, p);
    k = max((2.0f) / (r2), .0274f);
    p *= k;
    scale *= k;
}

float l = length(p.xy);
float rxy = l - 4.0f;
float n = 1.0f * p.z;
rxy = max(rxy, -(n) / 4);
float apodist = (rxy) / abs(scale);
apodist = 2 - apodist * 2;
float dist = (float)(noise - apodist + 0.01 * (pos.y - 67));
dist *= 0.7f;
dist = smin(dist, pos.y + cFcParams.y, 15f);

tex2 = p.z / scale; // <-- Texture Secondary
return smin(dist, pos.y - cFcParams.z, -15); // <-- Distance
}

float FCT_MZGN(float3 pos, float3 Params) {
// Notable Params: A lot
float3 cFcParams = Params;
float3 p = pos * 0.01f; // <-- 3D Position
float4 q = float4(p, 1); // <-- 4D Map onto Position

float4 qd;

for (int i = 0; i < 12; i++) { // <-- Primary iteration
    q.xyz = abs(q.xyz + 1.0f + cFcParams.y) - 1.0f;
    qd = q;
    qd.w = qd.w * 0.002f;
    q /= clamp(dot(qd, qd), 0.0f, 0.8f);
   // q.xyz *= cFcRot;
    q *= 1.567f + cFcParams.x;
}

tex = q.w; // <-- Texture Primary
tex2 = q.x; // <-- Texture Secondary
return (length(q.xyz) - 3.5f) / q.w * 100; // <-- Distance
}

float FCT_PIPES(float3 pos) {
float3 p = pos * 0.002f; // <-- 3D Position
float4 q = float4(p - 1.0f, 1); // <-- 4D Map onto Position

for (int i = 0; i < 11; i++) { // <-- Primary iteration
    q.xyz = abs(q.xyz + 1.0f) - 1.0f;
    q /= clamp(dot(q.xyz, q.xyz), 0.12f, 1.0f);
    q *= 1.837f;
}

tex = q.y; // <-- Texture Primary
tex2 = q.w; // <-- Texture Secondary
return (length(q.xz) - 1.2f) / q.w * 500f; // <-- Distance
}

float FCT_APOP(float3 pos, float3 Params) {
// Notable Params: (1.11 4.48 0.07)
float3 cFcParams = Params;
float3 p = pos * 0.01f; // <-- 3D Position
float scale = 1.0f;// <-- Scale

for (int i = 0; i < 10; i++) { // <-- Primary iteration
    p = -1.0f + 2.0f * fract(0.5f * p + 0.5f);
    float r2 = dot(p, p);
    float k = cFcParams.x / r2;
    p *= k;
    scale *= k;
    //p *= cFcRot;
}

dist = length(p) / scale;

return dist * 25; // <-- Distance
}

float FCT_APO(float3 pos) {
float3 CSize = float3(1.7f, 1.7f, 1.3f); // <-- CSize Constant
float3 p = pos.xzy; // <-- 3D Position
float scale = 1; // <-- Scale

for (int i = 0; i < 12; i++) { // <-- Primary iteration
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    float r2 = dot(p, p);
    float k = max((2f) / (r2), .027f);
    p *= k;
    scale *= k;
}

float l = length(p.xyz);
float rxy = l - 4.0f;
float n = l * p.z;
rxy = max(rxy, -(n) / 4f);
dist = (rxy) / abs(scale);

return max(dist, pos.x); // <-- Distance
}

float FCT_HTVT(float3 pos, float3 Params) {
// Notable Params: (0.067 0.75 -0.02), (0.09 2.72 -0.51)
float3 cFcParams = Params;
float3 CSize = float3(1f, 1f, 1.3f); // <-- CSize Constant
float3 p = pos.yzx; // <-- 3D Position
float scale = 1.0f;// <-- Scale

for (int i = 0; i < 4; i++) { // <-- Primary iteration
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    float r2 = dot(p, p);
    float k = max((2f) / (r2), cFcParams.x);
    p *= k;
    p.xyz = float3(1.0f * p.z, 1.0f * p.x, -1.0f * p.y); // Line present only in this loop
    scale *= k;
}

p = p.zxy;
//p *= cFcRot;

CSize = float3(1.2f, 0.4f, 1.4f); // <-- CSize Constant
tex2 = p.x; // <-- Texture Secondary

for (int i = 0; i < 8; i++) { // <-- Secondary iteration
    p = 2.0f * clamp(p, -CSize, CSize) - p;
    float r2 = dot(p, p);
    float k = max((1.6f) / (r2), cFcParams.y);
    p *= k;
    scale *= k;
}

float l = length(p.xyz);
float rxy = l - 1.4f;
float n = 1.0f * p.z;
rxy = max(rxy, -(n) / 4f);
dist = (rxy) / abs(scale) - 0.0005f;

tex = min(scale, 150f); // <-- Texture Primary
return dist * 1.4f; // <-- Distance
}

// Knighty's Pseudo Kleinian Fractal
float FCT_KNKL(float3 pos, float3 Params) {
// Notable Params: (-0.05 -0.37 0.07), (-0.84 0.97 0.13)
float3 cFcParams = Params;
float3 CSize = float3(0.97478f, 1.4202f, 0.97478f); // <-- CSize Constant
float3 p = pos.xzy * 0.01f; // <-- 3D Position
float3 C = cFcParams;

float DEfactor = 1;
float3 ap = p + float3(1,1,1);

for (int i = 0; i < 11; i++) { // <-- Primary iteration
    if (ap.x != p.x && ap.y != p.y && ap.z != p.z) {continue;}
    ap = p;
    p = 2 * clamp(p, -CSize, CSize) - p;
    float r2 = dot(p, p);
    float k = max(1.0f / r2, 1);
    p *= k;
    DEfactor *= k;
    p += C;
}

tex = p.y; // <-- Texture Primary
tex2 = DEfactor; // <-- Texture Secondary
return (0.5f * (p.z) / DEfactor) * 70; // <-- Distance
}

float FCT_KIFS(float3 pos, float3 Params) {
// No notable Params, just experiment!
float3 cFcParams = Params;
float3 p = pos * 0.01f; // <-- 3D Position

float l = 0.0f;
float3 Fold = float3(3,3,3);

for (int i = 0; i < 12; i++) { // <-- Primary iteration
    p.xyz = abs(p.xyz + Fold.xyz) - Fold.xyz;
    p = p * cFcParams.x;
    //p *= cFcRot;
}

l = length(p);

return (l * pow(cFcParams.x, -12) - 0.001f) * 100; // <-- Distance
}


float FCT_TEXT(float3 pos) {
float3 p = pos; // <-- 3D Position

float3 cell = domainRep3Idx(p, float3(12, 8, 12));
float3 h = float3(hash(cell.y), hash(cell.z + cell.x), hash(cell.x));
float rot = h.x + h.y + h.z * 6.3f;
float2 sc = float2(sin(rot), cos(rot));
p.xz = float2(sc.y * p.x - sc.x * p.z, sc.y * p.z + sc.x * p.x);
float3 cube = abs(p) - float3(3,3,3);
cube *= 0.9f;
dist = length(max(cube, 0.0f)) - 1.5f;
p += float3(5,5,5);
dist = max(dist, pos.y - 54);
if (dist < 2) {
    //float4 vol = texture3D(sVolumeMap, clamp(p.zxy * 0.1, 0., 1.)) * 0.9;
    float d = 2; //* (vol.r + vol.g + vol.b);
    dist = smax(dist - 1.5f, d * 1.0f - 0.1f, 2f);
} else {
    dist += 1f;
}

return min(dist, pos.y + 4f); // <-- Distance
}


float FCT_TEST(float3 pos) {
float3 pos1 = pos; // <-- 3D Position
float3 pos2 = pos; // <-- 3D Position

domainRep3(pos1, float3(45f, 45, 45));
domainRep3(pos2, float3(5, 5, 5));
dist = length(pos1) - 35;

tex = dist; // <-- Texture Primary
tex2 = length(pos2); // <-- Texture Secondary

dist = smin(dist, tex2 - 5, -4);

return min(dist, pos.y); // <-- Distance
}

    
    /*float sdPlane(Vector3 p, Vector4 n) {
        // n must be normalized
        // n is the normal, meaning n = float(0,1,0,0) is pointing up on the y
        return Vector3.Dot(p, new Vector3(n.x, n.y, n.z)) + n.w;
    }
    
    // TODO: fix this mess
    /*
    Vector3 getNormal(Vector3 p)
    {
        
        Vector2 offset = new Vector2(0.001f, 0.0f);
        Vector3 n = new Vector3(
            FCT_BBSK(new Vector3(p.x + offset.x, p.y + offset.y, p.z + offset.y)) - FCT_BBSK(p - offset.xyy).w,
            FCT_BBSK(p + offset.yxy).w - FCT_BBSK(p - offset.yxy).w,
            FCT_BBSK(p + offset.yyx).w - FCT_BBSK(p - offset.yyx).w);
        return Vector3.Normalize(n);
    }
    #1#
    
    float sdBox(Vector3 p, Vector3 b)
    {
        p.x = Mathf.Abs(p.x);
        p.y = Mathf.Abs(p.y);
        p.z = Mathf.Abs(p.z);
        Vector3 d = p - b;

        Vector3 e;
        e.x = Mathf.Max(d.x, 0);
        e.y = Mathf.Max(d.y, 0);
        e.z = Mathf.Max(d.z, 0);

        float length = (float) Mathf.Sqrt((e.x * e.x) + (e.y * e.y) + (e.z * e.z));
        
        return (float)(Mathf.Min(Mathf.Max(d.x, Mathf.Max(d.y, d.z)), 0) +
                       length);
    }
    
    float MengerSponge( Vector3 p )
    {
        float d = sdBox(p, new Vector3(50,50,50)); // size of box

        float s = 0.02f; // size of squares in box

        for( int m=0; m<5; m++ )
        {
            Vector3 a;
            a.x = (float)((p.x * s - 2.0f * Mathf.Floor(p.x * s / 2.0f)) - 1.0f);
            a.y = (float)((p.y * s - 2.0f * Mathf.Floor(p.y * s / 2.0f)) - 1.0f);
            a.z = (float)((p.z * s - 2.0f * Mathf.Floor(p.z * s / 2.0f)) - 1.0f);
            
            s *= 3.0f;

            Vector3 r;

            r.x = Mathf.Abs(1 - 3f * Mathf.Abs(a.x));
            r.y = Mathf.Abs(1 - 3f * Mathf.Abs(a.y));
            r.z = Mathf.Abs(1 - 3f * Mathf.Abs(a.z));

            float da = Mathf.Max(r.x,r.y);
            float db = Mathf.Max(r.y,r.z);
            float dc = Mathf.Max(r.z,r.x);
            float c = (float)(Mathf.Min(da,Mathf.Max(db,dc))-1.0f)/s;

            d = Mathf.Max(d,c);
        }

        return d;
    }
    
    float FCT_BBSK(Vector3 pos) {
        Vector3 cFcParams = new Vector3(2.18f, -0.18f, 0);
        Vector3 CSize = new Vector3(1.4f,0.87f, 1.1f);
        Vector3 p;
        p.x = 2 * pos.x;
        p.y = 2 * pos.z;
        p.z = 2 * pos.y;
        float scale = 1.0f;
    
        for( int i=0; i < 4;i++ ) 
        {
            p.x = 2 * Mathf.Clamp(p.x, -CSize.x, CSize.x) - p.x;
            p.y = 2 * Mathf.Clamp(p.y, -CSize.y, CSize.y) - p.y;
            p.z = 2 * Mathf.Clamp(p.z, -CSize.z, CSize.x) - p.z;
            float r2 = Vector3.Dot(p,p);
            //float r2 = dot(p,p+sin(p.z*.5)); //Alternate fractal
            float k = Mathf.Max((2f)/(r2), .17f);
            p *= k;
            //p *=rot;
            //p= p.yzx;
            p+=new Vector3(0.2f,0.2f,-0.5f);
            scale *= k;
        }
        
        p.x = 2 * Mathf.Clamp(p.x, -CSize.x * 4, CSize.x * 4) - p.x;
        p.y = 2 * Mathf.Clamp(p.y, -CSize.y * 4, CSize.y * 4) - p.y;
        p.z = 2 * Mathf.Clamp(p.z, -CSize.z * 4, CSize.x * 4) - p.z;

        for( int i=0; i < 8;i++ )
        {
            p.x = 2 * Mathf.Clamp(p.x, -CSize.x, CSize.x) - p.x;
            p.y = 2 * Mathf.Clamp(p.y, -CSize.y, CSize.y) - p.y;
            p.z = 2 * Mathf.Clamp(p.z, -CSize.z, CSize.x) - p.z;
            float r2 = Vector3.Dot(p,p);
            //float r2 = dot(p,p+sin(p.z*.3)); //Alternate fractal
            float k = Math.Max((cFcParams.x)/(r2),  0.027f);
            p     *= k;
            scale *= k;
            p.y += cFcParams.y;
        }

        float l = Mathf.Sqrt(Vector2.SqrMagnitude(new Vector2(p.x, p.y)));
        //l = mix(l,l2,0.5);
        float rxy = l - 4;
        float n = p.z;
        rxy = Math.Max(rxy, -(n) / 4);
        float dist = (rxy) / Math.Abs(scale);
        dist *= .75f;

        return dist;

    }*/


}
