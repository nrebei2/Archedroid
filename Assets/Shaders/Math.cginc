#ifndef MATH_CGINC
#define MATH_CGINC
#define PI 3.14159265358979

float rand(float3 r)
{
    return frac(sin(dot(r.xy, r.yz)));
}

inline float Mod(float a, float b)
{
    return frac(abs(a / b)) * abs(b);
}

inline float2 Mod(float2 a, float2 b)
{
    return frac(abs(a / b)) * abs(b);
}

inline float3 Mod(float3 a, float3 b)
{
    return frac(abs(a / b)) * abs(b);
}

inline float SmoothMin(float d1, float d2, float k)
{
    float h = exp(-k * d1) + exp(-k * d2);
    return -log(h) / k;
}

inline float Repeat(float pos, float span)
{
    if(span==0.0) return pos;
    return Mod(pos, span) - span * 0.5;
}

inline float2 Repeat(float2 pos, float2 span)
{
    return Mod(pos, span) - span * 0.5;
}

inline float3 Repeat(float3 pos, float3 span)
{

    return Mod(pos, span) - span * 0.5;
}

inline float3 Rotate(float3 p, float angle, float3 axis)
{
    float3 a = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float r = 1.0 - c;
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

inline float3 TwistY(float3 p, float power)
{
    float s = sin(power * p.y);
    float c = cos(power * p.y);
    float3x3 m = float3x3(
          c, 0.0,  -s,
        0.0, 1.0, 0.0,
          s, 0.0,   c
    );
    return mul(m, p);
}

inline float3 TwistX(float3 p, float power)
{
    float s = sin(power * p.y);
    float c = cos(power * p.y);
    float3x3 m = float3x3(
        1.0, 0.0, 0.0,
        0.0,   c,   s,
        0.0,  -s,   c
    );
    return mul(m, p);
}

inline float3 TwistZ(float3 p, float power)
{
    float s = sin(power * p.y);
    float c = cos(power * p.y);
    float3x3 m = float3x3(
          c,   s, 0.0,
         -s,   c, 0.0,
        0.0, 0.0, 1.0
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

// Union
float4 opUS(float4 d1, float4 d2, float k) {
	float h = clamp(0.5 + 0.5 * (d2.w - d1.w) / k, 0.0, 1.0);
	// Lerp between color and distance between objects
	float3 color = lerp(d2.rgb, d1.rgb, h);
	float dist =  lerp(d2.w, d1.w, h) - k * h * (1.0 - h);
	return float4(color, dist);
}

float4 qsqr( in float4 a ) // square a quaterion
{
    return float4( a.x*a.x - a.y*a.y - a.z*a.z - a.w*a.w,
                 2.0*a.x*a.y,
                 2.0*a.x*a.z,
                 2.0*a.x*a.w );
}

void sphereFold(inout float3 z, inout float dz)
{
	float r2 = dot(z,z);
	if (r2 < 0.5)
    { 
		float temp = 2.0;
		z *= temp;
		dz*= temp;
	}
    else if (r2 < 1.0)
    { 
		float temp = 1.0 / r2;
		z *= temp;
		dz*= temp;
	}
}

void boxFold(inout float3 z, inout float dz)
{
	z = clamp(z, -1.0, 1.0) * 2.0 - z;
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
    float rad = 0.0174532925 * angle;
    float s, c;
    sincos(rad, s, c);
    return float3(p.x, c*p.y + s*p.z, -s*p.y + c*p.z);
}
float3 RotateY(float3 v, float degree)
{
	float rad = 0.0174532925 * degree;
	float cosY = cos(rad);
	float sinY = sin(rad);
	return float3(cosY * v.x - sinY * v.z, v.y, sinY * v.x + cosY * v.z);
}
float3 RotateZ(float3 p, float angle)
{
    float rad = 0.0174532925 * angle;
    float s, c;
    sincos(rad, s, c);
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

float4 _OffsetPosition;
float4 _Scale;

float3 localize(float3 p)
{
    return mul(unity_WorldToObject, float4(p, 1)).xyz * _Scale.xyz; // + _OffsetPosition.xyz;
}

float fract(float x) {
 return x - floor(x);
}

float fract(float3 x) {
 return x - floor(x);
}

float mix(float x, float y, float z) {
 return x*(1. - z) + y*z;
}

void domainRep3(inout float3 p, float3 size) {
	p = modc(p + size,2. * size) - size;
}

float3 domainRep3Idx(inout float3 p, float3 size) {
	float3 idx = floor((p + size)/size/2.);
	p = modc(p + size, 2. * size) - size;
	return idx;
}

float smin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float smax( in float a, in float b, in float s ){
    float h = clamp( 0.5 + 0.5*(a-b)/s, 0.0, 1.0 );
    return mix(b, a, h) + h*(1.0-h)*s;
}

float hash(float h) {
	return fract(sin(h) * 43758.5453123);
}

float noise3d(float3 x) {
	float3 p = floor(x);
	float3 f = fract(x);
	f = f * f * (3.0 - 2.0 * f);

	float n = p.x + p.y * 157.0 + 113.0 * p.z;
	return mix(
			mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
					mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
			mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
					mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
}

#endif
