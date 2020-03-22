Shader "Raymarching/AllFractals"
{

Properties
{
    [Header(PBS)]
    _Color("Color", Color) = (1.0, 1.0, 1.0, 1.0)
    _Metallic("Metallic", Range(0.0, 1.0)) = 1
    _Glossiness("Smoothness", Range(0.0, 1.0)) = 0

    [Header(Pass)]
    [Enum(UnityEngine.Rendering.CullMode)] _Cull("Culling", Int) = 0

    [Toggle][KeyEnum(Off, On)] _ZWrite("ZWrite", Float) = 1

    [Header(Raymarching)]
    _Loop("Loop", Range(1, 1000)) = 100
    _MinDistance("Minimum Distance", Range(0.0001, 0.1)) = 0.0005
    _DistanceMultiplier("Distance Multiplier", Range(0.001, 2.0)) = 1.237

// @block Properties
[Header(Additional Parameters)]
_Texture("Texture", 2D) = "" {}
_TexturePattern("TexturePattern", 2D) = "" {}
_TextureColormap("TextureColormap", 2D) = "" {}
// @endblock
}

SubShader
{

Tags
{
    "RenderType" = "Opaque"
    "Queue" = "Geometry"
    "DisableBatching" = "True"
}

Cull [_Cull]

CGINCLUDE
#pragma multi_compile StaticMandelbulb DynamicMandelbulb Julia Juliabulb Sierpinski Mandelbox KaleidoscopicIFS Tglad Hartverdrahtet PseudoKleinian PseudoKnightyan Mandelbulb2 MengerSponge apo plane FCT_BBSK trinoise RecursiveTetrahedron TruchetTentacles FCT_PROTEIN FCT_ORBIT FCT_MNMT FCT_CRAB FCT_HUB FCT_HYPERAPO FCT_DLBT FCT_MZGN FCT_PIPES FCT_APOP FCT_APO FCT_HTVT FCT_KNKL FCT_KIFS FCT_TEXT FCT_TEST
//#define WORLD_SPACE

#define OBJECT_SHAPE_NONE

#define CAMERA_INSIDE_OBJECT

#define USE_RAYMARCHING_DEPTH

#define DISABLE_VIEW_CULLING

#define SPHERICAL_HARMONICS_PER_PIXEL

#define DISTANCE_FUNCTION DistanceFunction
#define PostEffectOutput SurfaceOutputStandard
#define POST_EFFECT PostEffect

#include "Assets\Shaders/Common.cginc"
#include "Assets/Shaders\DistanceFunctions.cginc"

uniform int _MaxIterations;
uniform float _Accuracy;
uniform float _maxDistance;

uniform fixed4 _GroundColor;
uniform float _ColorIntensity;

uniform float2 _ShadowDistance;
uniform float _ShadowIntensity, _ShadowPenumbra;

uniform int _ReflectionCount;
uniform float _ReflectionIntensity;
uniform float _EnvReflIntenisty;
uniform samplerCUBE _ReflectionCube;

uniform float3 _fractalPosition;
uniform float _fractalSmooth;
uniform float3 _fractaldegreeRotate;
uniform float3 _fractalScale;

uniform float3 _modInterval;

uniform float _fractalNumber;

uniform float _power;
uniform float3 Params;

uniform float4 cTexParam;

float4 distance;


// @block DistanceFunction
inline float DistanceFunction(float3 p)
{
    //TODO: Fix optimization, specifially in this part of the code
    
    p = RotateZ(RotateY(RotateX(p - _fractalPosition, _fractaldegreeRotate.x), _fractaldegreeRotate.y), _fractaldegreeRotate.z);
    p = ScaleX(p, _fractalScale.x);
    p = ScaleY(p, _fractalScale.y);
    p = ScaleZ(p, _fractalScale.z);
    
    p.x = Repeat(p.x, _modInterval.x);
    p.y = Repeat(p.y, _modInterval.y);
    p.z = Repeat(p.z, _modInterval.z);

    distance = sdfmap(p, _power, Params);
    
    return distance.x;
    
    //return apo(p, .0274, float3(1., 1., 1.3), float3(0., 0., 0.));
}
// @endblock


// Glow pattern needed functions
/*
float2 pattern(float2 p)
{
    p = frac(p);
    float r = 0.123;
    float v = 0.0, g = 0.0;
    r = frac(r * 9184.928);
    float cp, d;
    
    d = p.x;
    g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 1000.0);
    d = p.y;
    g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 1000.0);
    d = p.x - 1.0;
    g += pow(clamp(3.0 - abs(d), 0.0, 1.0), 1000.0);
    d = p.y - 1.0;
    g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 10000.0);

    const int ITER = 12;
    for(int i = 0; i < ITER; i ++)
    {
        cp = 0.5 + (r - 0.5) * 0.9;
        d = p.x - cp;
        g += pow(clamp(1.0 - abs(d), 0.0, 1.0), 200.0);
        if(d > 0.0) {
            r = frac(r * 4829.013);
            p.x = (p.x - cp) / (1.0 - cp);
            v += 1.0;
        }
        else {
            r = frac(r * 1239.528);
            p.x = p.x / cp;
        }
        p = p.yx;
    }
    v /= float(ITER);
    return float2(g, v);
}

inline float3 _GetDistanceFunctionNormal(float3 p)
{
    const float d = 1e-3;
    return normalize( float3(
        DistanceFunction(p+float3(  d,0.0,0.0))-DistanceFunction(p+float3( -d,0.0,0.0)),
        DistanceFunction(p+float3(0.0,  d,0.0))-DistanceFunction(p+float3(0.0, -d,0.0)),
        DistanceFunction(p+float3(0.0,0.0,  d))-DistanceFunction(p+float3(0.0,0.0, -d)) ));
}
*/


sampler2D _Texture;
float4 _Texture_ST;

sampler2D _TexturePattern;
sampler2D _TextureColormap;

float3 normal;
float glow = 0.0;

// @block PostEffect
inline void PostEffect(RaymarchInfo ray, inout PostEffectOutput o)
{
    // TODO: Fix and implement glow pattern (will use it during specific parts of gameplay)
    // glow pattern effect from Unity5Effects
    /*
    normal = _GetDistanceFunctionNormal(ray.endPos);
    float3 p3 = localize(ray.endPos);
    p3 *= 2.0;
    glow += max((modc(length(p3) - _Time.y*3, 15.0) - 12.0)*0.7, 0.0);
    float2 p2 = pattern(p3.xz*0.5);
    if(p2.x<1.3) { glow = 0.0; }
    glow += max(1.0-abs(dot(-GetCameraForward(), normal)) - 0.4, 0.0) * 1.0;
    float3 emission = float3(0.7, 0.7, 1.0)*glow*0.6;
    o.Emission += emission;
    */
    
    // find where exactly in the code brings up the black areas around the edges ---> its the raymarching algoritm itself
    float ao = 1.0 - 1.0 * ray.loop / ray.maxLoop;
    o.Occlusion *= ao;
    o.Emission *= ao;
    //o.Albedo.rgb = tex2D(_Texture, ray.endPos.xy * _Texture_ST.xy + _Texture_ST.zw);
    float3 txt = tex2D(_TexturePattern, ray.endPos.xz * cTexParam.x).rgb;
    float4 col = tex2D(_TextureColormap, fractalTexMap(ray.endPos,txt,distance.y,distance.z,distance.w, cTexParam));
    //col *= 1.05 - 10. * max(0.3-0.2,0.);
    o.Albedo.rgb = col;
    
    
    // Makes colors go sicko mode
    //o.Albedo = abs(1.0 - Mod(ray.endPos * 0.1 + _Time.y, 2.0));
   
}   

// @endblock

ENDCG

Pass
{
    Tags { "LightMode" = "ForwardBase" }

    ZWrite [_ZWrite]

    CGPROGRAM
    #include "Assets\Shaders/ForwardBaseStandard.cginc"
    #pragma target 3.0
    #pragma vertex Vert
    #pragma fragment Frag
    #pragma multi_compile_instancing
    #pragma multi_compile_fog
    #pragma multi_compile_fwdbase
    ENDCG
}

Pass
{
    Tags { "LightMode" = "ForwardAdd" }
    ZWrite Off 
    Blend One One

    CGPROGRAM
    #include "Assets\Shaders/ForwardAddStandard.cginc"
    #pragma target 3.0
    #pragma vertex Vert
    #pragma fragment Frag
    #pragma multi_compile_instancing
    #pragma multi_compile_fog
    #pragma skip_variants INSTANCING_ON
    #pragma multi_compile_fwdadd_fullshadows
    ENDCG
}

}

Fallback Off

CustomEditor "uShaderTemplate.MaterialEditor"

}