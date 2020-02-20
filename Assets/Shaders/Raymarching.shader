Shader "Raymarching/Forward_ModWorld_for_VR"
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

uniform float3 _modBool;

uniform float _fractalNumber;

uniform float _power;
uniform float3 Params;


// @block DistanceFunction
inline float DistanceFunction(float3 p)
{
    //TODO: Fix optimization, specifially in this part of the code
    
    p = RotateZ(RotateY(RotateX(p - _fractalPosition, _fractaldegreeRotate.x), _fractaldegreeRotate.y), _fractaldegreeRotate.z);
    p = ScaleX(p, _fractalScale.x);
    p = ScaleY(p, _fractalScale.y);
    p = ScaleZ(p, _fractalScale.z);
    
    if (_modBool.x == 1) {
        p.x = Repeat(p.x, _modInterval.x);
    }
    if (_modBool.y == 1) {
        p.y = Repeat(p.y, _modInterval.y);
    }
    if (_modBool.z == 1) {
        p.z = Repeat(p.z, _modInterval.z);
    }


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
            return apo(p, .0274, float3(1., 1., 1.3), float3(0., 0., 0.));
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
            break;  */
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
            return apo(p, .0274, float3(1., 1., 1.3), float3(0., 0., 0.));
            break;         
    }
    return 1;
    
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
    
    // find where exactly in the code brings up the black areas around the edges (Increasing "Loop" fixed that, though...)
    float ao = 1.0 - 1.0 * ray.loop / ray.maxLoop;
    o.Occlusion *= ao;
    o.Albedo.rgb = tex2D(_Texture, ray.endPos.xy * _Texture_ST.xy + _Texture_ST.zw);
    o.Emission *= ao;
    
    
    
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