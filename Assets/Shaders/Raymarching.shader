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

#define WORLD_SPACE

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

uniform float4 _fractal;
uniform float _fractalSmooth;
uniform float3 _fractaldegreeRotate;

uniform float3 _modInterval;

uniform float3 _modBool;

uniform float _fractalNumber;

uniform float _power;


// @block DistanceFunction
inline float DistanceFunction(float3 p)
{
    //TODO: Fix optimization, specifially in this part of the code
    
    p = RotateZ(RotateY(RotateX(p - _fractal.xyz, _fractaldegreeRotate.x), _fractaldegreeRotate.y), _fractaldegreeRotate.z);
    
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
            return FCT_BBSK(p);     
            break;                 
        default:
            return apo(p, .0274, float3(1., 1., 1.3), float3(0., 0., 0.));
            break;         
    }
    return 1;
    
    //return apo(p, .0274, float3(1., 1., 1.3), float3(0., 0., 0.));
}
// @endblock


sampler2D _Texture;
float4 _Texture_ST;
// @block PostEffect
inline void PostEffect(RaymarchInfo ray, inout PostEffectOutput o)
{
    // find where exactly in the code brings up the black areas around the edges (Increasing "Loop" fixed that, though...)
    float ao = 1.0 - 1.0 * ray.loop / ray.maxLoop;
    o.Occlusion *= ao;
    o.Emission += tex2D(_Texture, ray.endPos.xy * _Texture_ST.xy + _Texture_ST.zw);
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