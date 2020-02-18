using System;
using UnityEngine;
using UnityEngine.Rendering;
#if UNITY_EDITOR
using UnityEditor;
#endif
using System.Collections.Generic;
using Valve.VR.InteractionSystem;

[ExecuteInEditMode]
public class RaymarchingCamera : MonoBehaviour
{
    // Testing....
    [SerializeField] private Shader _shader;
    
    [Header("Setup")]
    [Range(1,1000)]
    public float Loop;
    [Range(0.0001f,1f)]
    public float MinimumDistance;
    [Range(0.001f, 2.0f)]
    public float DistanceMultiplier;

    /*
    [Header("Shading")]
    [Range(0, 4)]
    public float _ShadowIntensity;
    public Vector2 _ShadowDistance;
    [Range(1,128)]
    public float _ShadowPenumbra;

    [Header("Ambient Occlusion")]
    [Range(0.01f, 10.0f)]
    public float _AoStepsize;
    [Range(1,5)]
    public int _AoIterations;
    [Range(0,1)]
    public float _AoInstensity;

    [Header("Reflection")]
    [Range(0, 2)]
    public int _ReflectionCount;
    [Range(0, 1)]
    public float _ReflectionIntensity;
    [Range(0, 1)]
    public float _EnvReflIntenisty;

    public Cubemap _ReflectionCube;
    */

    
    [Header("Signed Distance Field")]
    public Fractal fractal;
    public Texture fractalTexutre;
    public Color Color;
    public float Metallic;
    public float Smoothness;
    
    public static float _fractalNumber;

    public static Vector3 _fractalPosition;
    public Vector3 fractalPosition;
    
    public float _fractalSmooth;
    public static Vector3 _fractaldegreeRotate;
    public Vector3 fractaldegreeRotate;
    
    public static Vector3 _fractalScale;
    public Vector3 fractalScale;
    
    public static float _power;
    public float power;
    
    public static Vector3 Params;
    public Vector3 Parameters;
    

    [Header("Make Object Repeat Indefinitely on Axis (Value 1 is True, Any Other Value is False)")]
    public Vector3 _modBool;
    
    [Header("ModInterval")]
    public Vector3 _modInterval;
    
    public Material material
    {
        get
        {
            if (!_raymarchMat && _shader)
            {
                _raymarchMat = new Material(_shader);
            }
            return _raymarchMat;
        }
    }

    private Material _raymarchMat;
    Renderer rend;
    private void Start()
    {
        rend = GetComponent<Renderer> ();
        foreach(MeshRenderer mat in GetComponents<MeshRenderer>())
        {
            mat.material = material;
        }
    }

    void Update()
    {
        _power = power;
        Parameters.x = LinearParamX.currentLinearMapping;
        Parameters.y = LinearParamY.currentLinearMapping;
        Parameters.z = LinearParamZ.currentLinearMapping;
        Params = Parameters;
        _fractalScale = fractalScale;
        _fractaldegreeRotate = fractaldegreeRotate;
        _fractalPosition = fractalPosition;

        rend.sharedMaterial.SetVector("_modBool", _modBool);
        rend.sharedMaterial.SetVector("_modInterval", _modInterval);
        
        rend.sharedMaterial.SetVector("_fractalPosition", _fractalPosition);
        rend.sharedMaterial.SetFloat("_fractalSmooth", _fractalSmooth);
        rend.sharedMaterial.SetVector("_fractaldegreeRotate", _fractaldegreeRotate);
        rend.sharedMaterial.SetVector("_fractalScale", _fractalScale);
        rend.sharedMaterial.SetFloat("_power", _power);
        rend.sharedMaterial.SetVector("Params", Params);
        
        /*
        material.SetFloat("_ShadowIntensity", _ShadowIntensity);
        material.SetFloat("_ShadowPenumbra", _ShadowPenumbra);
        material.SetVector("_ShadowDistance", _ShadowDistance);
        
        material.SetFloat("_AoStepsize", _AoStepsize);
        material.SetFloat("_AoInstensity", _AoInstensity);
        material.SetInt("_AoIterations", _AoIterations);

        material.SetInt("_ReflectionCount", _ReflectionCount);
        material.SetFloat("_ReflectionIntensity", _ReflectionIntensity);
        material.SetFloat("_EnvReflIntenisty", _EnvReflIntenisty);
        material.SetTexture("_ReflectionCube", _ReflectionCube);
        */
        
        rend.sharedMaterial.SetFloat("_Loop", Loop);
        rend.sharedMaterial.SetFloat("_MinDistance", MinimumDistance);
        rend.sharedMaterial.SetFloat("_DistanceMultiplier", DistanceMultiplier);
        rend.sharedMaterial.SetTexture("_Texture", fractalTexutre);
        rend.sharedMaterial.SetColor("_Color", Color);
        rend.sharedMaterial.SetFloat("_Metallic", Metallic);
        rend.sharedMaterial.SetFloat("_Smoothness", Smoothness);

        
        switch (fractal)
        {
            case Fractal.StaticMandelbulb:
                _fractalNumber = 0;
                break;
            case Fractal.DynamicMandelbulb:
                _fractalNumber = 1;
                break;
            case Fractal.Julia:
                _fractalNumber = 2;
                break;
            case Fractal.Juliabulb:
                _fractalNumber = 3;
                break;
            case Fractal.Sierpinski:
                _fractalNumber = 4;
                break;
            case Fractal.Mandelbox:
                _fractalNumber = 5;
                break;
            case Fractal.KaleidoscopicIFS:
                _fractalNumber = 6;
                break;
            case Fractal.Tglad:
                _fractalNumber = 7;
                break;
            case Fractal.Hartverdrahtet:
                _fractalNumber = 8;
                break;
            case Fractal.PseudoKleinian:
                _fractalNumber = 9;
                break;
            case Fractal.PseudoKnightyan:
                _fractalNumber = 10;
                break;
            case Fractal.Mandelbulb2:
                _fractalNumber = 11;
                break;
            case Fractal.MengerSponge:
                _fractalNumber = 12;
                break;
            case Fractal.apo:
                _fractalNumber = 13;
                break;
            case Fractal.plane:
                _fractalNumber = 14;
                break;
            case Fractal.FCT_BBSK:
                _fractalNumber = 15;
                break;
            case Fractal.trinoise:
                _fractalNumber = 16;
                break;
            case Fractal.RecursiveTetrahedron:
                _fractalNumber = 17;
                break;
            case Fractal.TruchetTentacles:
                _fractalNumber = 18;
                break;
            case Fractal.FCT_PROTEIN:
                _fractalNumber = 19;
                break;
            case Fractal.FCT_ORBIT:
                _fractalNumber = 20;
                break;
            case Fractal.FCT_MNMT:
                _fractalNumber = 21;
                break;
            case Fractal.FCT_CRAB:
                _fractalNumber = 22;
                break;
            case Fractal.FCT_HUB:
                _fractalNumber = 23;
                break;
            case Fractal.FCT_HYPERAPO:
                _fractalNumber = 24;
                break;
            case Fractal.FCT_DLBT:
                _fractalNumber = 25;
                break;
            case Fractal.FCT_MZGN:
                _fractalNumber = 26;
                break;
            case Fractal.FCT_PIPES:
                _fractalNumber = 27;
                break;
            case Fractal.FCT_APOP:
                _fractalNumber = 28;
                break;
            case Fractal.FCT_APO:
                _fractalNumber = 29;
                break;
            case Fractal.FCT_HTVT:
                _fractalNumber = 30;
                break;
            case Fractal.FCT_KNKL:
                _fractalNumber = 31;
                break;
            case Fractal.FCT_KIFS:
                _fractalNumber = 32;
                break;
            case Fractal.FCT_TEXT:
                _fractalNumber = 33;
                break;
            case Fractal.FCT_TEST:
                _fractalNumber = 34;
                break;
        }
        
        rend.sharedMaterial.SetFloat("_fractalNumber", _fractalNumber);
        //Debug.Log(_fractalNumber);
        

        

    }
}
    public enum Fractal{ StaticMandelbulb, DynamicMandelbulb, Julia, Juliabulb, Sierpinski, Mandelbox, 
        KaleidoscopicIFS, Tglad, Hartverdrahtet, PseudoKleinian, PseudoKnightyan, Mandelbulb2, MengerSponge, apo, plane,
        FCT_BBSK, trinoise, RecursiveTetrahedron, TruchetTentacles, FCT_PROTEIN, FCT_ORBIT, FCT_MNMT, FCT_CRAB, FCT_HUB, 
        FCT_HYPERAPO, FCT_DLBT, FCT_MZGN, FCT_PIPES, FCT_APOP, FCT_APO, FCT_HTVT, FCT_KNKL, FCT_KIFS, FCT_TEXT, FCT_TEST }