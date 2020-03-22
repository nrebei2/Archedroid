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
    public Fractal2 fractal;
    
    [Header("Texture and Lighting (Ambient light can be configured in lighting settings)")]
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

    [Header("Parameters")] 
    public GameObject ParamChanger;
    public static Vector3 Params;
    public Vector3 Parameters;
    
    
    [Header("ModInterval")]
    public Vector3 _modInterval;


    private string enumAsString;
    private string lastEnumAsString;
    
    
    public Texture texturePattern;
    public Texture textureColormap;

    [Header("cTexParam")]
    public float _patternScale;
    public float _colormapYscale;
    public float _patternInfluence;
    public float _fractalInfluence;
    private Vector4 cTexParam;
    
    
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
        
    }

    void Update()
    {
        foreach(MeshRenderer mat in GetComponents<MeshRenderer>())
        {
            mat.material = material;
        }
        
        _power = power;

        if (ParamChanger.activeSelf)
        {
            Parameters.x = LinearParamX.currentLinearMapping;
            Parameters.y = LinearParamY.currentLinearMapping;
            Parameters.z = LinearParamZ.currentLinearMapping;
        }

        Params = Parameters;
        _fractalScale = fractalScale;
        _fractaldegreeRotate = fractaldegreeRotate;
        _fractalPosition = fractalPosition;
        
        rend.sharedMaterial.SetVector("_modInterval", _modInterval);
        
        rend.sharedMaterial.SetVector("_fractalPosition", _fractalPosition);
        rend.sharedMaterial.SetFloat("_fractalSmooth", _fractalSmooth);
        rend.sharedMaterial.SetVector("_fractaldegreeRotate", _fractaldegreeRotate);
        rend.sharedMaterial.SetVector("_fractalScale", _fractalScale);
        rend.sharedMaterial.SetFloat("_power", _power);
        rend.sharedMaterial.SetVector("Params", Params);
        
        cTexParam = new Vector4(_patternScale, _colormapYscale, _patternInfluence, _fractalInfluence);
        rend.sharedMaterial.SetVector("cTexParam", cTexParam);

        
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
        rend.sharedMaterial.SetTexture("_TexturePattern", texturePattern);
        rend.sharedMaterial.SetTexture("_TextureColormap", textureColormap);
        rend.sharedMaterial.SetColor("_Color", Color);
        rend.sharedMaterial.SetFloat("_Metallic", Metallic);
        rend.sharedMaterial.SetFloat("_Glossiness", Smoothness);

        foreach (Fractal2 suit in (Fractal2[]) Enum.GetValues(typeof(Fractal2)))
        {
            if (fractal == suit)
            {
                _fractalNumber = (float) suit;
                enumAsString = Enum.GetName (typeof(Fractal2), (int)_fractalNumber);
                Debug.Log(enumAsString);
            }
        }
        
        if (lastEnumAsString != enumAsString)
        {
            material.DisableKeyword(lastEnumAsString);
            material.EnableKeyword(enumAsString);
            Debug.Log(lastEnumAsString + "--->" + enumAsString);
            lastEnumAsString = enumAsString;
        }

        rend.sharedMaterial.SetFloat("_fractalNumber", _fractalNumber);
        //Debug.Log(_fractalNumber);
    }
}
    public enum Fractal2{ StaticMandelbulb, DynamicMandelbulb, Julia, Juliabulb, Sierpinski, Mandelbox, 
        KaleidoscopicIFS, Tglad, Hartverdrahtet, PseudoKleinian, PseudoKnightyan, Mandelbulb2, MengerSponge, apo, plane,
        FCT_BBSK, trinoise, RecursiveTetrahedron, TruchetTentacles, FCT_PROTEIN, FCT_ORBIT, FCT_MNMT, FCT_CRAB, FCT_HUB, 
        FCT_HYPERAPO, FCT_DLBT, FCT_MZGN, FCT_PIPES, FCT_APOP, FCT_APO, FCT_HTVT, FCT_KNKL, FCT_KIFS, FCT_TEXT, FCT_TEST }