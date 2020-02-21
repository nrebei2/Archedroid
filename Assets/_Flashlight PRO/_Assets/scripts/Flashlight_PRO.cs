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
using Valve.VR;
using Valve.VR.InteractionSystem;

public class Flashlight_PRO : MonoBehaviour 
{
	[Space(10)]
	[SerializeField()] GameObject Lights; // all light effects and spotlight
	[SerializeField()] AudioSource switch_sound; // audio of the switcher
	[SerializeField()] ParticleSystem dust_particles; // dust particles



	private Light spotlight;
	private Material ambient_light_material;
	private Color ambient_mat_color;
	private bool is_enabled = false;

	public SteamVR_Action_Boolean turnOnAction;
	public SteamVR_Action_Single powerAction;
	private Interactable interactable;

	[SerializeField]
	float friction        = 0.3f;
    float angularFriction = 0.6f;
    float restitution     = 0.0f;

    private const float MAX_DIST = 10f;
    private const float MIN_DIST = 0.001f;

    public float STATIC_GRAVITY_MODIFIER  = 1.2f;
    private float BUERIED_GRAVITY_MODIFIER = 3f;

    private const int RaysToShoot = 30;
    float angle = 0;
    
    private Rigidbody rigidbody_;
    
    [Header("First object is head, second is body, third is feet")]
    public GameObject[] sphereCollisions;
    
    struct RaymarchingResult
    {
        public int     loop;
        public bool    isBuried;
        public float   distance;
        public float   length;
        public Vector3 direction;
        public Vector3 position;
        public Vector3 normal;
    }

    RaymarchingResult Raymarching(Vector3 dir, GameObject sphere, float radius)
    {
        var dist = 0f;
        var len  = 0f;
        var pos  = sphere.transform.position + radius * dir;
        var loop = 0;

        for (loop = 0; loop < 10; ++loop) {
            dist = DistanceFunction.CalcDistance(pos);
            len += dist;
            pos += dir * dist;
            if (dist < MIN_DIST || len > MAX_DIST) break;
        }

        var result = new RaymarchingResult();

        result.loop      = loop;
        result.isBuried  = DistanceFunction.CalcDistance(sphere.transform.position) < MIN_DIST;
        result.distance  = dist;
        result.length    = len;
        result.direction = dir;
        result.position  = pos;
        result.normal    = DistanceFunction.CalcNormal(pos);

        return result;
    }

    void FixedUpdate()
    {
	    foreach (GameObject sphere in sphereCollisions)
	    {
		    Collision(rigidbody_.velocity.normalized, sphere, sphere.GetComponent<SphereCollider>().radius / 10f);
		    Debug.Log(sphere.GetComponent<SphereCollider>().radius / 10f);
	    }
    }

    private bool Collision (Vector3 dir, GameObject sphere, float radius)
    {
        var ray = Raymarching(dir, sphere, radius);
        var v = rigidbody_.velocity;
        var g = Physics.gravity;

        if (ray.isBuried) {
            rigidbody_.AddForce((rigidbody_.mass * g.magnitude * BUERIED_GRAVITY_MODIFIER) * ray.normal);
            Debug.Log("buried!");
            return true;
        } else if (ray.length < MIN_DIST) {
            var prod = Vector3.Dot(v.normalized, ray.normal);
            var vv = (prod * v.magnitude) * ray.normal;
            var vh = v - vv;
            rigidbody_.velocity = vh * (1f - friction) + (-vv * restitution);
            // TODO: try implementing AddForceAtLocation
            rigidbody_.AddForce(-rigidbody_.mass * STATIC_GRAVITY_MODIFIER * g);
	        rigidbody_.AddTorque(-rigidbody_.angularVelocity * (1f - angularFriction));
            Debug.DrawLine(transform.position, ray.position, Color.black);
            Debug.DrawLine(ray.position, ray.position + ray.normal, Color.yellow);
            return true;
        }
        else
        {
            Debug.DrawLine(transform.position, ray.position, Color.blue);
            return false;
        }
    }

	// Use this for initialization
	void Start () 
	{
		rigidbody_ = GetComponent<Rigidbody>();
		// cache components
		spotlight = Lights.transform.Find ("Spotlight").GetComponent<Light> ();
		ambient_light_material = Lights.transform.Find ("ambient").GetComponent<Renderer> ().material;
		//ambient_mat_color = ambient_light_material.GetColor ("_TintColor");
		interactable = GetComponent<Interactable>();
	}

	/// <summary>
	/// changes the intensivity of lights from 0 to 100.
	/// call this from other scripts.
	/// </summary>
	private void Update()
	{
		if (interactable.attachedToHand != null)
		{
			float triggervalue = 100 * powerAction.GetAxis(SteamVR_Input_Sources.Any);
			Change_Intensivity(triggervalue);
			print(triggervalue);
		}
		
		if (interactable.attachedToHand != null && turnOnAction.stateDown)
		{
			Switch();
		}
	}

	public void Change_Intensivity(float percentage)
	{
		percentage = Mathf.Clamp (percentage, 0, 100);
		spotlight.intensity = (1.5f * percentage) / 100;
		//ambient_light_material.SetColor ("_TintColor", new Color(ambient_mat_color.r , ambient_mat_color.g , ambient_mat_color.b , percentage/2000)); 
	}




	/// <summary>
	/// switch current state  ON / OFF.
	/// call this from other scripts.
	/// </summary>
	public void Switch()
	{
		is_enabled = !is_enabled; 

		Lights.SetActive (is_enabled);

		if (switch_sound != null)
			switch_sound.Play ();
	}





	/// <summary>
	/// enables the particles.
	/// </summary>
	public void Enable_Particles(bool value)
	{
		if(dust_particles != null)
		{
			if(value)
			{
				dust_particles.gameObject.SetActive(true);
				dust_particles.Play();
			}
			else
			{
				dust_particles.Stop();
				dust_particles.gameObject.SetActive(false);
			}
		}
	}


}
