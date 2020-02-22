using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class ObjectSdfCollision : MonoBehaviour
{
    
	public float friction        = 0.3f;
    public float angularFriction = 0.6f;
    public float restitution     = 0.9f;

    private const float MAX_DIST = 10f;
    private const float MIN_DIST = 0.01f;

    public float STATIC_GRAVITY_MODIFIER  = 1.2f;
    public float BUERIED_GRAVITY_MODIFIER = 3f;

    private Rigidbody rigidbody_;
    
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
        result.isBuried  = DistanceFunction.CalcDistance(sphere.transform.position - 0.975f * radius * dir) < MIN_DIST;
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
		    Collision(rigidbody_.velocity.normalized, sphere, sphere.GetComponent<SphereCollider>().radius * this.transform.localScale.x);
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
            /*Debug.DrawLine(transform.position, ray.position, Color.black);
            Debug.DrawLine(ray.position, ray.position + ray.normal, Color.yellow);*/
            return true;
        }
        else
        {
            //Debug.DrawLine(transform.position, ray.position, Color.blue);
            return false;
        }
    }

	// Use this for initialization
    void Start()
    {
        rigidbody_ = GetComponent<Rigidbody>();
    }
}
