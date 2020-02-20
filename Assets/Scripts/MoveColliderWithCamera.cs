using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MoveColliderWithCamera : MonoBehaviour
{
    [Header("First object is head, second is body, third is feet")]
    public GameObject[] sphereCollisions;
    
    private Rigidbody rigidbody_;
    
    public GameObject _groundChecker;

    [Header("Change for VR vs. PC")]
    public GameObject Camera;

    public GameObject PlayerHeadlight;
    
    // Start is called before the first frame update
    void Start()
    {
        rigidbody_ = GetComponent<Rigidbody>();
    }

    // Update is called once per frame
    void Update()
    {
        GetComponent<CapsuleCollider>().center = new Vector3(Camera.transform.localPosition.x, Camera.transform.localPosition.y - 0.5f,
            Camera.transform.localPosition.z);
        
        _groundChecker.transform.position = new Vector3(Camera.transform.position.x, Camera.transform.position.y - 1.5f,
            Camera.transform.position.z);

        float offset = 0f;
        for (int i = 0; i < 3; i++)
        {
            sphereCollisions[i].transform.localPosition = new Vector3(Camera.transform.localPosition.x, Camera.transform.localPosition.y - offset,
                Camera.transform.localPosition.z);
            offset += 0.5f;
        }

        PlayerHeadlight.transform.position = Camera.transform.position;
        PlayerHeadlight.transform.rotation = Camera.transform.localRotation;
    }
}
