using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using Valve.VR;
using Valve.VR.InteractionSystem;
public class rigidbodycontroller : MonoBehaviour
{

    public SteamVR_Action_Vector2 input;
    
    public float Speed = 1f;
    public float JumpHeight = 2f;
    public LayerMask Ground;
    
    public float GroundDistance = 0.2f;
    private Rigidbody _body;
    private bool _isGrounded = true;
    
    public Transform _groundChecker;

    void Start()
    {
        _body = GetComponent<Rigidbody>();
    }

    void Update()
    {
        Vector3 direction = Player.instance.hmdTransform.TransformDirection(new Vector3(input.axis.x, 0, input.axis.y));
        _body.MovePosition(_body.position + Vector3.ProjectOnPlane(direction, Vector3.up) * Speed * Time.fixedDeltaTime);

        
        _isGrounded = Physics.CheckSphere(_groundChecker.position, GroundDistance, Ground, QueryTriggerInteraction.Ignore) || true;
        
        if (Input.GetButtonDown("Jump") && _isGrounded)
        {
            _body.AddForce(Vector3.up * Mathf.Sqrt(JumpHeight * -2f * Physics.gravity.y), ForceMode.VelocityChange);
        }
        
    }
}