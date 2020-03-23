using System.Collections;
using System.Collections.Generic;
using Valve.VR;
using Valve.VR.InteractionSystem;
using UnityEngine;

public class SwitchMovement : MonoBehaviour
{
    
    public SteamVR_Action_Boolean switchAction;
    private bool timeSwitch = false;

    // Update is called once per frame
    void Update()
    {
        if (switchAction.stateDown)
        {
            Switch();
        }
        /*if (timeSwitch) { GetComponent<Rigidbody>().Sleep(); }
        else { GetComponent<Rigidbody>().WakeUp(); }*/
    }

    private void Switch()
    {
         timeSwitch = !timeSwitch;
         //Switch to flying when timeswitch becomes true
         GetComponent<PlayerFlyingController>().enabled = timeSwitch;
         GetComponent<CharacterController>().enabled = timeSwitch;
         GetComponent<Rigidbody>().isKinematic = timeSwitch;
         
         //Switch to flying when timeswitch becomes false
         GetComponent<PlayerController>().enabled = !timeSwitch;
         GetComponent<MoveColliderWithCamera>().enabled = !timeSwitch;
         GetComponent<CapsuleCollider>().enabled = !timeSwitch;

         
    }
}
