using System.Collections;
using System.Collections.Generic;
using UnityEngine;



namespace Valve.VR.InteractionSystem
{
    public class LinearParamY : MonoBehaviour
    {
        
        public LinearMapping linearMapping;
        public static float currentLinearMapping = 0f;
        
        // Start is called before the first frame update
        void Start()
        {

        }

        // Update is called once per frame
        void Update()
        {
            currentLinearMapping = linearMapping.value;
        }
    }
}
