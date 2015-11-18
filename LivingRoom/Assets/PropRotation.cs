using UnityEngine;
using System.Collections;
using System.Net;
using System.Net.Sockets;
using System.Text;

public class PropRotation : MonoBehaviour {
    public Kinematics data;//using the UDP data data received in the Kinematics class instantiation in the drone object. Does not require sending in data again,less expensive. 
	//CharacterController 
    GameObject thisprop;
	GameObject thisprop2;
	GameObject thisprop3;
	GameObject thisprop4;
 
   float[] angle = { 0, 0, 0, 0 };
    
    float C = 40;

   // double t1, t2, t3, t4;
	// Use this for initialization
	void Start () {

      //  QualitySettings.antiAliasing = 50;

        // if (client.Available > 0)
        //   {

        //      // Received bytes
        //      IPEndPoint anyIP = new IPEndPoint(IPAddress.Any, port);
        //       byte[] data = client.Receive(ref anyIP);
        //
        //      double t1 = System.BitConverter.ToDouble(data, 48);
        //     double t2 = System.BitConverter.ToDouble(data, 56);
        //     double t3 = System.BitConverter.ToDouble(data, 64);
        //     double t4 = System.BitConverter.ToDouble(data, 72);
        //
        // }
        //thisprop = GetComponent<CharacterController>();

        thisprop = GameObject.Find("rotor2");
		thisprop2 = GameObject.Find("rotor1");
		thisprop3 = GameObject.Find("rotor3");
		thisprop4 = GameObject.Find("rotor4");

		angle[0] = 0f;
        angle[1] = 0f;
        angle[2] = 0f;
        angle[3] = 0f;
      
	}
	
	// Update is called once per frame
	void Update () {
        //angle = angle + 98.7f;

        //using known value of Thrust to relate to prop speed. We know from Blade Elemant Theory that rotational speed is proportional to Sqrt of Thrust. 
       // angle[0] = angle[0] + C*Mathf.Sqrt(data.getThrustValues(1));
       // angle[1] = angle[1] + C*Mathf.Sqrt(data.getThrustValues(2));
       // angle[2] = angle[2] + C*Mathf.Sqrt(data.getThrustValues(3));
       // angle[3] = angle[3] + C*Mathf.Sqrt(data.getThrustValues(4));
        //Debug.Log(data.getThrustValues(1));

        for(int i=0;i<4;i++)
        {
            angle[i] = angle[i] + 128.7f;
        }
      
        thisprop.transform.eulerAngles = new Vector3(0f,(float)angle[0] ,-0f);
	    thisprop2.transform.eulerAngles = new Vector3(0f,(float)angle[1] ,-0f);
	    thisprop3.transform.eulerAngles = new Vector3(0f,(float)angle[2] ,-0f);
	    thisprop4.transform.eulerAngles = new Vector3(0f,(float)angle[3] ,-0f);

        
        


	}
}
