/*
 * To change position
 * 
 * transform.localPosition = new Vector3(0, -5, 0);
 * 
 * To control Velocity
 * transform.Translate(Vx,Vy,Vz);
 * 
 * transform.Translate(moveSpeed*Input.GetAxis ("Horizontal")*Time.deltaTime,0f,moveSpeed*Input.GetAxis ("Vertical")*Time.deltaTime);
 * */


using UnityEngine;
using System.Collections;
using System.Net;
using System.Net.Sockets;
using System.Text;



public class Kinematics : MonoBehaviour {

	public float moveSpeed; 
	public float rotSpeed; 
	CharacterController uav ;
	Vector3 Forward;
	Vector3 Side;
	Vector3 Up;
	UdpClient client; //UDP Client Port
	public int port;
	double x;
	double y;
	double z;
	double p;
	double q;
	double r;
    double t1, t2, t3, t4;//declare thrusts for each rotor
	// infos
	public string lastReceivedUDPPacket="";
    AudioSource quad_sound;//new audiosource for quad's sounds.
    double velocity_squared;//will be summation of x,y,and z velocities
    double totalThrust;//will be average of total thrust,for now
    // Use this for initialization
    void Start () 
	{
		moveSpeed = 1f;
		rotSpeed = 1f*180/3.1415f;
		uav = GetComponent<CharacterController>();
		Forward = transform.TransformDirection(Vector3.forward);
		Side = transform.TransformDirection(Vector3.right);
		Up = transform.TransformDirection(Vector3.up);
		transform.localPosition = new Vector3(11f, 2.5f, 16.5f);
		//transform.localPosition = new Vector3(0f, 0f, 0f);
		print("UDPSend.init()");
		port = 25000;

        //Start of my edit
        quad_sound = (AudioSource)gameObject.AddComponent<AudioSource>();
        AudioClip staticSound;//create an audio clip variable that will take the audioclip file we need for the quad

        //staticSound = (AudioClip)Resources.Load("heli_500_100%");
       staticSound = (AudioClip)Resources.Load("Quadrotor_MidThrottle_firstTake");
        quad_sound.clip = staticSound;//the file "heli_500_100%" is now a clip in the AudioSource quad_sound
        quad_sound.loop = true; //enable looping of sound.
       
        //end of edit



        client = new UdpClient (port);
        quad_sound.Play();//play sound, hopefully should be static, contionously looping sound

    }

	// Update is called once per frame
	void Update () 
	{
		x = 0;
		y = 0;
		z = 0;
		p = 0;
		q = 0;
		r = 0;
       
			if (client.Available > 0) {
	
			// Received bytes
			IPEndPoint anyIP = new IPEndPoint (IPAddress.Any, port);
			byte[] data = client.Receive (ref anyIP);
			
			// Bytes encode the UTF8 encoding in the text format .
			//string text = Encoding.UTF8.GetString (data);
			
			// DDisplay the retrieved text .
			//print (">> " + text);
			
			// latest UDPpacket
			//lastReceivedUDPPacket = text;
			
			//x = (double)((float)System.BitConverter.ToInt32 (data, 0)/1024.0f);
			x = System.BitConverter.ToDouble (data, 0);
			y = System.BitConverter.ToDouble (data, 8);
			z = System.BitConverter.ToDouble (data, 16);
			p = System.BitConverter.ToDouble (data, 24);
			q = System.BitConverter.ToDouble (data, 32);
			r = System.BitConverter.ToDouble (data, 40);
           // t1 = System.BitConverter.ToDouble(data, 48);
          //  t2 = System.BitConverter.ToDouble(data, 56);
           // t3 = System.BitConverter.ToDouble(data, 64);
           // t4 = System.BitConverter.ToDouble(data, 72);
			//double.TryParse(text, out x);
			//print(">> x=" + x.ToString());
            
            
		}
		// For old quad
		//uav.Move(moveSpeed*((float)x*Forward + (float)y*Side - (float)z*Up)*Time.deltaTime);

		//For new quad
		uav.Move(moveSpeed*((float)x*Forward - (float)y*Side + (float)z*Up)*Time.deltaTime);
        velocity_squared = (x * x) + (y * y) + (z * z);//find summ of velocities squared
        //totalThrust = (t1 + t2 + t3 + t4)/4;//plain  summation for now
        quad_sound.volume = 0.7f+(float)(0.1*velocity_squared);
        quad_sound.pitch=0.9f+(float)(0.6 * Mathf.Sqrt((float)velocity_squared));
      
        //DEBUGGING//
       // Debug.Log(velocity_squared);
       //Debug.Log(quad_sound.pitch);
      // Debug.Log(quad_sound.volume);

		//transform.Rotate (rotSpeed * (float)q*Time.deltaTime, (float)r*rotSpeed*Time.deltaTime, -rotSpeed*(float)p*Time.deltaTime, Space.Self);
		//transform.localPosition = moveSpeed * ((float)x * Forward + (float)y * Side + ((float)z) * Up);

		//transform.eulerAngles = new Vector3((float)q*rotSpeed,(float)r*rotSpeed,-(float)p*rotSpeed);
		transform.eulerAngles = new Vector3(-(float)p*rotSpeed,-(float)r*rotSpeed,(float)q*rotSpeed);
        //transform.Translate = new Vector3((float)x, (float)y, (float)z);
        //transform.Translate((float)x, (float)y, (float)z);

       // quad_sound.Play();//play sound, hopefully should be static, contionously looping sound
    }


    //public float getThrustValues(int i)
    //{
    //    if (i == 1) { return (float)t1; }      
    //        
    //    else if (i == 2) { return (float)t2; }
          
    //    else if (i == 3) { return (float)t3; }
           
     //   else if (i == 4) { return (float)t4; }

    //    else { return 0; }
        
       
     // }
}
