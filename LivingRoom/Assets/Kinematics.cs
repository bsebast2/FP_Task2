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
	// infos
	public string lastReceivedUDPPacket="";


	// Use this for initialization
	void Start () 
	{
		moveSpeed = 1f;
		rotSpeed = 180/3.1415f;
		uav = GetComponent<CharacterController>();
		Forward = transform.TransformDirection(Vector3.forward);
		Side = transform.TransformDirection(Vector3.right);
		Up = transform.TransformDirection(Vector3.up);
		transform.localPosition = new Vector3(11f, 2.5f, 16.5f);
		//transform.localPosition = new Vector3(0f, 0f, 0f);
		print("UDPSend.init()");
		port = 25000;




			client = new UdpClient (port);
		
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
			//double.TryParse(text, out x);
			//print(">> x=" + x.ToString());
		}
		uav.Move(moveSpeed*((float)x*Forward + (float)y*Side - (float)z*Up)*Time.deltaTime);
		//transform.Rotate (rotSpeed * (float)q*Time.deltaTime, (float)r*rotSpeed*Time.deltaTime, -rotSpeed*(float)p*Time.deltaTime, Space.Self);
		//transform.localPosition = moveSpeed * ((float)x * Forward + (float)y * Side + ((float)z) * Up);

		//transform.eulerAngles = new Vector3((float)q*rotSpeed,(float)r*rotSpeed,-(float)p*rotSpeed);
		transform.eulerAngles = new Vector3((float)p*rotSpeed,-(float)r*rotSpeed,-(float)q*rotSpeed);
		//transform.Translate = new Vector3((float)x, (float)y, (float)z);
		//transform.Translate((float)x, (float)y, (float)z);
	}
}
