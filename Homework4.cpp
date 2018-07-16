#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
using namespace std;
#include <GL/glut.h>

#define PI 3.1415926

typedef struct{
	int num_vertices;
	int num_faces;
	float **vertex;
	int **face;
	float Kd, Ks, Exp_N, Odr, Odb, Odg;
	int n_edges;
	bool *IsFront;
    float **VertexNormal;
	float **VertexColor;
}Model;

GLsizei ww, wh;
float Iar, Iag, Iab;
float Ip[3];
float e_posx, e_posy, e_posz;
float L_posx[3], L_posy[3], L_posz[3];
float BGC_r, BGC_g, BGC_b; 
int num_light = 0;
int num_obj = 0;
int obj_i = -1;
float zNear, zFar, hFov;
float MM[4][4];
float EM[4][4];
float PM[4][4];
float **ZBuffer;
vector<float> **CBuffer;
Model *object;

void plot_pixel(int x, int y, float r, float g, float b){

	glColor3f(r, g, b);
	glBegin(GL_POINTS);	  
	  glVertex2i(x, y);
	glEnd();
 
}

void print_matrix(){
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			cout << setw(15) << MM[i][j] ;
		}
		cout << endl;
	}
}

void reset(){
	memset(MM, 0, sizeof(MM));
	for(int i=0; i<4; i++)
	MM[i][i]=1;
}

void translation(float x, float y, float z){
	float trans[4][4];
	float temp[4][4];
	memcpy(temp, MM, sizeof(temp));
	memset(trans, 0, sizeof(trans));

	 trans[0][3]=x;
	 trans[1][3]=y;
	 trans[2][3]=z;

	for(int i=0; i<4; i++)
		trans[i][i]=1;

	for(int i=0; i<4;i++){
		for(int j=0; j<4; j++){
      	   MM[i][j] = trans[i][0]*temp[0][j] + trans[i][1]*temp[1][j] + trans[i][2]*temp[2][j] + trans[i][3]*temp[3][j]; 
		}
	}
}

void scaling(float x, float y, float z){
	float scale[4][4];
	float temp[4][4];

	memcpy(temp, MM, sizeof(temp));
	memset(scale, 0, sizeof(scale));

	scale[0][0]=x;	scale[1][1]=y;	scale[2][2]=z;	scale[3][3]=1;

	for(int i=0; i<4;i++){
		for(int j=0; j<4; j++){
      	   MM[i][j] = scale[i][0]*temp[0][j] + scale[i][1]*temp[1][j] + scale[i][2]*temp[2][j] + scale[i][3]*temp[3][j]; 
		}
	}
}

void rotation(float angle_x, float angle_y, float angle_z){
	
	float rotate_x[4][4];
	float rotate_y[4][4];
	float rotate_z[4][4];
    float temp[4][4];

	if( angle_z != 0 ){
		float theta_z = angle_z*PI/180;

		memcpy(temp, MM, sizeof(temp));
    	memset(rotate_z, 0, sizeof(rotate_z));
		rotate_z[0][0]=cos(theta_z);  rotate_z[0][1]=-sin(theta_z);
		rotate_z[1][0]=sin(theta_z);  rotate_z[1][1]= cos(theta_z);
		rotate_z[2][2]= 1;
		rotate_z[3][3]= 1;

		for(int i=0; i<4;i++){
			for(int j=0; j<4; j++){
      		   MM[i][j] = rotate_z[i][0]*temp[0][j] + rotate_z[i][1]*temp[1][j] + rotate_z[i][2]*temp[2][j] + rotate_z[i][3]*temp[3][j]; 
			}
		}
	}

	if( angle_x != 0 ){
		float theta_x = angle_x*PI/180;

		memcpy(temp, MM, sizeof(temp));
    	memset(rotate_x, 0, sizeof(rotate_x));

	    rotate_x[1][1]=cos(theta_x);  rotate_x[1][2]=-sin(theta_x);
		rotate_x[2][1]=sin(theta_x);  rotate_x[2][2]= cos(theta_x);
		rotate_x[0][0]= 1;
		rotate_x[3][3]= 1;
		
		for(int i=0; i<4;i++){
			for(int j=0; j<4; j++){
      		   MM[i][j] = rotate_x[i][0]*temp[0][j] + rotate_x[i][1]*temp[1][j] + rotate_x[i][2]*temp[2][j] + rotate_x[i][3]*temp[3][j]; 
			}
		}
	}

	if( angle_y != 0 ){
		float theta_y = angle_y*PI/180;

		memcpy(temp, MM, sizeof(temp));
    	memset(rotate_y, 0, sizeof(rotate_y));

	    rotate_y[0][0]= cos(theta_y);  rotate_y[0][2]=sin(theta_y);
		rotate_y[2][0]=-sin(theta_y);  rotate_y[2][2]=cos(theta_y);
		rotate_y[1][1]= 1;
		rotate_y[3][3]= 1;
		
		for(int i=0; i<4;i++){
			for(int j=0; j<4; j++){
      		   MM[i][j] = rotate_y[i][0]*temp[0][j] + rotate_y[i][1]*temp[1][j] + rotate_y[i][2]*temp[2][j] + rotate_y[i][3]*temp[3][j]; 
			}
		}
	}

}

void init_buffer()
{
		// Initialize Z buffer array
		CBuffer = new vector<float> *[wh]; 
		ZBuffer = new float *[wh];
		for(int i=0; i<wh; i++){
			CBuffer[i] = new vector<float>[ww];
			ZBuffer[i] = new float[ww];
			for(int k=0; k<ww; k++){
				ZBuffer[i][k] = 5000;
				CBuffer[i][k].push_back(BGC_r);
				CBuffer[i][k].push_back(BGC_g);
				CBuffer[i][k].push_back(BGC_b);
			}
		}	
}

void Illumination()
{
	object[obj_i].VertexColor = new float *[object[obj_i].num_vertices];
	for(int i=0; i<object[obj_i].num_vertices; i++)
		object[obj_i].VertexColor[i] = new float[3];

	object[obj_i].VertexNormal = new float *[object[obj_i].num_vertices];
	for(int i=0; i<object[obj_i].num_vertices; i++)
		object[obj_i].VertexNormal[i] = new float[3];

	object[obj_i].IsFront = new bool[object[obj_i].num_faces];

	// ********** Compute normal vectors of the vertices ************
	vector<int> *face_store = new vector<int> [object[obj_i].num_vertices];   // store the near faces of each vetex

	float **FaceNormal = new float *[object[obj_i].num_faces];    // every normal vector of the face
	for(int i=0; i<object[obj_i].num_faces; i++)
		FaceNormal[i] = new float[3];

	for(int i=0; i<object[obj_i].num_faces; i++){
		// store the face which the vertex belong to
		for(int j=0; j<object[obj_i].n_edges; j++){
			int v_i;
			v_i = object[obj_i].face[i][j];   // v_i : index of the face 
			face_store[v_i-1].push_back(i);    			
		}
		// compute the normal vector for each face
		int f;
		float p1[3], p2[3], p3[3];
		for(int j=0; j<3; j++){         // read 3 vertices to compute the plane
			f = object[obj_i].face[i][j] -1;    
			if(j==0){
				p1[0] = object[obj_i].vertex[f][0];   p1[1] = object[obj_i].vertex[f][1];  p1[2] = object[obj_i].vertex[f][2];
			}
			else if(j==1){
				p2[0] = object[obj_i].vertex[f][0];   p2[1] = object[obj_i].vertex[f][1];  p2[2] = object[obj_i].vertex[f][2];
			}
			else if(j==2){
				p3[0] = object[obj_i].vertex[f][0];   p3[1] = object[obj_i].vertex[f][1];  p3[2] = object[obj_i].vertex[f][2];
			}
		}
			float V12[3], V13[3];
			float n[3], len_n;
			V12[0] =  p2[0] - p1[0];   V12[1] =  p2[1] - p1[1];   V12[2] =  p2[2] - p1[2]; 
			V13[0] =  p3[0] - p1[0];   V13[1] =  p3[1] - p1[1];   V13[2] =  p3[2] - p1[2]; 

			n[0] =  V12[1]*V13[2] - V12[2]*V13[1];    // n = v12 X v13
			n[1] =  V12[2]*V13[0] - V12[0]*V13[2];
			n[2] =  V12[0]*V13[1] - V12[1]*V13[0];

			len_n = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );   // normalization
			for( int j=0; j<3; j++ ){
				if(len_n != 0)
					FaceNormal[i][j] = -n[j]/len_n;
				else
					FaceNormal[i][j] = 0;
			}

            // check whether it's back face
			float eye_vector[] = { e_posx-p1[0], e_posy-p1[1], e_posz-p1[2] };	
            float NE;
			 NE = FaceNormal[i][0]*eye_vector[0] + FaceNormal[i][1]*eye_vector[1] + FaceNormal[i][2]*eye_vector[2];
			if( NE<0 )
				object[obj_i].IsFront[i] = false;
			else 
				object[obj_i].IsFront[i] = true;
	}


	for(int i=0; i<object[obj_i].num_vertices; i++){
		vector<int>::iterator v_nearf = face_store[i].begin(); // contain near face indices of the vertex
		float N[] = { 0, 0, 0 };  
		
		for(unsigned int j=0; j<face_store[i].size(); j++){  // sum up all the normal vectors of faces	
			for(int  k=0; k<3; k++){
				N[k] += FaceNormal[ *v_nearf ][k];
			}
			v_nearf++;
		}

		for(int  k=0; k<3; k++){
			N[k] = N[k] / face_store[i].size();
		}

		float N_len = sqrt( N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
		for(int  k=0; k<3; k++){      // Normalization
			if(N_len!=0)
				object[obj_i].VertexNormal[i][k] = N[k] / N_len;
			else
				object[obj_i].VertexNormal[i][k] = 0;
		}
	//	cout << object[obj_i].VertexNormal[i][0] << " " << object[obj_i].VertexNormal[i][1] <<" "<< object[obj_i].VertexNormal[i][2] <<" "<< endl << endl;
	}

	//  *************************   Illumination   ************************
	float Ka = 1.0, Fatt = 1.0;
	for(int i=0; i<object[obj_i].num_vertices; i++){
		float V[3], L[3], N[3], H[3], NL, RV,  L_len, V_len, H_len;

		object[obj_i].VertexColor[i][0] = Ka * Iar * object[obj_i].Odr;
		object[obj_i].VertexColor[i][1] = Ka * Iag * object[obj_i].Odg; 
		object[obj_i].VertexColor[i][2] = Ka * Iab * object[obj_i].Odb;
		
		//view vector
		V[0] = (e_posx - object[obj_i].vertex[i][0]);   V[1] = (e_posy - object[obj_i].vertex[i][1]);  V[2] = (e_posz - object[obj_i].vertex[i][2]);
		V_len = sqrt( V[0]*V[0] + V[1]*V[1] + V[2]*V[2] );
		for(int j=0; j<3; j++){      // Normalization
				if(V_len != 0)
					V[j] = V[j] / V_len;
				else
					V[j] = 0;
		}
		//vertex normal
		N[0] = object[obj_i].VertexNormal[i][0];   N[1] = object[obj_i].VertexNormal[i][1];  N[2] =  object[obj_i].VertexNormal[i][2]; 

		for(int li=0; li<num_light; li++){  // compound light
			//Light vector
			L[0] = (L_posx[li] - object[obj_i].vertex[i][0]);   L[1] = (L_posy[li] - object[obj_i].vertex[i][1]);  L[2] = (L_posz[li] - object[obj_i].vertex[i][2]); 
			L_len = sqrt( L[0]*L[0] + L[1]*L[1] + L[2]*L[2] );
			for(int j=0; j<3; j++){      // Normalization
				if(L_len != 0)
					L[j] = L[j] / L_len; 
				else
					L[j] = 0;
			}

			NL =  N[0]*L[0] + N[1]*L[1] + N[2]*L[2];
			if(NL<0)
				NL=0;

			H_len = sqrt( (L[0]+V[0])*(L[0]+V[0]) + (L[1]+V[1])*(L[1]+V[1]) + (L[2]+V[2])*(L[2]+V[2]) );
			for(int j=0; j<3; j++){
				if(H_len!=0){
					H[j] = (V[j]+L[j]) / H_len;
				}
				else
					H[j] = 0;
			}

			RV = H[0]*N[0] + H[1]*N[1] + H[2]*N[2];  // R.V is approximate to H.N 
			if(RV<0)
				RV=0;

			object[obj_i].VertexColor[i][0] += ( Fatt * object[obj_i].Kd * Ip[li] * NL ) * object[obj_i].Odr + Fatt * object[obj_i].Ks * Ip[li] * pow( RV, object[obj_i].Exp_N ); 
			object[obj_i].VertexColor[i][1] += ( Fatt * object[obj_i].Kd * Ip[li] * NL ) * object[obj_i].Odg + Fatt * object[obj_i].Ks * Ip[li] * pow( RV, object[obj_i].Exp_N) ; 
			object[obj_i].VertexColor[i][2] += ( Fatt * object[obj_i].Kd * Ip[li] * NL ) * object[obj_i].Odb + Fatt * object[obj_i].Ks * Ip[li] * pow( RV, object[obj_i].Exp_N ); 
		}
		//cout << object[obj_i].VertexColor[i][0] << "  "  << object[obj_i].VertexColor[i][1] << "  " << object[obj_i].VertexColor[i][2] << endl;
	}
	delete[](face_store);

	for(int i=0; i<object[obj_i].num_faces; i++)
		delete[](FaceNormal[i]);
	delete[](FaceNormal);
}

void sort_array(float **poly_vert, int i_max)
{
	float *temp_x = new float[object[obj_i].n_edges];   
	float *temp_y = new float[object[obj_i].n_edges];   
	float *temp_z = new float[object[obj_i].n_edges];   
		
			for(int j=0; j<object[obj_i].n_edges; j++){
				temp_x[j] = poly_vert[j][0];   
				temp_y[j] = poly_vert[j][1];   
				temp_z[j] = poly_vert[j][2];   
			}	
			int i_right, i_left;
			if( i_max - 1 < 0 )
				i_left = object[obj_i].n_edges-1;
			else
				i_left = i_max - 1;

			if( i_max + 1 >= object[obj_i].n_edges )
				i_right = 0;
			else
				i_right = i_max + 1 ;
	
			if( (temp_y[i_left]==temp_y[i_max]) && (temp_x[i_left]<temp_x[i_max]) ){	// check horizontal line
				for(int j=0; j<3; j++){
					poly_vert[j][0] = temp_x[i_left];    poly_vert[j][1] = temp_y[i_left];   poly_vert[j][2] = temp_z[i_left];			
					poly_vert[1][0] = temp_x[i_max];     poly_vert[1][1] = temp_x[i_max];    poly_vert[1][2] = temp_x[i_max];		
					poly_vert[2][0] = temp_x[i_right];   poly_vert[2][1] = temp_x[i_right];  poly_vert[2][2] = temp_x[i_right];
				}
				if( object[obj_i].n_edges == 4){
					if(i_left-1>=0){
						poly_vert[3][0] = temp_x[i_left-1];    poly_vert[3][1] = temp_y[i_left-1];     poly_vert[3][2] = temp_z[i_left-1]; 
					}
					else{
						poly_vert[3][0] = temp_x[object[obj_i].n_edges-1];    poly_vert[3][1] = temp_y[object[obj_i].n_edges-1];     poly_vert[3][2] = temp_z[object[obj_i].n_edges-1];
					}
				}
			}
			else{
				poly_vert[0][0] = temp_x[i_max];            poly_vert[0][1] = temp_y[i_max];           poly_vert[0][2] = temp_z[i_max];  //start point of scanline
				poly_vert[1][0] = temp_x[i_right];          poly_vert[1][1] = temp_y[i_right];         poly_vert[1][2] = temp_z[i_right];
				poly_vert[object[obj_i].n_edges-1][0] = temp_x[i_left];   poly_vert[object[obj_i].n_edges-1][1] = temp_y[i_left];  poly_vert[object[obj_i].n_edges-1][2] = temp_z[i_left];
				if( object[obj_i].n_edges == 4){
					if(i_left-1>=0){
						poly_vert[2][0] = temp_x[i_left-1];    poly_vert[2][1] = temp_y[i_left-1];     poly_vert[2][2] = temp_z[i_left-1]; 
					}
					else{
						poly_vert[2][0] = temp_x[object[obj_i].n_edges-1];    poly_vert[2][1] = temp_y[object[obj_i].n_edges-1];     poly_vert[2][2] = temp_z[object[obj_i].n_edges-1];
					}
				}
			}
			delete [](temp_z);   
			delete [](temp_y);   
        	delete [](temp_x);  
}

void zbuffer_shading() 
{
	    float **rgb_color = new float *[object[obj_i].n_edges];
		float **poly_vert = new float *[object[obj_i].n_edges];
		for(int i=0; i<object[obj_i].n_edges; i++){
			poly_vert[i] = new float[3];
		    rgb_color[i] = new float[3];
		}

		// calculate the depth of every polygon
		for(int i=0; i<object[obj_i].num_faces; i++){ 
  			for(int j=0; j<object[obj_i].n_edges; j++){
				int f;
				f = (object[obj_i].face[i][j]) - 1;  
				if( object[obj_i].IsFront[i]==false )
					break;
				else{
					poly_vert[j][0] = object[obj_i].vertex[f][0];       poly_vert[j][1] = object[obj_i].vertex[f][1];       poly_vert[j][2] = object[obj_i].vertex[f][2];
					rgb_color[j][0] = object[obj_i].VertexColor[f][0];	rgb_color[j][1] = object[obj_i].VertexColor[f][1];  rgb_color[j][2] = object[obj_i].VertexColor[f][2];
				}
			}

			// If it's hidden surface, skip this face 
			if( object[obj_i].IsFront[i] == false )
				continue;
            
			// find maximum y and its index to start scanline
			float ymax, ymin;
			float z_ymax, z_ymin, x_ymax, x_ymin;
			int i_min=0, i_max=0;
			for(int j=0; j<object[obj_i].n_edges; j++){
				if(j==0){
				   ymin = poly_vert[j][1];   x_ymin = poly_vert[j][0];  z_ymin = poly_vert[j][2];
				   ymax = poly_vert[j][1];   x_ymax = poly_vert[j][0];  z_ymax = poly_vert[j][2];
				}
				else{
				   if( poly_vert[j][1] < ymin ){
					  ymin = poly_vert[j][1];  x_ymin = poly_vert[j][0];  z_ymin = poly_vert[j][2]; 
					  i_min = j;
				   }
				   if( poly_vert[j][1] > ymax ){
					  ymax = poly_vert[j][1];  x_ymin = poly_vert[j][0];  z_ymax = poly_vert[j][2];
					  i_max = j;
				   }
				}
			}

			//sort the order of the vertices and colors
			sort_array(poly_vert, i_max);
			sort_array(rgb_color, i_max);
						
			// calculate the plane eq.
			float p1[3], p2[3], p3[3];   // find 3 vertices of the plane
			for(int j=0; j<3; j++){
			   p1[j] = poly_vert[0][j];
			   p2[j] = poly_vert[1][j];
			   p3[j] = poly_vert[2][j];
			}

			float V12[3], V13[3];
			float n[3], len_n;
			V12[0] =  p2[0] - p1[0];   V12[1] =  p2[1] - p1[1];   V12[2] =  p2[2] - p1[2]; 
			V13[0] =  p3[0] - p1[0];   V13[1] =  p3[1] - p1[1];   V13[2] =  p3[2] - p1[2]; 

			n[0] =  V12[1]*V13[2] - V12[2]*V13[1];    // n1 = v12 X v13
			n[1] =  V12[2]*V13[0] - V12[0]*V13[2];
			n[2] =  V12[0]*V13[1] - V12[1]*V13[0];

			len_n = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );   // normalization
			for( int j=0; j<3; j++ ){
				if(len_n!=0)
					n[j] = -n[j]/len_n;
				else
					n[j] = 0;
			}

            // find D of the eq.
			float D;
			D = -( n[0]*p1[0] + n[1]*p1[1] + n[2]*p1[2] );

			// scanline
			float xa, xb, za, zb; 
			float ra, ga, ba, rb, gb, bb;
			for(int ys= (int)ymax; ys>(int)ymin; ys--){
				if ( object[obj_i].n_edges == 4 ){
		   			if( ys > poly_vert[3][1] ){
						xa = (ys-poly_vert[3][1])/(poly_vert[0][1]-poly_vert[3][1]) * (poly_vert[0][0]-poly_vert[3][0]) + poly_vert[3][0];
						za = (ys-poly_vert[3][1])/(poly_vert[0][1]-poly_vert[3][1]) * (poly_vert[0][2]-poly_vert[3][2]) + poly_vert[3][2]; // Z-depth
						ra = (ys-poly_vert[3][1])/(poly_vert[0][1]-poly_vert[3][1]) * (rgb_color[0][0]-rgb_color[3][0]) + rgb_color[3][0]; // R 
						ga = (ys-poly_vert[3][1])/(poly_vert[0][1]-poly_vert[3][1]) * (rgb_color[0][1]-rgb_color[3][1]) + rgb_color[3][1]; // G
						ba = (ys-poly_vert[3][1])/(poly_vert[0][1]-poly_vert[3][1]) * (rgb_color[0][2]-rgb_color[3][2]) + rgb_color[3][2]; // B
					}
					else{
						xa = (ys-poly_vert[2][1])/(poly_vert[3][1]-poly_vert[2][1]) * (poly_vert[3][0]-poly_vert[2][0]) + poly_vert[2][0];
						za = (ys-poly_vert[2][1])/(poly_vert[3][1]-poly_vert[2][1]) * (poly_vert[3][2]-poly_vert[2][2]) + poly_vert[2][2]; // Z-depth
						ra = (ys-poly_vert[2][1])/(poly_vert[3][1]-poly_vert[2][1]) * (rgb_color[3][0]-rgb_color[2][0]) + rgb_color[2][0]; // R 
						ga = (ys-poly_vert[2][1])/(poly_vert[3][1]-poly_vert[2][1]) * (rgb_color[3][1]-rgb_color[2][1]) + rgb_color[2][1]; // G
						ba = (ys-poly_vert[2][1])/(poly_vert[3][1]-poly_vert[2][1]) * (rgb_color[3][2]-rgb_color[2][2]) + rgb_color[2][2]; // B
					}
				}
				else if( object[obj_i].n_edges == 3 ){
					if( ys > poly_vert[2][1] ){
						xa = (ys-poly_vert[2][1])/(poly_vert[0][1]-poly_vert[2][1]) * (poly_vert[0][0]-poly_vert[2][0]) + poly_vert[2][0];
						za = (ys-poly_vert[2][1])/(poly_vert[0][1]-poly_vert[2][1]) * (poly_vert[0][2]-poly_vert[2][2]) + poly_vert[2][2]; // Z-depth
						ra = (ys-poly_vert[2][1])/(poly_vert[0][1]-poly_vert[2][1]) * (rgb_color[0][0]-rgb_color[2][0]) + rgb_color[2][0]; // R 
						ga = (ys-poly_vert[2][1])/(poly_vert[0][1]-poly_vert[2][1]) * (rgb_color[0][1]-rgb_color[2][1]) + rgb_color[2][1]; // G
						ba = (ys-poly_vert[2][1])/(poly_vert[0][1]-poly_vert[2][1]) * (rgb_color[0][2]-rgb_color[2][2]) + rgb_color[2][2]; // B
					}
					else{
						xa = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (poly_vert[2][0]-poly_vert[1][0]) + poly_vert[1][0];
						za = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (poly_vert[2][2]-poly_vert[1][2]) + poly_vert[1][2]; // Z-depth
						ra = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (rgb_color[2][0]-rgb_color[1][0]) + rgb_color[1][0]; // R 
						ga = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (rgb_color[2][1]-rgb_color[1][1]) + rgb_color[1][1]; // G
						ba = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (rgb_color[2][2]-rgb_color[1][2]) + rgb_color[1][2]; // B
					}
				}			

				if( ys > poly_vert[1][1] ){
					xb = (ys-poly_vert[1][1])/(poly_vert[0][1]-poly_vert[1][1]) * (poly_vert[0][0]-poly_vert[1][0]) + poly_vert[1][0];
					zb = (ys-poly_vert[1][1])/(poly_vert[0][1]-poly_vert[1][1]) * (poly_vert[0][2]-poly_vert[1][2]) + poly_vert[1][2];  // Z-depth
					rb = (ys-poly_vert[1][1])/(poly_vert[0][1]-poly_vert[1][1]) * (rgb_color[0][0]-rgb_color[1][0]) + rgb_color[1][0];  // R 
					gb = (ys-poly_vert[1][1])/(poly_vert[0][1]-poly_vert[1][1]) * (rgb_color[0][1]-rgb_color[1][1]) + rgb_color[1][1];  // G
					bb = (ys-poly_vert[1][1])/(poly_vert[0][1]-poly_vert[1][1]) * (rgb_color[0][2]-rgb_color[1][2]) + rgb_color[1][2];  // B		
				}
				else {
					if( object[obj_i].n_edges == 4 ){
						xb = (ys-poly_vert[2][1])/(poly_vert[1][1]-poly_vert[2][1]) * (poly_vert[1][0]-poly_vert[2][0]) + poly_vert[2][0];
						zb = (ys-poly_vert[2][1])/(poly_vert[1][1]-poly_vert[2][1]) * (poly_vert[1][2]-poly_vert[2][2]) + poly_vert[2][2]; // Z-depth
					    rb = (ys-poly_vert[2][1])/(poly_vert[1][1]-poly_vert[2][1]) * (rgb_color[1][0]-rgb_color[2][0]) + rgb_color[2][0];  // R 
					    gb = (ys-poly_vert[2][1])/(poly_vert[1][1]-poly_vert[2][1]) * (rgb_color[1][1]-rgb_color[2][1]) + rgb_color[2][1];  // G
					    bb = (ys-poly_vert[2][1])/(poly_vert[1][1]-poly_vert[2][1]) * (rgb_color[1][2]-rgb_color[2][2]) + rgb_color[2][2];  // B
					}
					else{
						xb = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (poly_vert[2][0]-poly_vert[1][0]) + poly_vert[1][0];
						zb = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (poly_vert[2][2]-poly_vert[1][2]) + poly_vert[1][2]; // Z-depth
					    rb = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (rgb_color[2][0]-rgb_color[1][0]) + rgb_color[1][0];  // R 
					    gb = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (rgb_color[2][1]-rgb_color[1][1]) + rgb_color[1][1];  // G
					    bb = (ys-poly_vert[1][1])/(poly_vert[2][1]-poly_vert[1][1]) * (rgb_color[2][2]-rgb_color[1][2]) + rgb_color[1][2];  // B
					}
				}

				float rs, gs, bs, zs;
				for(int x=(int)xa; x<(int)xb; x++){
					if( (xb-xa)!=0 ){
						zs = (x-xa)/(xb-xa) * (zb-za) + za;
					    rs = (x-xa)/(xb-xa) * (rb-ra) + ra;  gs = (x-xa)/(xb-xa) * (gb-ga) + ga;  bs = (x-xa)/(xb-xa) * (bb-ba) + ba;
					}
					else{
						zs = za;  rs=ra;   gs=ga;  bs= ba;
					}

				    if( zs < ZBuffer[ys][x] ){
						ZBuffer[ys][x] = zs;
						CBuffer[ys][x].clear();
						CBuffer[ys][x].push_back(rs);	CBuffer[ys][x].push_back(gs);	CBuffer[ys][x].push_back(bs);
					}
				}
		}
	}
	
	for(int i=0; i<object[obj_i].n_edges; i++){
		delete[](poly_vert[i]);
		delete[](rgb_color[i]);
	}
	delete[](poly_vert);
	delete[](rgb_color);
}

void EyeSpace( float px, float py, float pz, float cx, float cy, float cz, float tilt ){
	float mirror[4][4];
	float GRM[4][4];
	float trans[4][4];

	float vz[3]={cx-px, cy-py, cz-pz};  // v3 = vz
	float vt[3]={0.0, 1.0, 0.0};
	float v1[3];
	float v2[3];
	float len_v1, len_v2, len_v3;

	len_v3 = sqrt( vz[0]*vz[0] + vz[1]*vz[1] + vz[2]*vz[2] );
    for( int i=0; i<3; i++ ){
		vz[i] = vz[i]/len_v3;
	}

	v1[0] =  vz[1]*vt[2]-vz[2]*vt[1];    // v1 = vt X vz
	v1[1] =  vz[2]*vt[0]-vz[0]*vt[2];
	v1[2] =  vz[0]*vt[1]-vz[1]*vt[0];

	len_v1 = sqrt( v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] );
	for( int i=0; i<3; i++ ){
		v1[i] = -v1[i]/len_v1;
	}

	v2[0] =  v1[1]*vz[2]-v1[2]*vz[1];    // v2 = v3 X v1
	v2[1] =  v1[2]*vz[0]-v1[0]*vz[2];
	v2[2] =  v1[0]*vz[1]-v1[1]*vz[0];

	len_v2 = sqrt( v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2] );
	for( int i=0; i<3; i++ ){
		v2[i] = -v2[i]/len_v2; 
	}

	for(int i=0; i<4; i++){            // generate rotation matrix
		for(int j=0; j<4; j++){
			if( i==0 && j<3 )
				GRM[i][j] = v1[j];
			else if( i==1 && j<3 )
				GRM[i][j] = v2[j];
			else if( i==2 && j<3 )
				GRM[i][j] = vz[j];
			else if( i==3 && j==3 )
				GRM[i][j] = 1.0;
			else
				GRM[i][j] = 0.0;
		}
	}

	memset(mirror, 0, sizeof(mirror));
	mirror[0][0]= -1;
	for(int i=1; i<4; i++)
		mirror[i][i]= 1;

	memset(trans, 0, sizeof(trans));
	trans[0][3]=-px;
	trans[1][3]=-py;
	trans[2][3]=-pz;
	for(int i=0; i<4; i++)
		trans[i][i]=1;

	for(int i=0; i<4;i++){
			for(int j=0; j<4; j++){
      		   EM[i][j] = GRM[i][0]*trans[0][j] + GRM[i][1]*trans[1][j] + GRM[i][2]*trans[2][j] + GRM[i][3]*trans[3][j]; 
			}
	}	

	float temp[4][4];
	memcpy(temp, EM, sizeof(temp));
	for(int i=0; i<4;i++){
			for(int j=0; j<4; j++){
      		   EM[i][j] = mirror[i][0]*temp[0][j] + mirror[i][1]*temp[1][j] + mirror[i][2]*temp[2][j] + mirror[i][3]*temp[3][j]; 
			}
	}	

	float rotate_z[4][4];
	if( tilt != 0 ){
		float theta = tilt*PI/180;

		memcpy(temp, EM, sizeof(temp));
    	memset(rotate_z, 0, sizeof(rotate_z));
		rotate_z[0][0]=cos(theta);  rotate_z[0][1]=-sin(theta);
		rotate_z[1][0]=sin(theta);  rotate_z[1][1]= cos(theta);
		rotate_z[2][2]= 1;
		rotate_z[3][3]= 1;

		for(int i=0; i<4;i++){
			for(int j=0; j<4; j++){
      		   EM[i][j] = rotate_z[i][0]*temp[0][j] + rotate_z[i][1]*temp[1][j] + rotate_z[i][2]*temp[2][j] + rotate_z[i][3]*temp[3][j]; 
			}
		}
	}

	cout << endl << "EM:"<< endl;
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			cout << setw(15) << EM[i][j] ;
		}
		cout << endl;
	}
    
}

void ProjSpace( float znear, float zfar, float fov ){

	float AR = ww/wh;
	//float AR = (view_right - view_left)/(view_up - view_bottom);
	float thetaF = fov*PI/180;

	memset(PM, 0, sizeof(PM));
	PM[0][0] = 1;
	PM[1][1] = AR;
	PM[2][2] = (zfar/(zfar-znear)) * tan(thetaF);
	PM[2][3] = ((zfar*znear)/(znear-zfar)) * tan(thetaF);
	PM[3][2] = tan(thetaF);

	cout << endl << "PM:"<< endl;
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			cout << setw(15) << PM[i][j] ;
		}
		cout << endl;
	}
	cout << endl;
}

void ScreenSpace(){
		
   	float tmp_vertex[4]; 
	float tmp_MM[4][4]; 

	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			tmp_MM[i][j] = MM[i][j];
		}
	}
	obj_i = 0;
	while( obj_i < num_obj ){        // object loop

		// Vertex Illumination
		Illumination();

		for(int i=0; i<object[obj_i].num_vertices; i++){      // v' = EM*v
				for(int j=0; j<3; j++)
					tmp_vertex[j] = object[obj_i].vertex[i][j]; 
				
				tmp_vertex[3] = 1;

				for(int j=0; j<4; j++)
					object[obj_i].vertex[i][j] = EM[j][0]*tmp_vertex[0] + EM[j][1]*tmp_vertex[1] + EM[j][2]*tmp_vertex[2] + EM[j][3]*tmp_vertex[3];  
		}

		for(int i=0; i<object[obj_i].num_vertices; i++){      // v' = PM*v
				for(int j=0; j<3; j++)
					tmp_vertex[j] = object[obj_i].vertex[i][j]; 
	
				tmp_vertex[3] = 1;

				for(int j=0; j<4; j++)
					object[obj_i].vertex[i][j] = PM[j][0]*tmp_vertex[0] + PM[j][1]*tmp_vertex[1] + PM[j][2]*tmp_vertex[2] + PM[j][3]*tmp_vertex[3];  
		}

		//Perspective Divide
		for(int i=0; i<object[obj_i].num_vertices; i++){      
			for(int j=0; j<4; j++)
				object[obj_i].vertex[i][j] = object[obj_i].vertex[i][j] / object[obj_i].vertex[i][3]; 
		}

		// WVM ( Utilize MM to produce WVM )
	    reset();
	    translation( -1, -1, 0 );
	    scaling( ww/2, wh/2, 1 );
	    translation( ww, wh, 0 );

		for(int i=0; i<object[obj_i].num_vertices; i++){      // v' = WVM*v
			for(int j=0; j<3; j++)
				tmp_vertex[j] = object[obj_i].vertex[i][j]; 
	
			tmp_vertex[3] = 1;

			for(int j=0; j<4; j++)
				object[obj_i].vertex[i][j] = MM[j][0]*tmp_vertex[0] + MM[j][1]*tmp_vertex[1] + MM[j][2]*tmp_vertex[2] + MM[j][3]*tmp_vertex[3];  
		}

		// Shading
		zbuffer_shading();

		obj_i++;
	}

		// restore MM
		for(int i=0; i<4; i++){
				for(int j=0; j<4; j++){
					MM[i][j] = tmp_MM[i][j];
				}
		}
}

void load_model(string modelname){
    
	ifstream modelfile(modelname);
    
	if( !modelfile ){
		cerr << "Model File can't be opened!!" << endl;
		exit(1);
	}
	modelfile >> object[obj_i].num_vertices >> object[obj_i].num_faces ;

	object[obj_i].vertex = new float *[object[obj_i].num_vertices];
	for(int i=0; i<object[obj_i].num_vertices; i++)
		object[obj_i].vertex[i] = new float[4];

	for(int i=0; i<object[obj_i].num_vertices; i++){
		modelfile >> object[obj_i].vertex[i][0] 
				  >> object[obj_i].vertex[i][1]
				  >> object[obj_i].vertex[i][2];
		object[obj_i].vertex[i][3] = 1;   // w
	}

	object[obj_i].face = new int *[object[obj_i].num_faces];
	for(int i=0; i<object[obj_i].num_faces; i++){
		modelfile >> object[obj_i].n_edges ;
		object[obj_i].face[i] = new int[object[obj_i].n_edges];

		for(int j=0; j<object[obj_i].n_edges; j++)
			modelfile >> object[obj_i].face[i][j];
	}
}

void load_file( string filename ){
	ifstream file_1( filename );
	ifstream file_2( filename );
	
	if( !file_1 || !file_2){
	  cerr << "File can't open!!" << endl;
	  exit(1);
	}

	string line;
	while( getline(file_1, line)){

		 if(line.substr(0, 5) == "light"){
			num_light++;
		 }
		 if(line.substr(0, 6) == "object"){
			 num_obj++;                 // count the number of objects 
		 }
		 else if(line.substr(0, 10) == "background"){
			istringstream strval(line.substr(11));
			strval >> BGC_r >> BGC_g >> BGC_b;  
		 }
	}

	reset();							//  reset model matrix
	object = new Model[num_obj];		//  create object arrary
	bool init_row = true;               //  first line for Screen Size

	while( getline(file_2, line) ){
		if( line.substr(0,1) == "#" )
			continue;
		else if( init_row ){
			istringstream strval(line.substr(0));
			strval >> ww >> wh;
			cout << "Screen width: "<< ww << endl
				 << "Screen height: " << wh << endl;
			init_buffer();
			init_row = false;
		}
		else if(line.substr(0, 7) == "ambient"){
            float Ia;
			istringstream strval(line.substr(8));
			strval >> Ia;
			Iar = Ia; Iag = Ia;  Iab = Ia;
			cout << "Ambient: "<< Ia << endl;
		}
		else if(line.substr(0, 5) == "light"){
			int light_i=0;
			istringstream strval(line.substr(6));
			strval >> light_i ;
			strval >> Ip[light_i-1]  >> L_posx[light_i-1] >> L_posy[light_i-1] >> L_posz[light_i-1];
			cout << "Light: "<< "No." << light_i  << "  Intensity: "<< Ip[light_i-1] << "  Postion: ( " << L_posx[light_i-1]  << ", " << L_posy[light_i-1] << ", "<< L_posz[light_i-1] << " )" << endl;
		}
		else if(line.substr(0, 5) == "reset"){      
		  cout << "Reset" << endl;
		  reset();
		  print_matrix();
		}
		else if(line.substr(0, 5) == "scale"){
			float scale_x, scale_y, scale_z;
			istringstream strval(line.substr(6));
			strval >> scale_x;
			strval >> scale_y;
			strval >> scale_z;
			cout << "Scale: ( "<< scale_x <<", "<< scale_y << ", "<< scale_z <<" )" << endl;
			scaling( scale_x, scale_y, scale_z );
			print_matrix();
		}
		else if(line.substr(0, 9) == "translate"){
			float mov_x, mov_y, mov_z;
			istringstream strval(line.substr(10));
			strval >> mov_x;		
			strval >> mov_y;
			strval >> mov_z;
			cout << "Translate: ( "<< mov_x<<", "<< mov_y << ", " << mov_z << " )" << endl;
			translation(mov_x, mov_y, mov_z);
			print_matrix();
		}
		else if(line.substr(0, 6) == "rotate"){
			float angle_x, angle_y, angle_z;
			istringstream strval(line.substr(7));
			strval >> angle_x;
			strval >> angle_y;
			strval >> angle_z;
			cout << "Rotate: ( "<< angle_x << ", " << angle_y << ", " <<angle_z <<" )" << endl;
			if(angle_x >= 360)
				angle_x = (int)angle_x % 360;
			else if(angle_x < 0)
				angle_x += 360;

			if(angle_y >= 360)
				angle_y = (int)angle_y % 360;
			else if(angle_y < 0)
				angle_y += 360;

			if(angle_z >= 360)
				angle_z = (int)angle_z % 360;
			else if(angle_z < 0)
				angle_z += 360;

			rotation( angle_x, angle_y, angle_z );
			print_matrix();

		}
		else if(line.substr(0, 6) == "object"){
			obj_i++;
			
			string obj_file;
			istringstream strval(line.substr(7));
			strval >> obj_file
			       >> object[obj_i].Odr >> object[obj_i].Odg >> object[obj_i].Odb
		           >> object[obj_i].Kd  >> object[obj_i].Ks  >> object[obj_i].Exp_N;
			
			load_model(obj_file);
			
			float tmp_vertex[4]; 
			for(int i=0; i<object[obj_i].num_vertices; i++){
				for(int j=0; j<3; j++){
					tmp_vertex[j] = object[obj_i].vertex[i][j]; 
				}
				tmp_vertex[3] = 1;

				for(int j=0; j<4; j++)
					object[obj_i].vertex[i][j] = MM[j][0]*tmp_vertex[0] + MM[j][1]*tmp_vertex[1] + MM[j][2]*tmp_vertex[2] + MM[j][3]*tmp_vertex[3];  
			}
		}
		else if(line.substr(0, 8) == "observer"){
			float cx, cy, cz, tilt;
			istringstream strval(line.substr(9));
			strval >> e_posx >> e_posy >> e_posz >> cx >> cy >> cz >> tilt >> zNear >> zFar >> hFov;
			EyeSpace( e_posx, e_posy, e_posz, cx, cy, cz, tilt );	
			ProjSpace( zNear, zFar, hFov );
		}
		else if(line.substr(0, 7) == "display"){
			ScreenSpace();
		}
		else if(line.substr(0, 3) == "end"){
			break;
		}
	}
}

void display() {

    glClearColor( 0.0f, 0.0f, 0.0f, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, ww, 0, wh);

	float r, g, b;
	for(int i=0; i<ww; i++){
		for(int j=0; j<wh; j++){
			vector<float>::iterator color = CBuffer[j][i].begin();
			r = *color;  color++;
			g = *color;  color++;
			b = *color;
			plot_pixel(i, j, r, g, b);
		}
	}
	glFlush();
}

void reshape(int w, int h){
    glViewport(0, 0, w, h);
	 ww = w;
	 wh = h;
}
void release_memory()
{
	for(int j=0; j<num_obj; j++){
		for(int i=0; i<object[j].num_vertices; i++)
			delete[](object[j].vertex[i]);
		delete[](object[j].vertex);

		for(int i=0; i<object[j].num_faces; i++)
			delete[](object[j].face[i]);
		delete[](object[j].face);

		for(int i=0; i<object[j].num_vertices; i++)
			delete[](object[j].VertexColor[i]);
	    delete[](object[j].VertexColor);

		for(int i=0; i<object[j].num_vertices; i++)
			delete[](object[j].VertexNormal[i]);
	    delete[](object[j].VertexNormal);

		delete[](object[j].IsFront);
	}
	delete[](object);

	for(int i=0; i<wh; i++){
		delete[](ZBuffer[i]);
	    delete[](CBuffer[i]);
	}
	delete[](ZBuffer);
	delete[](CBuffer);
}

void main(int argc, char **argv) {
	load_file( "Hw4E.in" );
	glutInit(&argc, argv);
  	glutInitWindowPosition(0, 0);
	glutInitWindowSize(ww, wh);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutCreateWindow("Homework4");
	
	//glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutMainLoop();

    release_memory();
}