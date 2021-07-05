#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#include <bits/stdc++.h>
#include <GL/glut.h>

#include "bitmap_image.hpp"

#include "1605064.h"
using namespace std;
#include<string.h>
#define pi (2 * acos(0.0))

vector<Object*> objects;
vector<Light*> lights;
void capture();

int drawaxes;
double angle;
double angle_rot;
point pos, u, r, l;
int k;

int recursion_level, pixels, total_obj, total_lights;
vector<Object*>::iterator iterO, endO;
vector<Light*>::iterator iterL, endL;

void printPoint(point p){
        cout<< "x:"<< p.x << " " << "y:"<< p.y << " " << "z:"<< p.z <<endl;
}

void drawAxes()
{
	if (drawaxes == 1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);
		{
			glVertex3f(100, 0, 0);
			glVertex3f(-100, 0, 0);

			glVertex3f(0, -100, 0);
			glVertex3f(0, 100, 0);

			glVertex3f(0, 0, 100);
			glVertex3f(0, 0, -100);
		}
		glEnd();
	}
}


struct point Cross_multiply(struct point a, struct point b)
{
    struct point c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}
struct point Change_Point(struct point a, struct point b,bool flag)
{
    struct point c;
    if (flag == true){
        c.x = a.x * cos(angle_rot) + b.x * sin(angle_rot);
        c.y = a.y * cos(angle_rot) + b.y * sin(angle_rot);
        c.z = a.z * cos(angle_rot) + b.z * sin(angle_rot);
    }
    else{
        c.x = a.x * cos(angle_rot) - b.x * sin(angle_rot);
        c.y = a.y * cos(angle_rot) - b.y * sin(angle_rot);
        c.z = a.z * cos(angle_rot) - b.z * sin(angle_rot);
    }
    return c;
}
void look_left()
{
    //rotation along u
    point c = Cross_multiply(u,r);
    r = Change_Point(r,c,false);
    point d = Cross_multiply(u,l);
    l = Change_Point(l,d,false);
}
void look_right()
{
    //rotation along u
    point c = Cross_multiply(u,r);
    r = Change_Point(r,c,true);
    point d = Cross_multiply(u,l);
    l = Change_Point(l,d,true);
}
void look_up()
{
    //rotation along r
    point c = Cross_multiply(r,l);
    l = Change_Point(l,c,false);
    point d = Cross_multiply(r,u);
    u = Change_Point(u,d,false);
}
void look_down()
{
    //rotation along r
    point c = Cross_multiply(r,l);
    l = Change_Point(l,c,true);
    point d = Cross_multiply(r,u);
    u = Change_Point(u,d,true);
}
void tilt_clockwise()
{
    //rotation along l
    point c = Cross_multiply(l,u);
    u = Change_Point(u,c,false);
    point d = Cross_multiply(l,r);
    r = Change_Point(r,d,false);
}
void tilt_anti_Clockwise()
{
    //rotation along l
    point c = Cross_multiply(l,u);
    u = Change_Point(u,c,true);
    point d = Cross_multiply(l,r);
    r = Change_Point(r,d,true);
}

void keyboardListener(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '0':
		capture();
		break;
	case '1':
        look_left();
		break;
    case '2':
        look_right();
		break;
    case '3':
        look_up();
		break;
    case '4':
        look_down();
		break;
    case '5':
        tilt_clockwise();
		break;
    case '6':
        tilt_anti_Clockwise();
		break;
	case '7':
		//Free memory
		for(iterO = objects.begin(), endO = objects.end() ; iterO != endO; ++iterO) delete *iterO;
		objects.clear();	
		for(iterL = lights.begin(), endL = lights.end() ; iterL != endL; ++iterL) delete *iterL;	
		lights.clear();
		cout<<"Memory freed"<<endl;
		exit(0);
		break;
	default:
		break;
	}
}

void specialKeyListener(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_DOWN: //down arrow key
		//cameraHeight -= 3.0;
		pos.x -= k*l.x;
		pos.y -= k*l.y;
		pos.z -= k*l.z;
		break;
	case GLUT_KEY_UP: // up arrow key
		//cameraHeight += 3.0;
		pos.x += k*l.x;
		pos.y += k*l.y;
		pos.z += k*l.z;
		break;

	case GLUT_KEY_RIGHT:
		//cameraAngle += 0.03;
		pos.x -= k*r.x;
		pos.y -= k*r.y;
		pos.z -= k*r.z;
		break;
	case GLUT_KEY_LEFT:
		//cameraAngle -= 0.03;
		pos.x += k*r.x;
		pos.y += k*r.y;
		pos.z += k*r.z;
		break;

	case GLUT_KEY_PAGE_UP:
	    pos.x += k*u.x;
		pos.y += k*u.y;
		pos.z += k*u.z;
		break;
	case GLUT_KEY_PAGE_DOWN:
	    pos.x -= k*u.x;
		pos.y -= k*u.y;
		pos.z -= k*u.z;
		break;

	case GLUT_KEY_INSERT:
		break;

	case GLUT_KEY_HOME:
		break;
	case GLUT_KEY_END:
		break;

	default:
		break;
	}
}

void mouseListener(int button, int state, int x, int y)
{ //x, y is the x-y of the screen (2D)
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
	    exit(0);
		break;

	case GLUT_RIGHT_BUTTON:
		//........
		if (state == GLUT_DOWN)
		{ // 2 times?? in ONE click? -- solution is checking DOWN or UP
			drawaxes = 1 - drawaxes;
		}

		break;

	case GLUT_MIDDLE_BUTTON:
		//........
		break;

	default:
		break;
	}
}

void display()
{

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0, 0, 0, 0); //color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x, pos.y, pos.z,  pos.x + l.x, pos.y + l.y, pos.z + l.z,   u.x, u.y, u.z);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);

	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	//drawGrid();
    //drawFloor();
	//glColor3f(1,0,0);
	//drawSquare(10);
    //glTranslatef(-95, 0, 0);
    //if(isShot) shoot();
    //drawSquare(10);
    //shoot();
	//drawSS();

	// drawCircle(30, 24);
	// Vector3D center(40,0,10);
	// double color[3] = {0,1,0};
	// double r = 20.0;
	// drawSphereFull(center,r,50,40,color);

	vector<Object*>::iterator iter, end;
	
	for(iter = objects.begin(), end = objects.end() ; iter != end; ++iter) {
		(*iter)->draw();	
	}
	for(int i=0; i<lights.size(); i++){
		lights[i]->draw();
	}
	// delete s;

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}

void animate()
{
	angle += 0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

double clipColor(double color){
	if (color > 1.0f) color = 1.0f;
	if (color < 0.0f) color = 0.0f;
	return color;
}

void capture()
{
	cout<<"Capturing Image..."<<endl;
	double planeDistance, windowWidth, windowHeight, viewAngle, du, dv, imageWidth, imageHeight;
	point eye, topleft, curPixel;
	windowWidth = 500;
	windowHeight = 500;
	viewAngle = 80.0 * pi / 180.0;
	imageHeight = pixels;
	imageWidth = pixels;
	eye = pos;
	cout<< "Camera Position: ";
	printPoint(eye);
	
	planeDistance = (windowHeight/2.0) / tan(viewAngle/2.0);
	topleft = eye + l*planeDistance - r*(windowWidth/2) + u*(windowHeight/2);
	du = windowWidth/imageWidth;
	dv = windowHeight/imageHeight;

	// Choose middle of the grid cell
	topleft = topleft + r*(0.5*du) - u*(0.5*dv);

	int nearest;
	double t, tMin=INT_MAX;
	int cnt=0;


	//test code started
	// point r_dir = {-0.528493, 0.848938, 0};
	// double *dummyColor = new double[3];
	// Ray *ray = new Ray(eye, r_dir);
	// for(int i=0; i<objects.size(); i++){
	// 	t = objects[i]->intersect(ray, dummyColor, 0);
	// 	cout<<i<<" t: "<<t<<endl;
	// }

	// point st = {0, 100, 10};
	// point r_dir = {0, 1, 0};
	// double *dummyColor = new double[3];
	// Ray *ray = new Ray(st, r_dir);
	// for(int i=0; i<objects.size(); i++){
	// 	t = objects[i]->intersect(ray, dummyColor, 0);
	// 	cout<<i<<" t: "<<t<<endl;
	// }
	//test code ended
	//ofstream fout("output.txt");
	bitmap_image image(imageWidth,imageHeight);
	// set background color
    for(int i=0;i<imageWidth;i++){
        for(int j=0;j<imageHeight;j++){
            image.set_pixel(i,j,0,0,0);
        }
        
    }
	int test_cnt=0;
	for(int i=0; i< imageWidth; i++){
		for(int j=0; j<imageHeight; j++){
			cnt++;
			tMin = INT_MAX;
			curPixel = topleft + r*(i*du) - u*(j*dv);
			//cout<<cnt++<<" ";
			//printPoint(curPixel);
			Ray *ray = new Ray(eye, curPixel-eye);
			double *color = new double[3];
			double *dummyColor = new double[3];
			Object *nearest_obj = 0;
			bool changed = false;
			//if(i==0 && j==0) ray->print(); //----------------
			for(int i=0; i<objects.size(); i++){
				//cout<<"level 0: ";
				t = objects[i]->intersect(ray, dummyColor, 0);
				
				if(t > 0 && t < tMin){
					tMin = t;
					nearest_obj = objects[i];
					changed = true;
					test_cnt++;
					//fout<<cnt<<" "<<t<<endl; //--------------
				}
			}
			if(changed){
				//cout<<"level 1: ";
				tMin = nearest_obj->intersect(ray, color, 1);
				double color_r = clipColor(color[0]) * 255;
				double color_g = clipColor(color[1]) * 255;
				double color_b = clipColor(color[2]) * 255;
				//cout<< color_r << " " << color_g << " " << color_b << endl; //------------
				image.set_pixel(i,j,color_r,color_g,color_b);
				//delete nearest_obj;
			}
			
			
			
			delete ray;
			delete color;
			delete dummyColor;
			//delete nearest_obj;
			
		}
	}
	
	//fout.close();
    image.save_image("out.bmp");;
	cout << "Captured!" << endl;
	//cout<<test_cnt<<endl;
}

void init()
{
	//codes for initialization
	drawaxes = 0;
	//clear the screen
	glClearColor(0, 0, 0, 0);
	angle = 0;
    angle_rot = pi/180.0* 2;
	k = 2;
	pos = {100, 100, 50};
	//point pos = {150, 250, 130};
	u = {0, 0, 1};
	r = {-1 / sqrt(2), 1 / sqrt(2), 0};
	l = {-1 / sqrt(2), -1 / sqrt(2), 0};

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80, 1, 1, 1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}
void loadData()
{
    FILE *scene;

	Object *temp = new Floor(1000, 20);
	objects.push_back(temp);
    //printf("Value of\n");
    if ((scene = fopen("scene.txt","r")) == NULL){
       printf("Error! opening file");
       exit(1);
    }
    fscanf(scene,"%d", &recursion_level);
    fscanf(scene,"%d", &pixels);
    fscanf(scene,"%d", &total_obj);
    //printf("%d %d %d\n\n", recursion_level,pixels,total_obj);
    char command[10];
	for(int i=0; i<total_obj; i++){

		fscanf(scene, "%s", command);
        //printf("%s\n",command);
        if(!strcmp(command, "sphere")){
            //printf("Sphere found\n");
            double x,y,z,r,amb,diff,spec,coeff;
			double color[3];
            int shine;
            fscanf(scene,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",&x,&y,&z,&r,&color[0],&color[1],&color[3],&amb,&diff,&spec,&coeff,&shine);
            Vector3D center(x,y,z);
			Object *sphere;
            sphere = new Sphere(center,r);
            sphere->setColor(color[0],color[1],color[3]);
            sphere->setCoEfficients(amb,diff,spec,coeff);
            sphere->setShine(shine);
            //sphere->print();
			objects.push_back(sphere);
			
        }
		else if(!strcmp(command, "triangle")){
            //printf("Triangle found\n");
			Vector3D v1,v2,v3;
            double amb,diff,spec,coeff;
			double color[3];
            int shine;
			fscanf(scene,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&v1.x,&v1.y,&v1.z,&v2.x,&v2.y,&v2.z,&v3.x,&v3.y,&v3.z);
            fscanf(scene,"%lf %lf %lf %lf %lf %lf %lf %d",&color[0],&color[1],&color[3],&amb,&diff,&spec,&coeff,&shine);
            
			Object *temp;
            temp = new Triangle(v1,v2,v3);
            temp->setColor(color[0],color[1],color[3]);
            temp->setCoEfficients(amb,diff,spec,coeff);
            temp->setShine(shine);
            //temp->print();
			objects.push_back(temp);
			
        }
		else if(!strcmp(command, "general")){
            //printf("general found\n");
			double a,b,c,d,e,f,g,h,i,j;
			Vector3D center;
			double length,width,height;
            double amb,diff,spec,coeff;
			double color[3];
            int shine;
			fscanf(scene,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&a,&b,&c,&d,&e,&f,&g,&h,&i,&j);
			fscanf(scene,"%lf %lf %lf %lf %lf %lf",&center.x,&center.y,&center.z,&length,&width,&height);
            fscanf(scene,"%lf %lf %lf %lf %lf %lf %lf %d",&color[0],&color[1],&color[3],&amb,&diff,&spec,&coeff,&shine);
            
			Object *temp;
            temp = new Quadratic(a,b,c,d,e,f,g,h,i,j,center,length,width,height);
            temp->setColor(color[0],color[1],color[3]);
            temp->setCoEfficients(amb,diff,spec,coeff);
            temp->setShine(shine);
            //temp->print();
			objects.push_back(temp);
			
        }

    }
	
	fscanf(scene,"%d", &total_lights);
	for(int i=0; i<total_lights; i++){
		Light *temp;
		Vector3D v;
		double color[3];
		fscanf(scene,"%lf %lf %lf",&v.x,&v.y,&v.z);
		fscanf(scene,"%lf %lf %lf",&color[0],&color[1],&color[2]);
		temp = new Light(v,color[0],color[1],color[2]);
		//temp->print();
		lights.push_back(temp);
	}
}



int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); //Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();
    printf("loadin data..\n");
    loadData();
	printf("data loaded..\n");
	glEnable(GL_DEPTH_TEST); //enable Depth Testing

	glutDisplayFunc(display); //display callback function
	glutIdleFunc(animate);	  //what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);
	//cout<<"2"<<endl;
	glutMainLoop(); //The main loop of OpenGL

	//Free memory
	vector<Object*>::iterator iter, end;
	for(iter = objects.begin(), end = objects.end() ; iter != end; ++iter) delete *iter;
	objects.clear();	
	vector<Light*>::iterator iterL, endL;
	for(iterL = lights.begin(), endL = lights.end() ; iterL != endL; ++iterL) delete *iterL;	
	lights.clear();
	//cout<<"End"<<endl;
	return 0;
}
