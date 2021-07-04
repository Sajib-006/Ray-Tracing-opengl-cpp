#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#include <bits/stdc++.h>
#include <GL/glut.h>

#include "bitmap_image.hpp"

#include "1605064.h"
//#include <fstream>
using namespace std;
#include<string.h>
#define pi (2 * acos(0.0))

vector<Object*> objects;
vector<Light*> lights;
void capture();

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
int drawfloor;
double angle;
double angle_rot;
double angle_q,angle_e,angle_a,angle_d,angle_ea,speed;
int k;
bool isShot;
long shotCnt;

int recursion_level, pixels, total_obj, total_lights;



point shotPoints[5000];
point pos = {100, 100, 50};
//point pos = {150, 250, 130};
point u = {0, 0, 1};
point r = {-1 / sqrt(2), 1 / sqrt(2), 0};
point l = {-1 / sqrt(2), -1 / sqrt(2), 0};

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

void drawGrid()
{
	int i;
	if (drawgrid == 1)
	{
		glColor3f(0.6, 0.6, 0.6); //grey
		glBegin(GL_LINES);
		{
			for (i = -8; i <= 8; i++)
			{

				if (i == 0)
					continue; //SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i * 20, -500, 0);
				glVertex3f(i * 20, 500, 0);

				//lines parallel to X-axis
				glVertex3f(-500, i * 20, 0);
				glVertex3f(500, i * 20, 0);
			}
		}
		glEnd();
	}
}
void drawFloor()
{
	//int i,j;
	if (drawfloor == 1)
	{
		//glColor3f(0.6, 0.6, 0.6); //grey
		for (int i = -25; i < 25; i++)
        {
            for(int j = -25; j < 25; j++)
            {
                if((i+j)%2 ==0 ) glColor3f(1,1,1);
                else glColor3f(0,0,0);
                glBegin(GL_QUADS);
                {
                    glVertex3f(i * 20, j * 20, 0);
                    glVertex3f(i * 20 + 20, j * 20, 0);
                    glVertex3f(i * 20 + 20, j * 20 + 20, 0);
                    glVertex3f(i * 20, j * 20 + 20, 0);
                }
                glEnd();
            }

        }


	}
}
void drawSquare(double a)
{
	//glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);
	{
		glVertex3f(a, a, 2);
		glVertex3f(a, -a, 2);
		glVertex3f(-a, -a, 2);
		glVertex3f(-a, a, 2);
	}
	glEnd();
}

void drawCircle(double radius, int segments)
{
	int i;
	struct point points[100];
	glColor3f(0.7, 0.7, 0.7);
	//generate points
	for (i = 0; i <= segments; i++)
	{
		points[i].x = radius * cos(((double)i / (double)segments) * 2 * pi);
		points[i].y = radius * sin(((double)i / (double)segments) * 2 * pi);
	}
	//draw segments using generated points
	for (i = 0; i < segments; i++)
	{
		glBegin(GL_LINES);
		{
			glVertex3f(points[i].x, points[i].y, 0);
			glVertex3f(points[i + 1].x, points[i + 1].y, 0);
		}
		glEnd();
	}
}

void drawCone(double radius, double height, int segments)
{
	int i;
	double shade;
	struct point points[100];
	//generate points
	for (i = 0; i <= segments; i++)
	{
		points[i].x = radius * cos(((double)i / (double)segments) * 2 * pi);
		points[i].y = radius * sin(((double)i / (double)segments) * 2 * pi);
	}
	//draw triangles using generated points
	for (i = 0; i < segments; i++)
	{
		//create shading effect
		if (i < segments / 2)
			shade = 2 * (double)i / (double)segments;
		else
			shade = 2 * (1.0 - (double)i / (double)segments);
		glColor3f(shade, shade, shade);

		glBegin(GL_TRIANGLES);
		{
			glVertex3f(0, 0, height);
			glVertex3f(points[i].x, points[i].y, 0);
			glVertex3f(points[i + 1].x, points[i + 1].y, 0);
		}
		glEnd();
	}
}

void drawSphere(double radius, int slices, int stacks,int part)
{
	struct point points[100][100];
	int i, j;
	double h, r;
	//generate points
	for (i = 0; i <= stacks; i++)
	{
		h = radius * sin(((double)i / (double)stacks) * (pi / 2));
		r = radius * cos(((double)i / (double)stacks) * (pi / 2));
		for (j = 0; j <= slices; j++)
		{
			points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}
	//draw quads using generated points
	for (i = 0; i < stacks; i++)
	{
		//glColor3f((double)i / (double)stacks, (double)i / (double)stacks, (double)i / (double)stacks);
		for (j = 0; j < slices; j++)
		{
		    if(j%2==0)
                glColor3f(1,1,1);
            else glColor3f(0,0,0);
			glBegin(GL_QUADS);
			{
				//upper hemisphere
				if(part==1){
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
				}

				//lower hemisphere
				else if (part==0){
                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
				}

			}
			glEnd();
		}
	}
}

void drawSphereFull(Vector3D center, double radius, int slices, int stacks, double* color)
{
	struct point points[100][100];
	int i, j;
	double h, r;
	glColor3f(color[0],color[1],color[2]);
	
	//generate points
	for (i = 0; i <= stacks; i++)
	{
		h = radius * sin(((double)i / (double)stacks) * (pi / 2));
		r = radius * cos(((double)i / (double)stacks) * (pi / 2));
		for (j = 0; j <= slices; j++)
		{
			points[i][j].x = center.x + r * cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = center.y + r * sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}
	//draw quads using generated points
	for (i = 0; i < stacks; i++)
	{
		//glColor3f((double)i / (double)stacks, (double)i / (double)stacks, (double)i / (double)stacks);
		for (j = 0; j < slices; j++)
		{
			
			glBegin(GL_QUADS);
			{
				//upper hemisphere
				
                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z+center.z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z+center.z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z+center.z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z+center.z);

				//lower hemisphere
				
                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z+center.z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z+center.z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z+center.z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z+center.z);


			}
			glEnd();
		}
	}
}

void drawCylinder(double radius, int slices, int stacks)
{
	struct point points[100][100];
	int i, j;
	double h, r, thetaC;
	//generate points
	for (i = 0; i <= stacks; i++)
	{
	    thetaC = atan(stacks/radius);
		h = radius * tan(((double)i / (double)stacks) * thetaC);
		//r = radius * cos(((double)i / (double)stacks) * (pi / 2));
		for (j = 0; j <= slices; j++)
		{
			points[i][j].x = radius * cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = radius * sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}
	//draw quads using generated points
	for (i = 0; i < stacks; i++)
	{
		//glColor3f((double)i / (double)stacks, (double)i / (double)stacks, (double)i / (double)stacks);
		for (j = 0; j < slices; j++)
		{
		    if(j%2==0)
                glColor3f(1,1,1);
            else glColor3f(0,0,0);
			glBegin(GL_QUADS);
			{
				//upper hemisphere

                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);


				//lower hemisphere

                glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);


			}
			glEnd();
		}
	}
}
void drawCylinderTail(double radius, int slices, int stacks)
{
	struct point points[100][100];
	int i, j;
	double h, r, r1;
	//generate points
	for (i = 0; i <= stacks; i++)
	{
		h = radius * sin(((double)i / (double)stacks) * (pi / 2));
		r = radius * cos(((double)i / (double)stacks) * (pi / 2));
		r1 = 2 * radius - r;
		for (j = 0; j <= slices; j++)
		{
			points[i][j].x = r1 * cos(((double)j / (double)slices) * 2 * pi);
			points[i][j].y = r1 * sin(((double)j / (double)slices) * 2 * pi);
			points[i][j].z = h;
		}
	}
	//draw quads using generated points
	for (i = 0; i < stacks; i++)
	{
		//glColor3f((double)i / (double)stacks, (double)i / (double)stacks, (double)i / (double)stacks);
		for (j = 0; j < slices; j++)
		{
		    if(j%2!=0)
                glColor3f(1,1,1);
            else glColor3f(0,0,0);
			glBegin(GL_QUADS);
			{
				//upper hemisphere

                glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);


			}
			glEnd();
		}
	}
}
void drawSS()
{
	/*glColor3f(1, 0, 0);
	drawSquare(20);

	glRotatef(angle, 0, 0, 1);
	glTranslatef(110, 0, 0);
	glRotatef(2 * angle, 0, 0, 1);
	glColor3f(0, 1, 0);
	drawSquare(15);

	glPushMatrix();
	{
		glRotatef(angle, 0, 0, 1);
		glTranslatef(60, 0, 0);
		glRotatef(2 * angle, 0, 0, 1);
		glColor3f(0, 0, 1);
		drawSquare(10);
	}
	glPopMatrix();

	glRotatef(3 * angle, 0, 0, 1);
	glTranslatef(40, 0, 0);
	glRotatef(4 * angle, 0, 0, 1);
	glColor3f(1, 1, 0);
	drawSquare(5);*/
	//glRotatef(90.0, 0, 0, 1);
	//glTranslatef(40, 0, 0);
	glPushMatrix();
	{
	    glTranslatef(-260, 0, 0);
        glRotatef(-90.0, 0, 1, 0);
        glColor3f(2, 2, 2);
        drawSquare(100);
	}
    glPopMatrix();
    glRotatef(angle_q, 0, 0, 1);
	glPushMatrix();
	{
	    glRotatef(90.0, 0, 1, 0);
        drawSphere(30,48,20,1);
	}
    glPopMatrix();
    glRotatef(angle_e, 0, 1, 0);
    glPushMatrix();
	{
	    glRotatef(90.0, 0, 1, 0);
        drawSphere(30,48,20,0);
	}
    glPopMatrix();
    glRotatef(angle_a, 0, 1, 0);
    glRotatef(angle_d, 1, 0, 0);
	glPushMatrix();
	{
	    glTranslatef(-45, 0, 0);
	    glRotatef(90.0, 0, 1, 0);
	    drawSphere(15,48,20,1);
	}
    glPopMatrix();
    glPushMatrix();
	{
	    glTranslatef(-95, 0, 0);
	    glRotatef(90.0, 0, 1, 0);
        drawCylinder(15,48,50);
	}
    glPopMatrix();
    glPushMatrix();
	{
	    glTranslatef(-145, 0, 0);
        glRotatef(-90.0, 0, 1, 0);
        drawCylinderTail(15,48,20);
	}
    glPopMatrix();
	/*glBegin(GL_POINTS); //starts drawing of points
      glVertex3f(20.0f,10.0f,2.0f);//upper-right corner
      //glVertex3f(-10.0f,-10.0f,0.0f);//lower-left corner
    glEnd();//end drawing of points*/
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
void rotate_gun_Z_anti_Clockwise()
{
    if(angle_q < 65.0) angle_q += speed ;
}
void rotate_gun_Z_Clockwise()
{
    if(angle_q > -65.0) angle_q -= speed ;
}
void rotate_gun_Y_anti_Clockwise()
{
    if(angle_e < 65.0) angle_e += speed ;
    angle_ea += speed;
}
void rotate_gun_Y_Clockwise()
{
    if(angle_e > -65.0) angle_e -= speed ;
    angle_ea -= speed;
}
void rotate_cylinder_Y_anti_Clockwise()
{
    if(angle_a < 65.0) angle_a += speed ;
    angle_ea += speed;
}
void rotate_cylinder_Y_Clockwise()
{
    if(angle_a > -65.0) angle_a -= speed ;
    angle_ea -= speed;
}
void rotate_cylinder_by_axis_anti_Clockwise()
{
    if(angle_d < 105.0) angle_d += speed ;
}
void rotate_cylinder_by_axis_Clockwise()
{
    if(angle_d > -105.0) angle_d -= speed ;
}
void keyboardListener(unsigned char key, int x, int y)
{
	switch (key)
	{
	case '0':
		capture();
		break;
	case '1':
		drawgrid=1-drawgrid;
		drawfloor=1-drawfloor;
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
    case 'q':
        rotate_gun_Z_anti_Clockwise();
		break;
    case 'w':
        rotate_gun_Z_Clockwise();
		break;
    case 'e':
        rotate_gun_Y_anti_Clockwise();
        //angle_ea = angle_e;
		break;
    case 'r':
        rotate_gun_Y_Clockwise();
        //angle_ea = angle_e;
		break;
    case 'a':
        rotate_cylinder_Y_anti_Clockwise();
        //angle_ea = angle_a;
		break;
    case 's':
        rotate_cylinder_Y_Clockwise();
        //angle_ea = angle_a;
		break;
    case 'd':
        rotate_cylinder_by_axis_anti_Clockwise();
		break;
    case 'f':
        rotate_cylinder_by_axis_Clockwise();
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
void shoot()
{


    for(int i=1; i<=shotCnt; i++)
    {
        glPushMatrix();
        {
            //printf("%d %d %d\n",shotPoints[i].x, shotPoints[i].y, shotPoints[i].z);
            glTranslatef(shotPoints[i].x, shotPoints[i].y, shotPoints[i].z);
            glRotatef(-90.0, 0, 1, 0);
            glColor3f(1, 0, 0);
            drawSquare(5);
        }
        glPopMatrix();

    }


    //printf("shoot %d\n",shotCnt);

}
void mouseListener(int button, int state, int x, int y)
{ //x, y is the x-y of the screen (2D)
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
	    if (state == GLUT_DOWN)
		{

            point vertex = {-160*cos(angle_q*pi/180.0),-160*sin(angle_q*pi/180.0),160*sin(angle_ea*pi/180.0)};
            //printf("%d %d %d\n",vertex.x, vertex.y,vertex.z);
            if(vertex.y > -100.0 && vertex.y < 100.0 && vertex.z > -100.0 && vertex.z < 100.0){
                shotCnt++;
                shotPoints[shotCnt] = {vertex.x-100.0, vertex.y,vertex.z};
            }

		}
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
	drawGrid();
    drawFloor();
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
		
		// if((s = dynamic_cast<Sphere*>(*iter)) != nullptr){
		// 	drawSphereFull(s->reference_point,s->length,30,24,s->color);
		// }
		(*iter)->draw();
			
	}
	for(int i=0; i<lights.size(); i++){
		lights[i]->draw();
	}
	// delete s;


	//drawCone(20, 50, 24);

	

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
	

    // for(int i=400;i<450;i++){
    //     for(int j=50;j<150;j++){
    //         image.set_pixel(i,j,255,0,0);
    //     }
    //     for(int j=200;j<300;j++){
    //         image.set_pixel(i,j,0,255,255);
    //     }
    // }
	
	double planeDistance, windowWidth, windowHeight, viewAngle, du, dv, imageWidth, imageHeight;
	point eye, topleft, curPixel;
	windowWidth = 500;
	windowHeight = 500;
	viewAngle = 80.0 * pi / 180.0;
	imageHeight = pixels;
	imageWidth = pixels;
	eye = pos;
	printPoint(eye);
	cout<< "angles---"<< viewAngle << tan(viewAngle/2.0) << tan(40)<<endl;
	planeDistance = (windowHeight/2.0) / tan(viewAngle/2.0);
	topleft = eye + l*planeDistance - r*(windowWidth/2) + u*(windowHeight/2);
	du = windowWidth/imageWidth;
	dv = windowHeight/imageHeight;
	cout<<"------"<<endl;
	printPoint(l*planeDistance);
	printPoint(r*(windowWidth/2));
	printPoint(u*(windowHeight/2));

	printPoint(l);
	printPoint(r);
	printPoint(u);

	
	cout<<"------"<<endl;
	Vector3D p(1,1,3),q(2,1,2),s(3,2,2);
	//s = p - q;
	//s.print();
	cout<<determinant(p,q,s)<<endl;
	// Choose middle of the grid cell
	printPoint(topleft);
	topleft = topleft + r*(0.5*du) - u*(0.5*dv);
	cout<<"topleft ";
	printPoint(topleft);
	int nearest;
	double t, tMin=INT_MAX;
	int cnt=0;
	
	cout<< "plane dis: "<<planeDistance<<" "<<" "<<du<<" "<<dv<<" "<<endl;

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
	ofstream fout("output.txt");
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
			if(i==0 && j==0) ray->print(); //----------------
			for(int i=0; i<objects.size(); i++){
				//cout<<"level 0: ";
				t = objects[i]->intersect(ray, dummyColor, 0);
				
				if(t > 0 && t < tMin){
					tMin = t;
					nearest_obj = objects[i];
					changed = true;
					test_cnt++;
					fout<<cnt<<" "<<t<<endl; //--------------
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
	
	fout.close();
    image.save_image("out.bmp");;
	cout << "image done" << endl;
	cout<<test_cnt<<endl;
}

void init()
{
	//codes for initialization
	drawgrid = 0;
	drawaxes = 0;
	drawfloor = 0;
	cameraHeight = 150.0;
	cameraAngle = 1.0;
	angle = 0;
    angle_rot = pi/180.0* 2;
    k = 2;
    angle_q = 0.0;
    angle_e = 0.0;
    angle_a = 0.0;
    angle_d = 0.0;
    angle_ea = 0.0;
    speed = 5.0;
    isShot = false;
    shotCnt = 0;
	//clear the screen
	glClearColor(0, 0, 0, 0);

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
    /*ifstream fin;
    fin.open("scene.txt");
    fin >> recursion_level;
    fin >> pixels;
    fin>> total_obj;*/
    FILE *scene;

	Object *temp = new Floor(1000, 20);
	objects.push_back(temp);
    printf("Value of");
    if ((scene = fopen("scene.txt","r")) == NULL){
       printf("Error! opening file");

       // Program exits if the file pointer returns NULL.
       exit(1);
    }
    fscanf(scene,"%d", &recursion_level);
    fscanf(scene,"%d", &pixels);
    fscanf(scene,"%d", &total_obj);
    printf("%d %d %d\n\n", recursion_level,pixels,total_obj);
    char command[10];
	for(int i=0; i<total_obj; i++){

		fscanf(scene, "%s", command);
        printf("%s\n",command);
        if(!strcmp(command, "sphere")){
            printf("Sphere found\n");
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
            sphere->print();
			objects.push_back(sphere);
			
        }
		else if(!strcmp(command, "triangle")){
            printf("Triangle found\n");
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
            temp->print();
			objects.push_back(temp);
			
        }
		else if(!strcmp(command, "general")){
            printf("general found\n");
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
            temp->print();
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
		temp->print();
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
    printf("hello\n");
    loadData();
	printf("hello2\n");
	glEnable(GL_DEPTH_TEST); //enable Depth Testing

	glutDisplayFunc(display); //display callback function
	glutIdleFunc(animate);	  //what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop(); //The main loop of OpenGL

	return 0;
}
