//#include <iostream>
#include<bits/stdc++.h>
#define pi (2 * acos(0.0))
using namespace std;


class Vector3D{

public:
    double x,y,z;
    Vector3D(){

    }
    Vector3D(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;

    }
    Vector3D operator + (const Vector3D& vec) {
        Vector3D temp;
        temp.x = this->x + vec.x;
        temp.y = this->y + vec.y;
        temp.z = this->z + vec.z;
        return temp;
    }
    Vector3D operator - (const Vector3D& vec) {
        Vector3D temp;
        temp.x = this->x - vec.x;
        temp.y = this->y - vec.y;
        temp.z = this->z - vec.z;
        return temp;
    }

    double dotProduct(Vector3D b){
        double val;
        val = this->x * b.x + this->y * b.y + this->z * b.z;
        return val;
    }
    Vector3D crossProduct(Vector3D b){
        Vector3D res;
        res.x = this->y * b.z - b.y * this->z;
        res.y = this->z * b.x - b.z * this->x;
        res.z = this->x * b.y - b.x * this->y;
        return res;

    }
    Vector3D normalize(){
        double val;
        val = sqrt(x*x + y*y + z*z);
        this->x/=val;
        this->y/=val;
        this->z/=val;
    }
    void print(){
        cout<< "Vector3D: x:"<<x<<" y:"<<y<<" z:"<<z<<endl;
    }
};




class point{
public:

	double x, y, z;
	point(double px=0, double py=0, double pz=0){
        x = px;
        y = py;
        z = pz;
    }
    
	point operator+(const point& p) const
    {
        return point(p.x+x, p.y+y , p.z+z);
    }
	point operator-(const point& p) const
    {
        return point(x-p.x, y-p.y , z-p.z);
    }
	point operator*(const point& p)
    {
        return point((p.x)*x, (p.y)*y , (p.z)*z);
        // point temp;
        // temp.x = x * p.x;
        // temp.y = y * p.y;
        // temp.z = z * p.z;
        // return temp;
    }
    point operator+(const double n)
    {
        return point(n+x, n+y , n+z);
    }
    point operator-(const double n)
    {
        return point(x-n, y-n , z-n);
    }
    point operator*(const double n)
    {
        return point(n*x, n*y , n*z);
    }
    
};

class Ray{
public:
    Vector3D start;
    Vector3D dir; // normalized

    Ray(Vector3D start, Vector3D dir){
        this->start = start;
        this->dir = dir;
        this->dir.normalize();
    }
    Ray(point start, point dir){
        this->start.x = start.x;
        this->start.y = start.y;
        this->start.z = start.z;
        this->dir.x = dir.x;
        this->dir.y = dir.y;
        this->dir.z = dir.z;
        this->dir.normalize();
    }
    void print(){
        cout<<"Ray printing.."<<endl;
        start.print();
        dir.print();
    }
};

Vector3D get_intersect_point(Ray *r,double t){
    Vector3D res;
    res.x = r->start.x + r->dir.x * t;
    res.y = r->start.y + r->dir.y * t;
    res.z = r->start.z + r->dir.z * t;
    return res;
}

class Object{
public:
    Vector3D reference_point;
    double height, width, length;
    double color[3];
    double coEfficients[4];
    int shine;


    Object(){
         //cout("Object claas constructor\n");

    }
    virtual void draw(){
        //cout("Object claas draw()\n");
    }
    void setColor(double r, double g, double b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
    void setShine(int shine){
        this->shine = shine;
    }
    void setCoEfficients(double ka, double kd, double ks, double kr){
        this->coEfficients[0] = ka;
        this->coEfficients[1] = kd;
        this->coEfficients[2] = ks;
        this->coEfficients[3] = kr;
    }
    virtual double intersect(Ray *r, double *color, int level){
        return -1.0;
    }
    virtual void print(){}

    double* getAmbientColor(){
        double *ambient_color = new double[3];
        ambient_color[0] = this->color[0] * this->coEfficients[0];
        ambient_color[1] = this->color[1] * this->coEfficients[0];
        ambient_color[2] = this->color[2] * this->coEfficients[0];
        return ambient_color;
    }
};

class Sphere : public Object{
public:
    int slices, stacks;
    
    Sphere(Vector3D center, double radius){
        reference_point = center;
        length = radius;
        slices = 30;
        stacks = 24;
    }
    void setdrawparameter(int slices, int stacks){
        this->slices = slices;
        this->stacks = stacks;
    }
    void draw(){
        
        struct point points[100][100];
        int i, j;
        double h, r;
        glColor3f(color[0],color[1],color[2]);
        
        //generate points
        for (i = 0; i <= stacks; i++)
        {
            h = length * sin(((double)i / (double)stacks) * (pi / 2));
            r = length * cos(((double)i / (double)stacks) * (pi / 2));
            for (j = 0; j <= slices; j++)
            {
                points[i][j].x = reference_point.x + r * cos(((double)j / (double)slices) * 2 * pi);
                points[i][j].y = reference_point.y + r * sin(((double)j / (double)slices) * 2 * pi);
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
                    
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z+reference_point.z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z+reference_point.z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z+reference_point.z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z+reference_point.z);

                    //lower hemisphere
                    
                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z+reference_point.z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z+reference_point.z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z+reference_point.z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z+reference_point.z);


                }
                glEnd();
            }
        }
    }
    // double intersect(Ray *r, double *color, int level){
    //     Vector3D R0;
    //     double val, tp, d_sqr, t1, t;
    //     cout<< "In intersect(): ";
    //     r->dir.print();

    //     R0 = r->start - reference_point;

    //     tp = R0.dotProduct(r->dir);                 //tp = R0.Rd
    //     val = R0.dotProduct(R0) - length * length;  //val = R0.R0-r^2
    //     cout<<"val "<<val<<" tp:"<<tp<<endl;
    //     if(val > 0.0f && tp > 0.0f) return -1.0;    //here tp is negative value, so if tp>0 then creates obtuse angle

    //     d_sqr = R0.dotProduct(R0) - tp * tp;
        
    //     t1 = sqrt(length*length -d_sqr);
    //     cout<<"d^2: "<<d_sqr<<"t1: "<<t1<<endl;
    //     if(val > 0.0f) t = tp - t1;                 // if eye outside
    //     else if(val < 0.0f) t = tp + t1;            // if eye inside
    //     else t = tp;

    //     color = getAmbientColor();
    //     return t;
    // }

    double intersect(Ray *r, double *color_out, int level){
        Vector3D R0;
        double val, b, d_sqr, t1, t, discrim;
        //cout<< "In intersect(): ";
        //r->dir.print();

        R0 = r->start - reference_point;
        
        b = 2 * R0.dotProduct(r->dir);                 //tp = R0.Rd
        val = R0.dotProduct(R0) - length * length;  //val = R0.R0-r^2
        
        discrim = b * b - 4 * val;
        //cout<<"val "<<val<<" b:"<<b<<"discrim: "<<discrim<<endl;
        //R0.print();
        //r->dir.print();
        //cout<<"discrim: "<<discrim<<endl;
        
        if(discrim < 0.0f ) return -1.0;    //here tp is negative value, so if tp>0 then creates obtuse angle
        discrim = sqrt(discrim);
        t = min((-b+discrim)/2.0, (-b-discrim)/2.0 );
        //color_out = getAmbientColor();
        color_out[0] = color[0];
        color_out[1] = color[1];
        color_out[2] = color[2];

        
        return t;
    }
    void print(){
        printf("-----Sphere------\n");
        printf("Center: %lf %lf %lf\n",reference_point.x,reference_point.y,reference_point.z);
        printf("Radius: %lf\n",length);
        printf("Color: %lf %lf %lf\n",color[0],color[1],color[2]);
        printf("Coefficients: %lf %lf %lf %lf\n",coEfficients[0],coEfficients[1],coEfficients[2],coEfficients[3]);
        printf("Shine: %d\n\n",shine);
    }

};
double determinant(Vector3D a, Vector3D b, Vector3D c){
    double det;
    det = a.x*(b.y*c.z - b.z*c.y) - b.x*(a.y*c.z - a.z*c.y) + c.x*(a.y*b.z - a.z*b.y);
    return det;
}

class Triangle: public Object{
    Vector3D vertex[3];
public:
    Triangle(Vector3D v1, Vector3D v2, Vector3D v3){
        this->vertex[0] = v1;
        this->vertex[1] = v2;
        this->vertex[2] = v3;
    }
    void draw(){
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3d(vertex[0].x, vertex[0].y, vertex[0].z);
            glVertex3d(vertex[1].x, vertex[1].y, vertex[1].z);
            glVertex3d(vertex[2].x, vertex[2].y, vertex[2].z);
        }
        glEnd();
    }
    double intersect(Ray *r, double *color_out, int level){
        double b, y, t, det;
        det = determinant(vertex[0]-vertex[1], vertex[0]-vertex[2], r->dir);
        b = determinant(vertex[0]-r->start, vertex[0]-vertex[2], r->dir)/det;
        y = determinant(vertex[0]-vertex[1], vertex[0]-r->start, r->dir)/det;
        t = determinant(vertex[0]-vertex[1], vertex[0]-vertex[2], vertex[0]-r->start)/det;

        
        color_out[0] = color[0];
        color_out[1] = color[1];
        color_out[2] = color[2];
        if(b>0 && y>0 && b+y<1 && t>0) {
            //cout<<"----------------"<<b<<" "<<y<<" "<<t<<endl;
            return t;
        }
        return -1.0;
    }
};

class Floor: public Object{
public:
    Floor(double floorWidth, double tileWidth){
        reference_point = Vector3D(-floorWidth/2,-floorWidth/2,0);
        length = tileWidth;
    }
    void draw(){
        //cout<<"in floor draw()"<<endl;
        int row = -(reference_point.x / length);
        int col = -(reference_point.y / length);
        //cout<<row<<" "<<col<<endl;
        for (int i = -row; i < row; i++)
        {
            for(int j = -col; j < col; j++)
            {
                if((i+j)%2 ==0 ) glColor3f(1,1,1);
                else glColor3f(0,0,0);
                glBegin(GL_QUADS);
                {
                    glVertex3f(i * length, j * length, 0);
                    glVertex3f(i * length + length, j * length, 0);
                    glVertex3f(i * length + length, j * length + length, 0);
                    glVertex3f(i * length, j * length + length, 0);
                }
                glEnd();
            }

        }
    }
    double intersect(Ray *r, double *color_out, int level){
        Vector3D n(0,0,1);
        double temp, temp2, t, D;
        D = 0;
        temp = (D + r->start.dotProduct(n));
        temp2 = r->dir.dotProduct(n);
        t = -temp/temp2;
        Vector3D int_point = get_intersect_point(r,t);
        if(int_point.x >= reference_point.x && int_point.x <= -reference_point.x && int_point.y >= reference_point.y && int_point.y <= -reference_point.y){
            int val = int(ceil(int_point.x/length)) + int(ceil(int_point.y/length));
            if( val%2 == 0 ) {
                color_out[0] = 1;
                color_out[1] = 1;
                color_out[2] = 1;
            }
            else {
                color_out[0] = 0;
                color_out[1] = 0;
                color_out[2] = 0;
            } 
        }
        else{
            color_out[0] = 0;
            color_out[1] = 0;
            color_out[2] = 0;
        }
        return t;
        //return -1.0;
    }
};
class Light{
    Vector3D light_pos;
    double color[3];
public:
    Light(Vector3D light_pos){
        this->light_pos = light_pos;
    }
    Light(Vector3D light_pos, double r, double g, double b){
        this->light_pos = light_pos;
        this->color[0] = r;
        this->color[1] = g;
        this->color[2] = b;
    }
    void setColor(double r, double g, double b){
        this->color[0] = r;
        this->color[1] = g;
        this->color[2] = b;
    }

};


