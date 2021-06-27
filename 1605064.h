//#include <iostream>
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
};

class Ray{
    Vector3D start;
    Vector3D dir; // normalized
public:
    Ray(Vector3D start, Vector3D dir){
        this->start = start;
        this->dir = dir;
    }
};


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
};

class Sphere : public Object{
public:
    Sphere(Vector3D center, double radius){
        reference_point = center;
        length = radius;
    }
    void draw(){

    }
    double intersect(Ray *r, double *color, int level){
        return -1.0;
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
class Triangle: public Object{
    Vector3D vertex[3];
public:
    Triangle(Vector3D v1, Vector3D v2, Vector3D v3){
        this->vertex[0] = v1;
        this->vertex[1] = v2;
        this->vertex[2] = v3;
    }
    void draw(){

    }
    double intersect(Ray *r, double *color, int level){
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

    }
    double intersect(Ray *r, double *color, int level){
        return -1.0;
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


