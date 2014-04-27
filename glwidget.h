
//These two lines are header guiardians against multiple includes
#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <vector>
#include "CVector.h"
#include <QImage>

using namespace std;

struct material
{
    float reflection;   //refelection coefficient ranging between 0 ~ 1
    float red, green, blue;
    float diffussPower, specularPower;  //used to control the size of bright spots
};

struct sphere
{
    CVector3D pos;
    float size;
    material mat;
};

struct light
{
    CVector3D pos;
    float red, green, blue;
};

struct ray
{
    CVector3D start;
    CVector3D dir;
};

//This is our OpenGL Component we built it on top of QGLWidget
class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    //Constructor for GLWidget
    GLWidget(QWidget *parent = 0);

    //Destructor for GLWidget
    ~GLWidget();

protected:
    //Initialize the OpenGL Graphics Engine
    void initializeGL();

    //All our painting stuff are here
    void paintGL();

    //When user resizes main window, the scrollArea will be resized and it will call this function from
    //its attached GLWidget
    void resizeGL(int width, int height);

private:
    //Makes each face of the polygon
    void makePolygon(int a, int b, int c, int d);
    float getMin(const float &n1, const float &n2);
    float getMax(const float &n1, const float &n2);
    void normalize(CVector3D &v);
    bool hitSphere(const ray &r, const sphere &s, float &t0, float &t1);

    //It is much more efficient if you create all the objects in a scene
    //in one place and group them in a list. Then this will be plugged into graphics hardware at runtime.
    GLuint makeObject();

    //Reference to object. created using makeObject()
    GLuint object;

    //Holds the last mouse position
    QPoint lastPos;
};


#endif
