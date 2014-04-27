#include <QtGui>
#include <QtOpenGL>
#include <math.h>
#include "glwidget.h"

//ray tracer sampleing size (pictureSize by pictureSize)
#define pictureSize 512

//ray tracer camera position
#define eyePositionX pictureSize/2
#define eyePositionY pictureSize/2
#define eyePositionZ -300

#define Left	-2.0
#define Right	+2.0
#define Bottom	-2.0
#define Top	+2.0
#define Near    -10.0
#define Far     +10.0

GLfloat vertices [][3] = {{-1.0, -1.0, 1.0}, {-1.0, 1.0, 1.0}, {1.0, 1.0, 1.0},
                          {1.0, -1.0, 1.0}, {-1.0, -1.0, -1.0}, {-1.0, 1.0, -1.0},
                          {1.0, 1.0, -1.0}, {1.0, -1.0, -1.0}};

GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(parent)
{
    //initialize the object to zero
    object = 0;
}

GLWidget::~GLWidget()
{
    makeCurrent();

    //Delete newly created object
    glDeleteLists(object, 1);
}

bool GLWidget::hitSphere(const ray &r, const sphere &s, float &t0, float &t1)
{
    float discriminant, a, b, c;
    a = (r.dir).Dot(r.dir);
    b = ((r.start-s.pos).Dot(r.dir))*2;
    c = (r.start - s.pos).Dot((r.start - s.pos))-(s.size)*(s.size);
    discriminant = b*b - 4*a*c;

    if(discriminant < 0)
        return false;
    else
    {
        float q;
        if(b < 0)
            q = (-b + sqrt(discriminant))/2.0f;
        else
            q = (-b - sqrt(discriminant))/2.0f;
        t0 = q/a;
        t1 = c/q;
        return true;
    }
}

void GLWidget::normalize(CVector3D &v)
{
    float magnitude = v.Magnitude();
    v /= magnitude;
}

void GLWidget::initializeGL()
{
    //Background color will be white
    glClearColor(1.0, 1.0, 1.0, 1.0);
    object = makeObject();
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_DEPTH_TEST);
    GLubyte image[pictureSize][pictureSize][3];

    //set up eye position. you may adjust the location by modifying eyePositionX, eyePositionY, eyePositionZ
    //note: z is positive going into the screen and negative going out of the screen
    CVector3D eye = CVector3D(eyePositionX, eyePositionY, eyePositionZ);

    //limit the number reflections. if refectionLevelMax is set to zero or a negative number, then no reflection will be computed
    int refectionLevelMax = 1;  //have 2 levels of reflections for now... you may change it if you want

    //create material
    //material = {reflectionCoefficient, rgb, diffussPower, specularPower}
    //note: diffussPower and specularPower are used to control the size of the bright areas
    material m0 = {0.3f, 1.0f, 0.0f, 0.0f,      5.0f,   128.0f};
    material m1 = {0.4f, 1.0f, 0.8431f, 0.0f,   3.0f,   128.0f};
    material m2 = {0.4f, 0.4f, 0.90f, 0.90f,    2.0f,   128.0f};
    material m3 = {0.4f, 0.3f, 0.2f, 0.6f,      2.0f,   64.0f};

    //create spheres
    //note: (pictureSize/2.0f) is the center. if pictureSize is even then the center is off by one
    //sphere = {location, radius, material}
    sphere s0 = {CVector3D(pictureSize/2.0f, pictureSize/2.0f-50.0f, 200.0f),           100.0f, m0};
    sphere s1 = {CVector3D(pictureSize/2.0f-200.0f, pictureSize/2.0f+100, 400.0f),      180.0f, m2};
    sphere s2 = {CVector3D(pictureSize/2.0f+200.0f, pictureSize/2.0f-10, 100.0f),       50.0f,  m1};
    sphere s3 = {CVector3D(pictureSize/2.0f-250.0f, pictureSize/2.0f-60, 150.0f),       50.0f,  m3};
    sphere s4 = {CVector3D(pictureSize/2.0f+100, pictureSize/2.0f, 50.0f),              30.0f,  m2};

    //add five spheres to the scene
    vector<sphere> spheres;
    spheres.push_back(s0);
    spheres.push_back(s1);
    spheres.push_back(s2);
    spheres.push_back(s3);
    spheres.push_back(s4);

    //create light sources
    //light = {location, rgb}
    light l0 = {CVector3D(pictureSize, pictureSize/3, -100.0f),     1.0f, 1.0f, 1.0f};
    light l1 = {CVector3D(pictureSize/2, pictureSize/2+100, 0.0f),  1.0f, 1.0f, 1.0f};

    //add two light sources to the scene
    vector<light> lights;
    lights.push_back(l0);
    lights.push_back(l1);

    //create texture...i = height, j = width
    int i,j;
    for(i = 0; i < pictureSize; i++)
    {
        for(j = 0; j < pictureSize; j++)
        {
            //rgb initialized to 0, 0, 0
            float red = 0.0f, green = 0.0f, blue = 0.0f;

            //flag for ambient light...we only want to add it once for each point
            bool addAmbientLight = true;

            //for each light source, compute refection...
            for(unsigned int lightIndex = 0; lightIndex < lights.size(); lightIndex++)
            {
                //set up initial view ray, v
                ray viewRay, lightRay;
                viewRay.start = eye;
                viewRay.dir = (CVector3D(j, i, 0.0f) - eye);
                normalize(viewRay.dir);

                //reset reflection coefficient and refectionLevel
                float reflectionCoef = 1.0f;
                int refectionLevel = 0;

                //do....while(refectionLevel < refectionLevelMax)
                do
                {
                    //t0, t1 are used for computing the closest intersection
                    float t0 = 0.0f, t1 = 0.0f, tNear = 0.0f;

                    //find the closest intersection
                    bool firstIteration = true;
                    int indexClosestSphere = -1;
                    for(unsigned int index0 = 0; index0 < spheres.size(); index0++)
                    {
                        if(hitSphere(viewRay, spheres.at(index0), t0, t1) && getMin(t0, t1) > 0.0f)
                        {
                            if(firstIteration)  //keep track of the closest intersection
                            {
                                tNear = getMin(t0, t1);  //initialize tNear
                                firstIteration = false;
                            }
                            if(getMin(t0, t1) <= tNear)
                            {
                                tNear = getMin(t0, t1);
                                indexClosestSphere = index0;
                            }
                        }
                    }
                    //the closest sphere is now "(spheres.at(indexClosestSphere))", break if viewRay misses
                    if(indexClosestSphere == -1)
                        break;

                    //restore t0, t1 value
                    hitSphere(viewRay, spheres.at(indexClosestSphere), t0, t1);

                    //p, n, l, r
                    CVector3D pointOfIntersection, normal, lightDirection, reflectDirection;

                    //compute p
                    pointOfIntersection = viewRay.start + (getMin(t0, t1)*(viewRay.dir));

                    //compute n
                    normal = (pointOfIntersection - (spheres.at(indexClosestSphere)).pos);
                    normalize(normal);

                    //compute l
                    lightDirection = lights.at(lightIndex).pos - pointOfIntersection;
                    normalize(lightDirection);
                    lightRay.start = pointOfIntersection;
                    lightRay.dir = lightDirection;

                    //compute r
                    reflectDirection = ((2*(lightRay.dir).Dot(normal))*normal) - lightDirection;
                    normalize(reflectDirection);

                    //check if current light source is on the opposite side by testing the sign of dot product
                    if((lightRay.dir).Dot(normal) > 0.0f)
                    {
                        //check if point is in shadow
                        bool inShadow = false;
                        for(unsigned index1 = 0; index1 < spheres.size(); index1++)
                        {
                            if(hitSphere(lightRay, spheres.at(index1), t0, t1) && getMin(t0, t1) > 0.0f)
                            {
                                inShadow = true;
                                break;
                            }
                        }
                        //if not in shadow then accumulate diffuss and specular
                        if(!inShadow)
                        {
                            float diffuse = 0.0f, specular = 0.0f;
                            CVector3D specularColorControl = CVector3D(1.0f, 1.0f, 0.5f);   //bright spot will be yellow instead of white

                            //diffuse
                            diffuse = pow((lightRay.dir).Dot(normal), (spheres.at(indexClosestSphere)).mat.diffussPower);
                            red += reflectionCoef * diffuse * lights.at(lightIndex).red * (spheres.at(indexClosestSphere)).mat.red;
                            green += reflectionCoef * diffuse * lights.at(lightIndex).green * (spheres.at(indexClosestSphere)).mat.green;
                            blue += reflectionCoef * diffuse * lights.at(lightIndex).blue * (spheres.at(indexClosestSphere)).mat.blue;

                            //specular
                            specular = pow(getMax((-viewRay.dir).Dot(reflectDirection), 0.0f), (spheres.at(indexClosestSphere)).mat.specularPower);
                            red += reflectionCoef * specular * lights.at(lightIndex).red * specularColorControl.x;
                            green += reflectionCoef * specular * lights.at(lightIndex).green * specularColorControl.y;
                            blue += reflectionCoef * specular * lights.at(lightIndex).blue * specularColorControl.z;
                        }
                    }
                    //accumulate ambient light(only once for each point)
                    if(addAmbientLight)
                    {
                        CVector3D ambientColorControl = CVector3D(0.10f, 0.10f, 0.10f);
                        red += ambientColorControl.x * (spheres.at(indexClosestSphere)).mat.red;
                        green += ambientColorControl.y * (spheres.at(indexClosestSphere)).mat.green;
                        blue += ambientColorControl.z * (spheres.at(indexClosestSphere)).mat.blue;
                        addAmbientLight = false;
                    }
                    //compute new view ray for reflection
                    viewRay.start = pointOfIntersection;
                    viewRay.dir = reflectDirection;
                    reflectionCoef *= (spheres.at(indexClosestSphere)).mat.reflection;
                    refectionLevel++;
                }while(refectionLevel < refectionLevelMax+1);
            }
            //save the color into buffer
            image[i][j][0] = (unsigned char)getMin(red*255.0f, 255.0f);
            image[i][j][1] = (unsigned char)getMin(green*255.0f, 255.0f);
            image[i][j][2] = (unsigned char)getMin(blue*255.0f, 255.0f);
        }//for
    }//for
/*
    //Make sure the output image is not upsidedown
    GLubyte temp[pictureSize][pictureSize][3];
    for(int i = 0; i < pictureSize; i++)
    {
        for(int j = 0; j < pictureSize; j++)
        {
            temp[pictureSize-1-i][j][0] = image[i][j][0];
            temp[pictureSize-1-i][j][1] = image[i][j][1];
            temp[pictureSize-1-i][j][2] = image[i][j][2];
        }
    }

    //output image to file
    QImage output(((unsigned char *) temp), pictureSize, pictureSize, QImage::Format_RGB888);
    output.save("output.bmp", 0, -1 );
*/
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, pictureSize, pictureSize, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glScaled(2.0, 2.0, 2.0);
    glCallList(object);
}

void GLWidget::resizeGL(int width, int height)
{
    //Set projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(Left, Right, Bottom, Top, Near, Far);

    //Return to ModelView Matrix
    glMatrixMode(GL_MODELVIEW);
    glViewport(0, 0, width, height);
}

void GLWidget::makePolygon(int a, int b, int c, int d)
{
    glBegin(GL_POLYGON);
        glTexCoord2f(0.0, 0.0);
        glVertex3fv(vertices[a]);
        glTexCoord2f(0.0, 1.0);
        glVertex3fv(vertices[b]);
        glTexCoord2f(1.0, 1.0);
        glVertex3fv(vertices[c]);;
        glTexCoord2f(1.0, 0.0);
        glVertex3fv(vertices[d]);
    glEnd();
}

GLuint GLWidget::makeObject()
{
    //Creating and compiling a list of objects
    GLuint list = glGenLists(1);
    glNewList(list, GL_COMPILE);
    glShadeModel(GL_SMOOTH);
    makePolygon(0, 1, 2, 3);
    glEndList();
    return list;
}

float GLWidget::getMin(const float &n1, const float &n2)
{
    if(n1 > n2)
        return n2;
    else
        return n1;
}

float GLWidget::getMax(const float &n1, const float &n2)
{
    if(n1 < n2)
        return n2;
    else
        return n1;
}

