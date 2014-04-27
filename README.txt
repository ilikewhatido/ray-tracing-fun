This ray tracer is built on top of a opensource texturing framework. The area 
of interest is the "initializeGL" function in "glwidget.cpp". All the main raytraicng code is in this 
function. I stored the rendered scene as a texture file and used OpenGL to display it. There is also 
an output image called “output.bmp” which can be found in the ”debug” folder. The texture being
displayed and the output image should be identical to each other. 

The followings are structures used in my ray tracer. They can be found in "glwidget.h". Ideally these 
structures should be converted into classes. But since my ray tracer is a simple one there is no need 
to make things more complicated. 

struct material
{
    float reflection;
    float red, green, blue;
    float diffussPower, specularPower;
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

My ray tracer supports the following features:

- ambient light
- specular light
- diffuse light
- refection 
- hard shadow
- multiple spheres
- multiple light sources
- modifiable sampling size
- modifiable view point

I left lots of comments in the "initializeGL" function. Please refer to line 69 in "glwidget.cpp" 
for more details. 
