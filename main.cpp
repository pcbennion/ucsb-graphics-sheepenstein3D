
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "MatrixVectorOps.h"

using namespace std;

#define PI 3.14159265

// picture info
GLuint texture;
int height, width;
unsigned char* picture;

// parameters for GL functions
float zNear = 0.01, zFar = 10;

// coord array for walls
float wall[50][4][4];
float* buffer; // memory allocated storage for vector operations

// view info
float camerapos[3];
float yRotation;
int viewheight, viewwidth;
unsigned char* pixels;
float* zBuffer;

// matricies
float* modelview;
float* projection;
float* viewport;
float* fVector; float* sVector; float* uVector; float* upVector;

// keys held down
int i, j, k, l, J, L;
int OGL; // spacebar toggle

static unsigned int maze[] = {
 1, 0, 1, 1, 0,
1, 0, 1, 1, 1,
 0, 0, 0, 0, 0,
1, 1, 1, 0, 1,
 0, 1, 0, 0, 0,
1, 0, 0, 1, 1,
 1, 1, 1, 0, 0,
0, 0, 0, 0, 1,
 1, 1, 1, 1, 0,
0, 0, 0, 0, 0
};

//////////////////////////////////////////////////////////////////////////////////
// Sofware version of gluLookAt
//////////////////////////////////////////////////////////////////////////////////
void myLookAt(GLdouble eyeX,  GLdouble eyeY,  GLdouble eyeZ,  GLdouble centerX,  GLdouble centerY,  GLdouble centerZ,  GLdouble upX,  GLdouble upY,  GLdouble upZ)
{
    float mag; // magnitude, for normalization

    // initialize matrix
    getIdentityMatrix(modelview);

    // Begen calculations. From opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml
    fVector[0] = centerX-eyeX; fVector[1] = centerY-eyeY; fVector[2] = centerZ-eyeZ, fVector[3] = 1; // f vector
    mag = sqrt((fVector[0]*fVector[0])+(fVector[1]*fVector[1])+(fVector[2]*fVector[2]));
    fVector[0] = fVector[0]/mag; fVector[1] = fVector[1]/mag; fVector[2] = fVector[2]/mag; // normalize
    modelview[8] = -fVector[0]; modelview[9] = -fVector[1]; modelview[10] = -fVector[2]; // put in matrix

    upVector[0] = upX; upVector[1] = upY; upVector[2] = upZ; upVector[3] = 1; // up vector
    mag = sqrt((upVector[0]*upVector[0])+(upVector[1]*upVector[1])+(upVector[2]*upVector[2]));
    upVector[0] = upVector[0]/mag; upVector[1] = upVector[1]/mag; upVector[2] = upVector[2]/mag; // normalize

    sVector[0] = fVector[0]; sVector[1] = fVector[1]; sVector[2] = fVector[2]; sVector[3] = 1;
    vectorCrossProduct(sVector, upVector); // compute s vector
    modelview[0] = sVector[0]; modelview[1] = sVector[1]; modelview[2] = sVector[2]; // put in matrix

    mag = sqrt((sVector[0]*sVector[0])+(sVector[1]*sVector[1])+(sVector[2]*sVector[2]));
    uVector[0] = sVector[0]/mag; uVector[1] = sVector[1]/mag; uVector[2] = sVector[2]/mag; uVector[3] = 1; // normalize s vector for first part of u vector
    vectorCrossProduct(uVector, fVector); // finish computing u vector
    modelview[4] = -uVector[0]; modelview[5] = -uVector[1]; modelview[6] = -uVector[2]; // put in matrix

    // put camera coords into matrix, translated into global coords
    modelview[3] = -camerapos[2] * sin((PI/180) * yRotation) + camerapos[0] * cos((PI/180) * yRotation);
    modelview[7] = -camerapos[1];
    modelview[11] = camerapos[0] * sin((PI/180) * yRotation) + camerapos[2] * cos((PI/180) * yRotation);

    // convert matrix to column major format
    //matrixTranspose(modelview);

    /* DEBUG: print calculated matrix and current one in opengl
    gluLookAt(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ);
    matrixPrint(modelview, "Calculated Modelview");
    float oglmatrix[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, oglmatrix);
    cout<<"OpenGL:"<<endl;
    int i, j;
    int index;
    for(i=0; i<4; i++) {
        for(j=0; j<4; j++) {
            index = 4*i+j;
            std::cout<<oglmatrix[index];
            if((4*i+j)<15) std::cout<<", ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;*/

    // Load finished matrix into OpenGL
    //glLoadMatrixf(modelview);
    //matrixTranspose(modelview); // flip it back to row-major
}

//////////////////////////////////////////////////////////////////////////////////
// Sofware version of gluprojection
//////////////////////////////////////////////////////////////////////////////////
void myPerspective(GLdouble  fovy,  GLdouble  aspect,  GLdouble  zNear,  GLdouble  zFar)
{
    float f = cos(fovy/2 * PI/180)/sin(fovy/2 * PI/180); // cotangent = cos/sin

    // Compute matrix. Values from opengl.org/sdk/docs/man2/xhtml/gluprojection.xml
    getIdentityMatrix(projection);
    projection[0] = f / aspect;
    projection[5] = f;
    projection[10] = (zFar+zNear)/(zNear-zFar);
    projection[11] = (2*zFar*zNear)/(zNear-zFar);
    projection[14] = -1;
    projection[15] = 0;

    // convert matrix to column major format
    //matrixTranspose(projection);

    /*DEBUG: print calculated matrix and current one in opengl
    gluPerspective(fovy, aspect, zNear, zFar);
    matrixPrint(projection, "Calculated Projection");
    float oglmatrix[16];
    glGetFloatv(GL_PROJECTION_MATRIX, oglmatrix);
    cout<<"OpenGL:"<<endl;
    int i, j;
    int index;
    for(i=0; i<4; i++) {
        for(j=0; j<4; j++) {
            index = 4*i+j;
            std::cout<<oglmatrix[index];
            if((4*i+j)<15) std::cout<<", ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;*/

    // Load finished matrix into OpenGL
    //glLoadMatrixf(projection);
    //matrixTranspose(projection); // flip it back to row-major
}

//////////////////////////////////////////////////////////////////////////////////
// Initializes GL options and loads image to OpenGL
//////////////////////////////////////////////////////////////////////////////////
void init(void)
{
   glClearColor (0.0, 0.0, 0.0, 0.0);
   glShadeModel(GL_FLAT);
   glEnable(GL_DEPTH_TEST);

   glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

   glGenTextures(1, &texture);
   glBindTexture(GL_TEXTURE_2D, texture);

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                   GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                   GL_NEAREST);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width,
                height, 0, GL_RGB, GL_UNSIGNED_BYTE,
                picture);
}

//////////////////////////////////////////////////////////////////////////////////
// In case the window gets resized
//////////////////////////////////////////////////////////////////////////////////
void resize (int w, int h)
{
    if (h==0)  h=1;// don't divide by 0
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    viewheight = h; viewwidth = w;
    viewport[0] = w/2; viewport[3] = (w)/2; viewport[5] = h/2; viewport[7] = (h)/2;
}

//////////////////////////////////////////////////////////////////////////////////
// Rendering
//////////////////////////////////////////////////////////////////////////////////
void drawTriangle(float* a, float au, float av, float* b, float bu, float bv, float* c, float cu, float cv, char red, char green, char blue)
{
    // find min and max
    int xmin, xmax, ymin, ymax;
    float zval;
    float f;
    f = 100000; f = ((a[0]<f) ? a[0] : f); f = ((b[0]<f) ? b[0] : f); f = ((c[0]<f) ? c[0] : f);
    xmin = (int)f;
    f = 100000; f = ((a[1]<f) ? a[1] : f); f = ((b[1]<f) ? b[1] : f); f = ((c[1]<f) ? c[1] : f);
    ymin = (int)f;
    f = 0; f = ((a[0]>f) ? a[0] : f); f = ((b[0]>f) ? b[0] : f); f = ((c[0]>f) ? c[0] : f);
    xmax = (int)f;
    f = 0; f = ((a[1]>f) ? a[1] : f); f = ((b[1]>f) ? b[1] : f); f = ((c[1]>f) ? c[1] : f);
    ymax = (int)f;

    // Loop through all pixels in bounding box
    float alpha, beta, gamma;
    float alphaW, betaW, gammaW, d;
    int u, v;
    for(int x=xmin; x<=xmax; x++)
    {
        for(int y=ymin; y<=ymax; y++)
        {
            // get barycentric coords
           alpha = ((c[1]-b[1])*x    + (b[0]-c[0])*y    + c[0]*b[1] - b[0]*c[1])
                  /((c[1]-b[1])*a[0] + (b[0]-c[0])*a[1] + c[0]*b[1] - b[0]*c[1]);
            beta = ((a[1]-c[1])*x    + (c[0]-a[0])*y    + a[0]*c[1] - c[0]*a[1])
                  /((a[1]-c[1])*b[0] + (c[0]-a[0])*b[1] + a[0]*c[1] - c[0]*a[1]);
           gamma = ((a[1]-b[1])*x    + (b[0]-a[0])*y    + a[0]*b[1] - b[0]*a[1])
                  /((a[1]-b[1])*c[0] + (b[0]-a[0])*c[1] + a[0]*b[1] - b[0]*a[1]);

           // find average z of triangle
           // because no intersecting triangles are allowed, using this value for all pixels in
           // triangle will still let z-buffer work properly
           zval = alpha*a[2] + beta*b[2] + gamma*c[2];

            // if within triangle, draw pixel
            if(alpha>0 && beta>0 && gamma>0)
            {
                if(zval<zBuffer[viewwidth*y + x])
                {
                    //cout<<a[3]<<", "<<b[3]<<", "<<c[3]<<endl;

                    // homogenize barycentric coords
                    d = b[3]*c[3] + c[3]*beta*(a[3]-b[3]) + b[3]*gamma*(a[3]-c[3]);
                    betaW = (a[3]*c[3]*beta)/d;
                    gammaW = (a[3]*b[3]*gamma)/d;
                    alphaW = 1-betaW-gammaW;

                    //cout<<d<<", "<<alphaW<<", "<<betaW<<", "<<gammaW<<endl;

                    // calculate texture coords
                    u = (int)((au*alphaW + bu*betaW + cu*gammaW)*width);
                    v = (int)((av*alphaW + bv*betaW + cv*gammaW)*height);

                    //cout<<u<<", "<<v<<endl;

                    // set pixels
                    pixels[3*viewwidth*y + 3*x + 0] = picture[3*width*v + 3*u + 0];
                    pixels[3*viewwidth*y + 3*x + 1] = picture[3*width*v + 3*u + 1];
                    pixels[3*viewwidth*y + 3*x + 2] = picture[3*width*v + 3*u + 2];

                    /*
                    // for solid color only
                    //cout<<"Draw: "<<x<<", "<<y<<endl;
                    pixels[3*viewwidth*y + 3*x + 0] = red;
                    pixels[3*viewwidth*y + 3*x + 1] = green;
                    pixels[3*viewwidth*y + 3*x + 2] = blue;*/

                    zBuffer[viewwidth*y + x] = zval;
                }
            }
        }
    }
}

void OGLRender()
{
    glEnable(GL_TEXTURE_2D);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glBindTexture(GL_TEXTURE_2D, texture);
    for(int index=0; index<50; index++)
    {
        if(maze[index] != 0)
        {
            glBegin(GL_TRIANGLES);
                  glColor3f(1.0, 0.0, 0.0); // red
                  glTexCoord2f(1, 0); glVertex3f(wall[index][0][0], wall[index][0][1], wall[index][0][2]);
                  glTexCoord2f(1, 1); glVertex3f(wall[index][1][0], wall[index][1][1], wall[index][1][2]);
                  glTexCoord2f(0, 0); glVertex3f(wall[index][2][0], wall[index][2][1], wall[index][2][2]);

                  glColor3f(0.0, 1.0, 0.0); // green
                  glTexCoord2f(1, 1); glVertex3f(wall[index][1][0], wall[index][1][1], wall[index][1][2]);
                  glTexCoord2f(0, 0); glVertex3f(wall[index][2][0], wall[index][2][1], wall[index][2][2]);
                  glTexCoord2f(0, 1); glVertex3f(wall[index][3][0], wall[index][3][1], wall[index][3][2]);
            glEnd();
        }
    }
    glFlush();
    glDisable(GL_TEXTURE_2D);

    // done with scene
    glutSwapBuffers();
}

// software-only rendering
void myRender()
{
    // clear pixel buffer
    int index;
    for(index=0; index<viewheight*viewwidth*3; index++) pixels[index] = 0;
    for(index=0; index<viewheight*viewwidth; index++) zBuffer[index] = 100000;

    int valid[4]; // tells renderer if triangle was culled
    float* a = new float[4]; float* b = new float[4]; float* c = new float[4]; float* d = new float[4];
    float* vectorptr;
    float homogeneous;
    //glEnable(GL_TEXTURE_2D);
    //glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    //glBindTexture(GL_TEXTURE_2D, texture);
    for(int index=0; index<50; index++)
    {
        if(maze[index] != 0)
        {
            // PERSPECTIVE DIVIDE:
            // check each vertex for this wall in sequence
            valid[0] = 1; valid[1] = 1; valid[2] = 1; valid[3] = 1;
            for(int vert = 0; vert<4; vert++)
            {
                // populate buffer for matrix multiplication
                buffer[0] = wall[index][vert][0]; buffer[1] = wall[index][vert][1]; buffer[2] = wall[index][vert][2]; buffer[3] = wall[index][vert][3]; buffer[4] = wall[index][vert][4];
                //if(index==0&&vert==0) vectorPrint(buffer, "original");

                // multiply by modelview
                multiplyMatrixVector(modelview, buffer);
                //if(index==0&&vert==0) matrixPrint(modelview, "modelview");
                //if(index==0&&vert==0) vectorPrint(buffer, "after modelview multiply");

                // multiply by projection
                multiplyMatrixVector(projection, buffer);
                //if(index==0&&vert==0) matrixPrint(projection, "projection");
                //if(index==0&&vert==0) vectorPrint(buffer, "after projection multiply");

                // Z-culling
                if(buffer[3]<zNear || buffer[3]>zFar)
                {
                    valid[vert] = 0;
                    //if(index==0) cout<<"Z-Culled: "<<vert<<endl;
                }

                // perspective divide
                homogeneous = buffer[3];
                buffer[0] = buffer[0]/buffer[3];
                buffer[1] = buffer[1]/buffer[3];
                buffer[2] = buffer[2]/buffer[3];
                buffer[3] = buffer[3]/buffer[3];

                //if(index==0&&vert==0) vectorPrint(buffer, "after perspective divide");

                // X/Y-culling
                if(buffer[0]<-1 || buffer[0]>1 || buffer[1]<-1 || buffer[1]>1)
                {
                    valid[vert] = 0;
                    //if(index==0) cout<<"XY-Culled: "<<vert<<endl;
                }

                // viewport transform
                multiplyMatrixVector(viewport, buffer);
                //if(index==0&&vert==0) matrixPrint(viewport, "viewport");
                //if(index==0&&vert==0) vectorPrint(buffer, "after viewport multiply");

                buffer[3] = homogeneous;

                // save results
                switch (vert)
                {
                    case 0: vectorptr = a; break;
                    case 1: vectorptr = b; break;
                    case 2: vectorptr = c; break;
                    case 3: vectorptr = d; break;
                }
                vectorptr[0] = buffer[0]; vectorptr[1] = buffer[1]; vectorptr[2] = buffer[2]; vectorptr[3] = buffer[3];
            }
            if(valid[1] && valid[2])
            {
                drawTriangle(a,1,0,b,1,1,c,0,0,255,0,0);
                drawTriangle(b,1,1,c,0,0,d,0,1,0,255,0);
            }
            /*
            // Draw using gl triangles
            glBegin(GL_TRIANGLES);
            if(valid[1] == 1 && valid[2] == 1) {
                  glColor3f(1.0, 0.0, 0.0); // red
                  glVertex2f(a[0], a[1]);
                  glVertex2f(b[0], b[1]);
                  glVertex2f(c[0], c[1]);
            }

            if(valid[1] == 1 && valid[2] == 1) {
                  glColor3f(0.0, 1.0, 0.0); // green
                  glVertex2f(b[0], b[1]);
                  glVertex2f(c[0], c[1]);
                  glVertex2f(d[0], d[1]);
            }
            glEnd();*/
        }
    }

    delete a, delete b, delete c, delete d;

    // Render picture
    glDrawPixels(viewwidth, viewheight, GL_RGB, GL_UNSIGNED_BYTE, pixels);

    // done with scene
    glutSwapBuffers();
}

void display (void)
{
    // clear buffers
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // matrix stuff
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(OGL == 1) gluPerspective(45, 1, zNear, zFar);
    else myPerspective(45, 1, zNear, zFar);
    glMatrixMode   (GL_MODELVIEW);
    glLoadIdentity ();
    if(OGL == 1) gluLookAt(camerapos[0], camerapos[1], camerapos[2], camerapos[0]+5*sin((PI/180) * yRotation), camerapos[1], camerapos[2]+5*cos((PI/180) * yRotation), 0, 1, 0);
    else myLookAt(camerapos[0], camerapos[1], camerapos[2], camerapos[0]+5*sin((PI/180) * yRotation), camerapos[1], camerapos[2]+5*cos((PI/180) * yRotation), 0, 1, 0);

    if(OGL==1) OGLRender();
    else myRender();
}

//////////////////////////////////////////////////////////////////////////////////
// Handles keyboard events
//////////////////////////////////////////////////////////////////////////////////
void keyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
      // the escape key and 'q' quit the program
      case 27: //esc
      case 'q':
          exit(0);
          break;
      // directional keys
      case 'i':
          i = 1;
          break;
      case 'j':
          j = 1;
          break;
      case 'k':
          k = 1;
          break;
      case 'l':
          l = 1;
          break;
      case 'I':
          i = 1;
          break;
      case 'J':
          J = 1;
          break;
      case 'K':
          k = 1;
          break;
      case 'L':
          L = 1;
          break;
      case 32: //spacebar
          OGL = (OGL == 1) ? 0 : 1;
          break;
  }
}
void keyboardUp(unsigned char key, int x, int y)
{
    switch (key)
    {
        // the escape key and 'q' quit the program
        case 27: //esc
        case 'q':
            exit(0);
            break;
        // directional keys
        case 'i':
            i = 0;
            break;
        case 'j':
            j = 0;
            J = 0;
            break;
        case 'k':
            k = 0;
            break;
        case 'l':
            l = 0;
            L = 0;
            break;
        case 'I':
            i = 0;
            break;
        case 'J':
            J = 0;
            j = 0;
            break;
        case 'K':
            k = 0;
            break;
        case 'L':
            L = 0;
            l = 0;
            break;
    }
}

//////////////////////////////////////////////////////////////////////////////////
// Called occasionally to update rendered geometry
//////////////////////////////////////////////////////////////////////////////////
void idle()
{
    // move camera if the appropriate keys are held down
    float xInc=0, zInc=0; // how much local x and z are incremented
    if(i==1) zInc+=0.035;
    if(k==1) zInc+=-0.035;
    if(J==1) xInc+=0.035;
    if(L==1) xInc+=-0.035;
    if(j==1) yRotation+=2;
    if(l==1) yRotation+=-2;

    // translate incrementations to local x and z coords to global coords
    camerapos[0] += xInc * cos((PI/180) * yRotation) + zInc * sin((PI/180) * yRotation);
    camerapos[2] += - xInc * sin((PI/180) * yRotation) + zInc * cos((PI/180) * yRotation);

    // call display to render changes
    glutPostRedisplay();
}

//////////////////////////////////////////////////////////////////////////////////
// Read in a raw PPM file of the "P6" style.
//
// Input: "filename" is the name of the file you want to read in
// Output: "pixels" will point to an array of pixel values
//         "width" will be the width of the image
//         "height" will be the height of the image
//
// The PPM file format is:
//
//   P6
//   <image width> <image height>
//   255
//   <raw, 8-bit binary stream of RGB values>
//
// Open one in a text editor to see for yourself.
//
//////////////////////////////////////////////////////////////////////////////////
void readPPM(const char* filename, unsigned char*& pixels, int& width, int& height)
{
  // try to open the file
  FILE* file;
  file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << " Couldn't open file " << filename << "! " << endl;
    exit(1);
  }

  // read in the image dimensions
  fscanf(file, "P6\n%d %d\n255\n", &width, &height);
  int totalPixels = width * height;

  // allocate three times as many pixels since there are R,G, and B channels
  pixels = new unsigned char[3 * totalPixels];
  fread(pixels, 1, 3 * totalPixels, file);
  fclose(file);

  // output some success information
  cout << " Successfully read in " << filename << " with dimensions: "
       << width << " " << height << endl;
}

int main (int argc, char **argv)
{
    readPPM("input.ppm", picture, width, height);

    // Initialize walls
    int index;
    for(int x=0; x<5; x++)
    {
        for(int z=0; z<5; z++)
        {
            // Wall along x axis
            index = x+(z*10);
            wall[index][0][0] = x; wall[index][0][1] = 1.0f; wall[index][0][2] = z; wall[index][0][3] = 1;
            wall[index][1][0] = x; wall[index][1][1] = 0.0f; wall[index][1][2] = z; wall[index][1][3] = 1;
            wall[index][2][0] = x+1; wall[index][2][1] = 1.0f; wall[index][2][2] = z; wall[index][2][3] = 1;
            wall[index][3][0] = x+1; wall[index][3][1] = 0.0f; wall[index][3][2] = z; wall[index][3][3] = 1;
            // Wall along z axis
            index += 5;
            wall[index][0][0] = x; wall[index][0][1] = 1.0f; wall[index][0][2] = z; wall[index][0][3] = 1;
            wall[index][1][0] = x; wall[index][1][1] = 0.0f; wall[index][1][2] = z; wall[index][1][3] = 1;
            wall[index][2][0] = x; wall[index][2][1] = 1.0f; wall[index][2][2] = z+1; wall[index][2][3] = 1;
            wall[index][3][0] = x; wall[index][3][1] = 0.0f; wall[index][3][2] = z+1; wall[index][3][3] = 1;
        }
    }
    buffer = new float[4];

    // Initialize camera position and rotation
    camerapos[0] = 0.5, camerapos[1] = 0.5, camerapos[2] = 3.5;
    yRotation = 90.0f;

    // matrix stuff
    modelview = new float[16];
    projection = new float[16];
    viewport = new float[16];
    fVector = new float[4]; sVector = new float[4]; uVector = new float[4]; upVector = new float[4];
    getIdentityMatrix(viewport);
    viewheight = 800; viewwidth = 800;
    viewport[0] = viewwidth/2; viewport[3] = (viewwidth)/2; viewport[5] = viewheight/2; viewport[7] = (viewheight)/2;
    OGL = 1;

    // pixel info
    pixels = new unsigned char[viewheight*viewwidth*3];
    for(index=0; index<viewheight*viewwidth*3; index++) pixels[index] = 0;
    zBuffer = new float[viewheight*viewwidth];
    for(index=0; index<viewheight*viewwidth; index++) zBuffer[index] = 100000;

    // glut loop
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
   glutInitWindowSize(800, 800);
   glutInitWindowPosition(100, 100);
   glutCreateWindow("CMPSC 180, Homework 1");
    init();

   glutDisplayFunc(display);
   glutKeyboardFunc(keyboard);
   glutKeyboardUpFunc(keyboardUp);
   glutIdleFunc(idle);
   glutReshapeFunc(resize);

   glutMainLoop();

   return 0;
}

