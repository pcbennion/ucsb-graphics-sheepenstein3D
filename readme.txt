CS180: Computer Graphics

Projects 1 - 3 : Sheepenstein 3D

Submission is designed to be run on the computers at
UCSB's CSIL tech lab. Makefile and includes are for a
Linux Fedora distribution, with GLU and GLUT installed.

This is a simple scanline graphics program. It takes an 
image from 'input.ppm' and renders a small maze. As
submitted, 'input.ppm' is a sample image of sheep
grazing in a field.

The controls are as follows:
- i and k: forward and backward, respectively.
- j and l: turn left and right.
- J and L: strafe left and right.
- SPACEBAR: toggle render mode - OpenGL and software.
- q and ESC: quit.

The program runs in two modes: OpenGL, where all
rendering is handled by integrated hardware; and
software-only, where the entire graphics pipeline is 
simulated within the program itself.

The software-only rendering mode consists of custom
functions to replace GLULookAt and GLUPerspective; a
rendering function to perform all necessary transforms
on the triangles in the scene; and a final drawTriangle
function for determining the colors of pixels that
comprise a textured triangle. drawTriangle uses
homogenous barycentric coordinates to correct for the
depth distortions that arise from  the frustrum 
transformations. In this mode, all objects that lie
partially outside the view frustrum are culled for
simplicity's sake.

Z-buffering only exists in the software-only mode. In
OpenGL, all objects will always display in draw order.