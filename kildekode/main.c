#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "GL/glut.h"

static int numx, numy;         /// numx*numy = number of nodes in equation system
static float **fvals;          /// solution of equation
static float min, max;         /// min and max of fvals
static float h;                /// grid step size
static GLuint display_list;    /// display list for scene

GLfloat rotV=0.0f;             /// defines rotation angle about vertical axis
GLfloat rotH=0.0f;             /// defines rotation angle about a horizontal axis
GLfloat vspeed=0.0f;           /// vertical rotation speed
GLfloat hspeed=0.0f;           /// horizontal rotation speed
GLfloat scale=1.0f;            /// scaling variable

void readFile(char* fileName) {
    /***************************************************************************************
    * PURPOSE:
    * Read file to initialize global variables 'numx', 'numy' and 'h', and
    * read function values into global matrix 'fvals'
    *
    * Date/version: 03.05.2019/1.0
    ***************************************************************************************/
	FILE* file;
	file = fopen(fileName, "r");
	if (!file) printf("Data file not found: %s", fileName);

	fscanf(file," %i %i", &numx, &numy); ///read numx and numy in from first line of file
    h = 1.0/(numx-1);                    ///set grid step-size h

	///Allocate memory to matrix 'fvals':
	///RETRIVED FROM: http://pleasemakeanote.blogspot.com/2008/06/2d-arrays-in-c-using-malloc.html
	fvals = (float**) malloc((numx)*sizeof(float*));
    for (int i = 0; i < numx; i++){
        fvals[i] = (float*) malloc((numx)*sizeof(float));
    }

	/**Read function values into matrix.
	* Function values are stored in file as a list.
	* Each line contains the calculated function value at point (x,y), starting at (0,0).
	* Function values are organized in COLUMN-MAJOR order in file.
	**/
	for(int n=0; n < numx; n++) {               /// Loop over columns first,
        for(int m=0; m < numy; m++) {           /// then loop over rows.
            fscanf(file,"%f", &fvals[n][m]);       /// Read values directly into matrix
                                                /// by exploiting ordering of file
            ///find max and min function values:
            float z=fvals[n][m];
            if ((n==0) & (m==0)) {
                max=z;
                min=z;
            }
            else if (z>max){
                max=z;
            }
            else if (z<min){
                min=z;
            }

        }
	}

	fclose(file);

}

void renderBitmapString(
    /*******************
    * PURPOSE:
    * Render a string starting at specified raster position
    *
    * RETRIEVED FROM: https://www.lighthouse3d.com/tutorials/glut-tutorial/bitmap-fonts/
    *
    *******************/
  float x,            /// x position
  float y,            /// y position
  float z,            /// z position
  void *font,         /// 'font' is chosen font
  char *string) {     /// 'string' string to render
  char *c;
  glRasterPos3f(x, y,z);
  for (c=string; *c != '\0'; c++) {
    glutBitmapCharacter(font, *c);
  }
}

void init() {
    /************
    * PURPOSE:
    * Setup for OpenGL.
    * Compile surface plot, grid and x-, y-, z-axes for later execution.
    *
    * Date/version: 03.05.2019/1.0
    *************/
    glClearColor(1.0, 1.0, 1.0, 0);                     /// set background to white
    glEnable(GL_BLEND);                                 /// enables blending
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  /// alpha blending for transparency

    ///set up axes and scene, and store result in display list:
    glColor3f(1.0,1.0,1.0);
    display_list = glGenLists(1);
    glNewList(display_list, GL_COMPILE);
        ///make surface:
        glColor4f(0,255,255,0.5);           ///cyan, somewhat transparent
        for (int i=0; i<numx-1; i++)    {
            for (int j=0; j<numy-1; j++)    {
                glBegin(GL_POLYGON);
                float x = i*h;              ///define proper x-,
                float y = j*h;              ///y-,
                float z = fvals[i][j];      ///and z-coordinates
                glVertex3f(x,y,z);
                x = i*h;
                y = (j+1)*h;
                z = fvals[i][j+1];
                glVertex3f(x,y,z);
                x = (i+1)*h;
                y = (j+1)*h;
                z = fvals[i+1][j+1];
                glVertex3f(x,y,z);
                x = (i+1)*h;
                y = (j)*h;
                z = fvals[i+1][j];
                glVertex3f(x,y,z);
            glEnd();
        }
    }

        ///make grid on surface:
        glColor3f(0.0f,0.0f,1.0f);             ///set grid color to pure blue
        for (int i=0; i<numx-1; i++)
    {
            for (int j=0;j<numy-1; j++)
        {
                glBegin(GL_LINE_LOOP);
                    float x = i*h;
                    float y = j*h;
                    float z = fvals[i][j];
                    glVertex3f(x,y,z);
                    x = i*h;
                    y = (j+1)*h;
                    z = fvals[i][j+1];
                    glVertex3f(x,y,z);
                    x = (i+1)*h;
                    y = (j+1)*h;
                    z = fvals[i+1][j+1];
                    glVertex3f(x,y,z);
                    x = (i+1)*h;
                    y = (j)*h;
                    z = fvals[i+1][j];
                    glVertex3f(x,y,z);
                glEnd();
        }
    }

        /// display values at (1,1,z), (0,1,z), (1,0,z)
        glColor3f(0.0f,0.0f,0.0f);
        char buffer[10]={'\0'};
        sprintf(buffer, "z = %f", fvals[numx-1][numy-1]);
        renderBitmapString(1.0f, 1.0f , fvals[numx-1][numy-1]*(1.0f+0.01f), GLUT_BITMAP_HELVETICA_12, buffer);
        sprintf(buffer, "z = %f", fvals[0][numy-1]);
        renderBitmapString(-0.22f, 1.22f , fvals[0][numy-1]*(1.0f+0.01f), GLUT_BITMAP_HELVETICA_12, buffer);
        sprintf(buffer, "z = %f", fvals[numx-1][0]);
        renderBitmapString(1.0f, -0.05f , fvals[numx-1][0]*(1.0f+0.01f), GLUT_BITMAP_HELVETICA_12, buffer);

        /// make vertical dashed lines from floor to displayed values:
        glLineStipple(1, 0x00FF);
        glEnable(GL_LINE_STIPPLE);
        glBegin(GL_LINE_STRIP);
        glVertex3f(0.0f,1.0f,0.0f);
        glVertex3f(0.0f,1.0f,fvals[0][numy-1]);
        glEnd();

        glBegin(GL_LINE_STRIP);
        glVertex3f(1.0f,0.0f,0.0f);
        glVertex3f(1.0f,0.0f,fvals[numx-1][0]);
        glEnd();

        glBegin(GL_LINE_STRIP);
        glVertex3f(1.0f,1.0f,0.0f);
        glVertex3f(1.0f,1.0f,fvals[numx-1][numy-1]);
        glEnd();

        glBegin(GL_LINE_LOOP);
        glVertex3f(0.0f,0.0f,0.0f);
        glVertex3f(0.0f,1.0f,0.0f);
        glVertex3f(1.0f,1.0f,0.0f);
        glVertex3f(1.0f,0.0f,0.0f);
        glEnd();

        glDisable(GL_LINE_STIPPLE);

        /// x-axis:
        glColor3f(0, 0, 0);
        glBegin(GL_LINE_STRIP);
            glVertex3f(0.0f, 0.0f, 0.0f);
            glVertex3f(1.3f, 0.0f, 0.0f);
        glEnd();

        glBegin(GL_TRIANGLES);
            glVertex3f(1.25f, 0.0f, 0.01f);
            glVertex3f(1.25f, 0.0f, -0.01f);
            glVertex3f(1.3f, 0.0f, 0.0f);
        glEnd();

        /// y-axis:
        glBegin(GL_LINE_STRIP);
            glVertex3f(0.0f, 0.0f, 0.0f);
            glVertex3f(0.0f, 1.3f, 0.0f);
        glEnd();

        glBegin(GL_TRIANGLES);
            glVertex3f(0.0f, 1.3f, 0.0f);
            glVertex3f(0.0f, 1.25f, 0.01f);
            glVertex3f(0.0f, 1.25f, -0.01f);
        glEnd();

        /// z-axis:
        glBegin(GL_LINE_STRIP);
            glVertex3f(0.0f, 0.0f, 0.0f);
            glVertex3f(0.0f, 0.0f, max);
        glEnd();

        glBegin(GL_TRIANGLES);
            glVertex3f(0.0f, 0.0f, max);
            glVertex3f(0.01f*max, 0.0f, 1.25f/1.3f*max);
            glVertex3f(-0.01f*max, 0.0f, 1.25f/1.3f*max);
        glEnd();

        glBegin(GL_TRIANGLES);
            glVertex3f(0.0f, 0.00f, max);
            glVertex3f(0.0f, -0.01f*max, 1.25f/1.3f*max);
            glVertex3f(0.0f, 0.01f*max, 1.25f/1.3f*max);
        glEnd();

        ///name x-, y- and z-axis:
        glColor3f(0, 0, 0);
        glRasterPos3f(1.25f, 0.05f, 0.0f);
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'x');
        glRasterPos3f(1.05f, 0.075f, 0.0f);
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, '1');
        glRasterPos3f(0.05f, 1.25f, 0.0f);
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'y');
        glRasterPos3f(0.075f, 1.05f, 0.0f);
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, '1');
        glRasterPos3f(0.0f, 0.05f, max);
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, 'z');


    glEndList();

}

void display()  {
    /************
    *PURPOSE:
    * Draws 3D-plot of solution as surface with grid
    *
    * Date/version: 03.05.2019/1.0
    *************/

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); /// clear screen and buffer
    glMatrixMode(GL_MODELVIEW);                         /// use modelview matrix
    glLoadIdentity();                                   /// reset matrix
    gluLookAt(-max, -max, max,                          /// eye position
              0.0f, 0.0f, (max-min)*1.0f/2,             /// reference point
              0.0f, 0.0f, 1.0f);                        /// up vector
    glScalef(scale,scale,scale);                        /// initial scale=1
    glRotatef(rotV,1,-1,0);                             /// initial v. rot.=0
    glRotatef(rotH,0,0,1);                              /// initial h. rot.=0


    ///draw scene:
    glPushMatrix();
    glCallList(display_list);
    glPopMatrix();

    ///increment horizontal/vertical rot. speed (user controllable):
    rotV+=vspeed;
    rotH+=hspeed;

    glutSwapBuffers();
}

void reshape(int w, int h){
    /***************
    * PURPOSE:
    * In case of change of window height/width, this function maintains the correct perspective
    *
    * Date/version: 03.05.2019/1.0
    ****************/

    if (h==0) h=1;                  /// prevents division by zero in next line
    float ratio = w*1.0 / h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45,ratio,0,1000);
    glMatrixMode(GL_MODELVIEW);

}

void keyPress(int key, int xx, int yy) {
    /*********************
    * PURPOSE:
    * Process keypresses to rotate and scale figure
    *
    * Date/version: 03.05.2019/1.0
    **********************/

	switch (key) {
	    case GLUT_KEY_UP:   /// rotate fig. "towards" camera
            vspeed+=0.005f;
            break;
	    case GLUT_KEY_DOWN: /// rotate fig. "away" from camera
            vspeed-=0.005f;
            break;
        case GLUT_KEY_LEFT: /// rotate clockwise about vertical axiz
            hspeed-=0.005f;
            break;
        case GLUT_KEY_RIGHT:/// rotate counterclockwise about vertical axis
            hspeed+=0.005f;
            break;
        case GLUT_KEY_F1:   /// scale down
            scale-=0.05;
            break;
        case GLUT_KEY_F2:   /// scale up
            scale+=0.05;
            break;
	}
}

int main(int argc, char** argv)
{   /**************************
    * PURPOSE:
    * The main function initializes and calls
    * functions in the right order for program to run
    * successfully
    *
    * Date/version: 03.05.2019/1.0
    *******************************/
    readFile("solution.DAT");
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB);
    glutInitWindowSize(800,600);
    glutInitWindowPosition(100,100);
    glutCreateWindow("Visualization of Poisson equation solution");
    glutSpecialFunc(keyPress);      /// process keypresses
    glutDisplayFunc(display);       /// callback function
    glutIdleFunc(display);          /// callback function
    glutReshapeFunc(reshape);       /// callback function

    init();
    glutMainLoop();                 ///event processing cycle

    return 0;
}
