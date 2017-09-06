#define NOMINMAX
#include <windows.h>
#include <gl/GL.h>
#include <gl/GLU.h>
#include <iostream>

#include "tools/constsandtypes.h"
#include "tools/vectortools.h"
#include "vis/ItemDraw_c.h"

void drawSphere(double r, int ndiv = 100, double colorR = 0, double colorG = 1, double colorB = 0)
{
	GLUquadric* quadGLObj = gluNewQuadric();
	gluQuadricNormals(quadGLObj, GLU_SMOOTH);
	glColor3d(colorR,colorG,colorB);
	gluSphere(quadGLObj, r, ndiv, ndiv);
}

void drawTorus(int numc, int numt, double radius, double thickness = 0.01, dblArray3_t color = { 0.75, 0.75, 0.75 })
{
	int sign[2]{ 1, 0 };
	double s{ 0.0 };
	double t{ 0.0 };
	double x{ 0.0 };
	double y{ 0.0 };
	double z{ 0.0 };
	double twopi{ 2 * M_PI };
	for (int iii = 0; iii < numc; iii++)
	{
		glBegin(GL_QUAD_STRIP);
		glColor3d(color[0], color[1], color[2]);
		for (int jjj = 0; jjj <= numt; jjj++)
		{
			for (unsigned int kkk = 0; kkk < 2; kkk++)
			{
				s = (iii + sign[kkk]) % numc + 0.5;
				t = jjj % numt;
				x = (radius + thickness*cos(s * twopi / numc))*cos(t * twopi / numt);
				y = (radius + thickness*cos(s * twopi / numc))*sin(t * twopi / numt);
				z = thickness * sin(s * twopi / numc);
				glVertex3d(x, y, z);
			}//end kkk
		}//end jjj
		glEnd();
	}//end iii
	//glEnd(); //in python, this is nested in the top loop - doesn't make sense to call glEnd() multiple times...maybe a mistake
}			   //but just in case, if it torus doesn't draw right, move up a level

//private
void ParticlePic::incrementInd()
{
	//circular buffer
	//if (!maketrail_m)
		//return;
	/*if (trailind_m < (retain_m * interval_m))
	{
		trailind_m++;
		return;
	}*/
	//trailind_m = 0;
	//trail_m.begin();
}

void ParticlePic::setPos(const dblArray3_t& P)
{
	/*if (interval_m == 1) //code on bottom sets index to 0 every time if interval is 1
	{
		trail_m[trailind_m] = P;
		return;
	}
	if ((trailind_m % interval_m) == 0)
	{
		trail_m[(trailind_m / interval_m)] = P;
		//std::cout << "setPos if condition true " << trailind_m << " " << interval_m << " " << trailind_m % interval_m << std::endl;
	}*/
	//trail_m[2] = nulldA3_t;
	trail_m.pop_back();
	trail_m.push_front(P);
}

//public
void ParticlePic::draw(const dblArray3_t& P)
{
	setPos(P);
	//std::cout << trail_m.size() << std::endl;
	for (std::list<dblArray3_t>::iterator it = trail_m.begin(); it != trail_m.end(); ++it)
	{
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		//std::cout << trail_m[iii] << "  " << trailind_m;
		dblArray3_t trans = *it;
		glTranslated(trans[0], trans[1], trans[2]);
		drawSphere(radius_m);
		glPopMatrix();
	}
	//incrementInd();
}

void WireCoilPic::draw(const dblArray3_t& dummy)
{
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	dblArray3_t xyz = { 0.0, 0.0, -1.0 };
	dblArray3_t c = cross3dVectorsAndMultByConst(xyz, axis_m);
	glTranslated(cent_m[0], cent_m[1], cent_m[2]);
	glRotated(90, c[0], c[1], c[2]);
	drawTorus(5, 30, radius_m, thickness_m);
	glPopMatrix();
}