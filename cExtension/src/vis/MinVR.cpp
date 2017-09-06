#include "vis/MinVR.h"

using namespace MinVR;

	// Callback for event handling, inherited from VREventHandler
void MyVRApp::onVREvent(const std::string &eventName, VRDataIndex *eventData) {
	//std::cout << "Event: " << eventName << std::endl;
	if (eventName == "/KbdEsc_Down") {
		_quit = true;
	}
	else if (eventName == "/MouseBtnLeft_Down") {
		_radius += 5.0 * _incAngle;
	}
	else if (eventName == "/MouseBtnRight_Down") {
		_radius -= 5.0 * _incAngle;
	}
	else if ((eventName == "/KbdLeft_Down") || (eventName == "/KbdLeft_Repeat")) {
		_horizAngle -= _incAngle;
	}
	else if ((eventName == "/KbdRight_Down") || (eventName == "/KbdRight_Repeat")) {
		_horizAngle += _incAngle;
	}
	else if ((eventName == "/KbdUp_Down") || (eventName == "/KbdUp_Repeat")) {
		_vertAngle -= _incAngle;
	}
	else if ((eventName == "/KbdDown_Down") || (eventName == "/KbdDown_Repeat")) {
		_vertAngle += _incAngle;
	}

	if (_horizAngle > 6.283185) _horizAngle -= 6.283185;
	if (_horizAngle < 0.0) _horizAngle += 6.283185;

	if (_vertAngle > 6.283185) _vertAngle -= 6.283185;
	if (_vertAngle < 0.0) _vertAngle += 6.283185;
}


void MyVRApp::onVRRenderContext(VRDataIndex *renderState, VRDisplayNode *callingNode) {
	if (!renderState->exists("IsConsole", "/")) {
	}
}

// Callback for rendering, inherited from VRRenderHandler
void MyVRApp::onVRRenderScene(VRDataIndex *renderState, VRDisplayNode *callingNode) {
	if (renderState->exists("IsConsole", "/")) {
		VRConsoleNode *console = dynamic_cast<VRConsoleNode*>(callingNode);
		console->println("Console output...");
	}
	else {
		dblArray3_t elecP = partBVP(partList_m[0]);
		glClear(GL_COLOR_BUFFER_BIT);
		int width = renderState->getValue("WindowWidth", "/");
		int height = renderState->getValue("WindowHeight", "/");
		double ratio = width / height;
		glViewport(0, 0, width, height);
		glClear(GL_COLOR_BUFFER_BIT);

		if (renderState->exists("ProjectionMatrix", "/")) {
			// This is the typical case where the MinVR DisplayGraph contains
			// an OffAxisProjectionNode or similar node, which sets the
			// ProjectionMatrix and ViewMatrix based on head tracking data.

			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			VRMatrix4 P = renderState->getValue("ProjectionMatrix", "/");
			glLoadMatrixd(P.m);

			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();

			// In the original adaptee.cpp program, keyboard and mouse commands
			// were used to adjust the camera.  Now that we are creating a VR
			// program, we want the camera (Projection and View matrices) to
			// be controlled by head tracking.  So, we switch to having the
			// keyboard and mouse control the Model matrix.  Guideline: In VR,
			// put all of the "scene navigation" into the Model matrix and
			// leave the Projection and View matrices for head tracking.
			VRMatrix4 M = VRMatrix4::translation(VRVector3(0.0, 0.0, -_radius)) *
				VRMatrix4::rotationX(_vertAngle) *
				VRMatrix4::rotationY(_horizAngle);

			VRMatrix4 V = renderState->getValue("ViewMatrix", "/");
			glLoadMatrixd((V*M).m);
			glTranslatef(0.0, 0.0, -2);
		}
		else {
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(60, width / height, 0.1, 50.0);
			glMatrixMode(GL_MODELVIEW);
			glLoadIdentity();
			glTranslatef(0.0, 0.0, -15);
		}
		
		drawPics(elecP);
		for (unsigned int iii = 0; iii < 3; iii++)
		{
			if ((elecP[iii] < -5) || (elecP[iii] > 5))
			{
				dblArray3_t updP{ -2.0, 0.0, 0.0 };
				dblArray3_t updV{ 1000.0, 1000.0, 1000.0 };
				partList_m[0]->updPmanual(updP);
				partList_m[0]->updPmanual(updV);
			}
		}
	}
}

void MyVRApp::run() {
	while (!_quit) {
		_vrMain->mainloop();
	}
}

void MyVRApp::drawPics(const dblArray3_t& P)
{
	for (unsigned int iii = 0; iii < picList_m.size(); iii++)
	{
		if (!(picList_m[iii]))
			continue;
		picList_m[iii]->draw(P);
	}
}