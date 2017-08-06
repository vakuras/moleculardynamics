///
/// NovaFX GL Header
///
/// Written by Vadim Kuras. 2010.
///

#ifndef NOVAFX_H
#define NOVAFX_H

#include <vector>
#include "NovaCamera.h"

using namespace std;

namespace NovaFX
{
	class NovaFX
	{
	private:
		HWND						m_hWnd;
		int							m_Bits;
		HDC							m_hDC;
		HGLRC						m_hRC;
		UINT						m_Height;
		UINT						m_Width;
		NovaCamera::Camera*			m_Camera;
		GLUquadricObj*				m_Quadric;
		vector<NSphere*>			m_Spheres;
		vector<NCylinder*>			m_Cylinders;
		bool						m_DrawSpheres;
		bool						m_DrawCylinders;

	public:
		void HideSpheres();
		void ShowSpheres();
		void HideCylinders();
		void ShowCylinders();
		NovaFX(HWND hWnd, int width, int height, int bits);
		~NovaFX();
		bool Initialize(GLfloat r, GLfloat g, GLfloat b, GLfloat a);
		bool Render();
		void Resize(int newWidth, int newHeight);
		size_t AddSphere(NSphere * sphere);
		NSphere * GetSphere(size_t index);
		size_t AddCylinder(NCylinder * sphere);
		NCylinder * GetCylinder(size_t index);
		void RemoveSphere(size_t index);
		void RemoveCylinder(size_t index);
		void RemoveAllSpheres();
		void RemoveAllCylinders();
		void SetCamera(NovaCamera::Camera * camera);
		NovaCamera::Camera * GetCamera();
	protected:
		bool Error(LPCTSTR msg);
	};
}

#endif