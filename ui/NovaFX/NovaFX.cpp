///
/// NovaFX GL
///
/// Written by Vadim Kuras. 2010.
///

#include "NovaFX.h"

namespace NovaFX
{
	NovaFX::NovaFX(HWND hWnd, int width, int height, int bits) : m_hWnd(hWnd), 
																 m_Width(width),
																 m_Height(height),
																 m_Bits(bits)
	{
		m_Quadric = NULL;
		m_DrawSpheres = true;
		m_DrawCylinders = true;
	}

	NovaFX::~NovaFX()
	{
		gluDeleteQuadric(m_Quadric);
		m_Quadric = NULL;

		for(vector<NSphere*>::iterator it = m_Spheres.begin(); it < m_Spheres.end(); it++)
			delete (*it);

		for(vector<NCylinder*>::iterator it = m_Cylinders.begin(); it < m_Cylinders.end(); it++)
			delete (*it);
	}

	void NovaFX::HideSpheres()
	{
		m_DrawSpheres = false;
	}

	void NovaFX::ShowSpheres()
	{
		m_DrawSpheres = true;
	}

	void NovaFX::HideCylinders()
	{
		m_DrawCylinders = false;
	}

	void NovaFX::ShowCylinders()
	{
		m_DrawCylinders = true;
	}

	bool NovaFX::Initialize(GLfloat r, GLfloat g, GLfloat b, GLfloat a)
	{
		GLuint pixelFormat;

		PIXELFORMATDESCRIPTOR pfd=
		{
			sizeof(PIXELFORMATDESCRIPTOR),
			1,
			PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,
			PFD_TYPE_RGBA,
			m_Bits,
			0, 0, 0, 0, 0, 0,
			0,
			0,
			0,
			0, 0, 0, 0,
			16,
			0,
			0,
			PFD_MAIN_PLANE,
			0,
			0, 0, 0
		};

		if (!(m_hDC=GetDC(m_hWnd)))
			return Error(TEXT("Can't Create A GL Device Context."));

		if (!(pixelFormat=ChoosePixelFormat(m_hDC,&pfd)))
			return Error(TEXT("Can't Find A Suitable PixelFormat."));

		if(!SetPixelFormat(m_hDC,pixelFormat,&pfd))
			return Error(TEXT("Can't Set The PixelFormat."));

		if (!(m_hRC=wglCreateContext(m_hDC)))
			return Error(TEXT("Can't Create A GL Rendering Context."));

		if(!wglMakeCurrent(m_hDC,m_hRC))
			return Error(TEXT("Can't Activate The GL Rendering Context."));

		glShadeModel(GL_SMOOTH);
		glClearColor(r, g, b, a);
		glClearDepth(1.0f);
		glEnable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
		glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
		glEnable(GL_CULL_FACE);
		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);

		Resize(m_Width, m_Height);

		return true;
	}

	void NovaFX::Resize(int newWidth, int newHeight)
	{
		m_Width = newWidth;
		m_Height = newHeight;

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		gluPerspective(45.0f,(GLfloat)m_Width/(GLfloat)m_Height,0.1f,2000.0f);

		glMatrixMode(GL_MODELVIEW);
		glViewport(0,0,m_Width,m_Height);
	}
	
	bool NovaFX::Render()
	{
		if (!m_Camera)
			return false;

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glLoadIdentity();

		m_Camera->Render();

		if (m_DrawSpheres)
		{
			for(vector<NSphere*>::iterator it = m_Spheres.begin(); it < m_Spheres.end(); it++)
			{
				glColor4f((*it)->color.r,(*it)->color.g,(*it)->color.b,(*it)->color.a);
				glPushMatrix();	
				glTranslatef((*it)->pos->x,(*it)->pos->y,(*it)->pos->z);
				gluQuadricNormals(m_Quadric, GLU_SMOOTH);
				gluSphere(m_Quadric, (*it)->pos->w, (*it)->slices, (*it)->stacks);
				glPopMatrix();
			}
		}

		if (m_DrawCylinders)
		{
			for(vector<NCylinder*>::iterator it = m_Cylinders.begin(); it < m_Cylinders.end(); it++)
			{
				glColor4f((*it)->color.r,(*it)->color.g,(*it)->color.b,(*it)->color.a);
				glPushMatrix();	
				glTranslatef((*it)->pos->x,(*it)->pos->y,(*it)->pos->z);
				glRotatef(1.0f,(*it)->left->x,(*it)->left->y,(*it)->left->z);
				glRotatef(1.0f,(*it)->up->x,(*it)->up->y,(*it)->up->z);
				glRotatef(1.0f,(*it)->forward->x,(*it)->forward->y,(*it)->forward->z);
				gluQuadricNormals(m_Quadric, GLU_SMOOTH);
				gluCylinder(m_Quadric, (*it)->baseRadius, (*it)->topRadius ,(*it)->pos->w, (*it)->slices, (*it)->stacks);
				glPopMatrix();
			}
		}

		glFlush();
		SwapBuffers(m_hDC);

		return true;		
	}

	bool NovaFX::Error(LPCTSTR msg)
	{
		MessageBox(NULL, msg, TEXT("Error"), MB_OK | MB_ICONEXCLAMATION);
		return false;
	}

	size_t NovaFX::AddSphere(NSphere * sphere)
	{
		if (!m_Quadric)
			m_Quadric = gluNewQuadric();

		m_Spheres.push_back(sphere);
		return m_Spheres.size()-1;
	}

	size_t NovaFX::AddCylinder(NCylinder * cylinder)
	{
		if (!m_Quadric)
			m_Quadric = gluNewQuadric();

		m_Cylinders.push_back(cylinder);
		return m_Cylinders.size()-1;
	}

	NSphere * NovaFX::GetSphere(size_t index)
	{
		if (m_Spheres.empty() || m_Spheres.size() < index)
			return NULL;

		return m_Spheres[index];
	}

	NCylinder * NovaFX::GetCylinder(size_t index)
	{
		if (m_Cylinders.empty() || m_Cylinders.size() < index)
			return NULL;

		return m_Cylinders[index];
	}

	void NovaFX::RemoveSphere(size_t index)
	{
		if (m_Spheres.empty() || m_Spheres.size() < index)
			return;

		delete m_Spheres[index];

		m_Spheres.erase(m_Spheres.begin() + index);

		if (m_Cylinders.empty() && m_Spheres.empty() && m_Quadric)
		{
			gluDeleteQuadric(m_Quadric);
			m_Quadric = NULL;
		}
	}

	void NovaFX::RemoveCylinder(size_t index)
	{
		if (m_Cylinders.empty() || m_Cylinders.size() < index)
			return;

		delete m_Cylinders[index];

		m_Cylinders.erase(m_Cylinders.begin() + index);

		if (m_Cylinders.empty() && m_Spheres.empty() && m_Quadric)
		{
			gluDeleteQuadric(m_Quadric);
			m_Quadric = NULL;
		}
	}

	void NovaFX::RemoveAllSpheres()
	{
		if (m_Spheres.empty())
			return;

		for(vector<NSphere*>::iterator it = m_Spheres.begin(); it<m_Spheres.end(); it++)
			delete *it;

		m_Spheres.clear();

		if (m_Quadric && m_Cylinders.empty())
		{
			gluDeleteQuadric(m_Quadric);
			m_Quadric = NULL;
		}
	}

	void NovaFX::RemoveAllCylinders()
	{
		if (m_Cylinders.empty())
			return;

		for(vector<NCylinder*>::iterator it = m_Cylinders.begin(); it<m_Cylinders.end(); it++)
			delete *it;

		m_Cylinders.clear();

		if (m_Quadric && m_Spheres.empty())
		{
			gluDeleteQuadric(m_Quadric);
			m_Quadric = NULL;
		}
	}

	void NovaFX::SetCamera(NovaCamera::Camera * camera)
	{
		m_Camera = camera;
	}

	NovaCamera::Camera * NovaFX::GetCamera()
	{
		return m_Camera;
	}
}