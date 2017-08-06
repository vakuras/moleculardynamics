///
/// NovaCamera Header
///
/// Written by Vadim Kuras. 2010.
///

#ifndef NOVACAMERA_H
#define NOVACAMERA_H

#include "NovaVariables.h"

namespace NovaFX
{
	namespace NovaCamera
	{
		class Camera
		{
		private:
			Nfloat3			m_Position;
			Nfloat3			m_ViewDirection;
			Nfloat3			m_RightVector;
			Nfloat3			m_UpVector;
			Nfloat3			m_Rotated;
			GLfloat				m_ForwardSpeed;
			GLfloat				m_StrafeSpeed;
			GLfloat				m_UpwardSpeed;
			GLfloat				m_UpwardDevide;
			GLfloat				m_StrafeDevide;
			GLfloat				m_ForwardDevide;
		public:
			Camera(GLfloat forwardDevide, GLfloat strafeDevide, GLfloat upwardDevide);
			void Move(Nfloat3 direction);
			void RotateX(GLfloat angle);
			void RotateY(GLfloat angle);
			void RotateZ(GLfloat angle);
			void Render();
			void MoveForward(GLfloat speed);
			void StrafeRight(GLfloat speed);
			void MoveUpward(GLfloat speed);
		};
	}
}

#endif