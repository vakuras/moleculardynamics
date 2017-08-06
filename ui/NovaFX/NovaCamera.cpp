///
/// NovaCamera implementation
///
/// Written by Vadim Kuras. 2010.
///

#include "NovaCamera.h"

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#define PIdiv180 (float)M_PI/180.0f

namespace NovaFX
{
	namespace NovaCamera
	{
		Camera::Camera(GLfloat forwardDevide, GLfloat strafeDevide, GLfloat upwardDevide) :
				m_ForwardDevide(forwardDevide),
				m_StrafeDevide(strafeDevide),
				m_UpwardDevide(upwardDevide)
		{
			m_Position = Nfloat3(0.0f, 0.0f, 0.0f);
			m_ViewDirection = Nfloat3(0.0f, 0.0f, -1.0f);
			m_RightVector = Nfloat3(1.0f, 0.0f, 0.0f);
			m_UpVector = Nfloat3(0.0f, 1.0f, 0.0f);
			m_Rotated = Nfloat3(0.0f, 0.0f, 0.0f);

			m_UpwardSpeed = m_StrafeSpeed = m_ForwardSpeed = 0;
		}

		void Camera::Move(Nfloat3 direction)
		{
			m_Position = m_Position + direction;
		}

		void Camera::RotateX(GLfloat angle)
		{
			m_Rotated.x += angle;
	
			m_ViewDirection = (m_ViewDirection*cosf(angle*PIdiv180) 
				+ m_UpVector*sinf(angle*PIdiv180)).Normalize();

			m_UpVector = m_ViewDirection.CrossProduct(m_RightVector)*-1;
		}

		void Camera::RotateY(GLfloat angle)
		{
			m_Rotated.y += angle;
	
			m_ViewDirection = (m_ViewDirection*cosf(angle*PIdiv180)
				- m_RightVector*sinf(angle*PIdiv180)).Normalize();

			m_RightVector = m_ViewDirection.CrossProduct(m_UpVector);
		}

		void Camera::RotateZ(GLfloat angle)
		{
			m_Rotated.z += angle;
	
			m_RightVector = (m_RightVector*cosf(angle*PIdiv180)
				+ m_UpVector*sinf(angle*PIdiv180)).Normalize();

			m_UpVector = m_ViewDirection.CrossProduct(m_RightVector)*-1;
		}

		void Camera::Render()
		{
			m_Position = m_Position + (m_ViewDirection*-m_ForwardSpeed);
			m_Position = m_Position + (m_RightVector*m_StrafeSpeed);
			m_Position = m_Position + (m_UpVector*m_UpwardSpeed);

			m_UpwardSpeed = m_UpwardSpeed / m_UpwardDevide;
			m_StrafeSpeed = m_StrafeSpeed / m_StrafeDevide;
			m_ForwardSpeed = m_ForwardSpeed / m_ForwardDevide;

			Nfloat3 viewPoint = m_Position + m_ViewDirection;

			gluLookAt(m_Position.x, m_Position.y, m_Position.z,
						viewPoint.x, viewPoint.y, viewPoint.z,
						m_UpVector.x, m_UpVector.y, m_UpVector.z);
		}

		void Camera::MoveForward(GLfloat speed)
		{
			m_ForwardSpeed = speed;
		}

		void Camera::StrafeRight(GLfloat speed)
		{
			m_StrafeSpeed = speed;
		}

		void Camera::MoveUpward(GLfloat speed)
		{
			m_UpwardSpeed = speed;
		}

	}
}