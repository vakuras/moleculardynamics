///
/// NovaFX Variables Header
///
/// Written by Vadim Kuras. 2010.
///

#ifndef NOVAVARS_H
#define NOVAVARS_H

#include <Windows.h>
#include <gl\GL.h>
#include <gl\GLU.h>
#include <cmath>

namespace NovaFX
{
	struct Nfloat3
	{
		GLfloat x;
		GLfloat y;
		GLfloat z;

		Nfloat3() {}
		Nfloat3(GLfloat x, GLfloat y, GLfloat z);
		Nfloat3 operator+(const Nfloat3 & other);
		Nfloat3 operator-(const Nfloat3 & other);
		Nfloat3 operator*(const GLfloat & other);
		GLfloat operator*(const Nfloat3 & other);
		GLfloat Length(Nfloat3 * v);
		Nfloat3 CrossProduct(const Nfloat3 & other);
		Nfloat3 Normalize();
	};

	struct Nfloat4
	{
		GLfloat x;
		GLfloat y;
		GLfloat z;
		GLfloat w;
	};

	struct GLRGBA
	{
		GLfloat r;
		GLfloat g;
		GLfloat b;
		GLfloat a;
		
		GLRGBA(){}
		GLRGBA(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha);
	};

	struct NSphere
	{
		Nfloat4 * pos;
		GLuint slices;
		GLuint stacks;
		GLRGBA color;

		NSphere();
		~NSphere();
	};

	struct NCylinder
	{
		Nfloat4 * pos;
		Nfloat3 * left;
		Nfloat3 * up;
		Nfloat3 * forward;
		GLuint slices;
		GLuint stacks;
		GLfloat baseRadius;
		GLfloat topRadius;
		GLRGBA color;

		NCylinder();
		~NCylinder();
	};
}

#endif