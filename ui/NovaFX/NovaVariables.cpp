///
/// NovaFX Variables implementation
///
/// Written by Vadim Kuras. 2010.
///

#include "NovaVariables.h"

#define SQR(x) (x*x)

namespace NovaFX
{
	NSphere::NSphere()
	{
		pos = new Nfloat4();	
	}

	NSphere::~NSphere()
	{
		if (pos)
			delete pos;
	}

	NCylinder::NCylinder()
	{
		pos = new Nfloat4();	
		left = new Nfloat3();
		up = new Nfloat3();
		forward = new Nfloat3();
	}

	NCylinder::~NCylinder()
	{
		if (pos)
			delete pos;

		if (left)
			delete left;

		if (up)
			delete up;

		if (forward)
			delete forward;
	}

	GLRGBA::GLRGBA(GLfloat red, GLfloat green, GLfloat blue, GLfloat alpha)
	{
		r = red;
		g = green;
		b = blue;
		a = alpha;
	}

	Nfloat3::Nfloat3(GLfloat x, GLfloat y, GLfloat z) : x(x), y(y), z(z)
	{}

	Nfloat3 Nfloat3::operator+(const Nfloat3 & other)
	{
		return Nfloat3(x + other.x, y + other.y, z + other.z);
	}

	Nfloat3 Nfloat3::operator-(const Nfloat3 & other)
	{
		return Nfloat3(x - other.x, y - other.y, z - other.z);
	}

	Nfloat3 Nfloat3::operator*(const GLfloat & other)
	{
		return Nfloat3(x * other, y * other, z * other);
	}
		
	GLfloat Nfloat3::operator*(const Nfloat3 & other)
	{
		return GLfloat(x * other.x + y * other.y + z * other.z);
	}

	GLfloat Nfloat3::Length(Nfloat3 * v)
	{
		return (GLfloat)(sqrt(SQR(v->x)+SQR(v->y)+SQR(v->z)));
	}

	Nfloat3 Nfloat3::CrossProduct(const Nfloat3 & other)
	{
		return Nfloat3(y*other.z - z*other.y, z*other.x - x*other.z,  x*other.y - y*other.x);
	}
			
	Nfloat3 Nfloat3::Normalize()
	{
		Nfloat3 res(0.0f, 0.0f, 0.0f);
		GLfloat len = Length(this);
		if (len == 0.0f) 
			return Nfloat3(0.0f, 0.0f, 0.0f);
		res.x = x / len;
		res.y = y / len;
		res.z = z / len;
		return res;
	}
}