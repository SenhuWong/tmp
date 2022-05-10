#include "vector3d.h"
#include <cmath>

Vector3D::Vector3D()
{
}

Vector3D::Vector3D(double *vin)
{
	for (int i = 0; i < d_dim; i++)
	{
		components[i] = vin[i];
	}
}

Vector3D::Vector3D(double c1, double c2, double c3)
{
	components[0] = c1;
	components[1] = c2;
	components[2] = c3;
}
Vector3D::Vector3D(const Vector3D &another)
{
	(*this)[0] = another[0];
	(*this)[1] = another[1];
	(*this)[2] = another[2];
}
double &Vector3D::operator[](int i)
{
	return components[i];
}
const double &Vector3D::operator[](int i) const
{
	return components[i];
}
Vector3D Vector3D::operator-(const Vector3D &another) const
{
	Vector3D result = Vector3D(components[0] - another.components[0],
							   components[1] - another.components[1],
							   components[2] - another.components[2]);
	return result;
}
double Vector3D::dot_product(Vector3D &another) const
{
	double result = 0;
	for (int i = 0; i < d_dim; i++)
	{
		result += this->components[i] * another.components[i];
	}
	return result;
}
double Vector3D::dot_product(const Vector3D &another) const
{
	double result = 0;
	for (int i = 0; i < d_dim; i++)
	{
		result += this->components[i] * another.components[i];
	}
	return result;
}
Vector3D &Vector3D::operator=(const Vector3D &another)
{
	(*this)[0] = another[0];
	(*this)[1] = another[1];
	(*this)[2] = another[2];
	return *this;
}

Vector3D &Vector3D::normalized()
{
	double norm = 0;
	for (int i = 0; i < 3; i++)
	{
		norm += std::pow(components[i], 2);
	}
	norm = sqrt(norm);
	for (int i = 0; i < 3; i++)
	{
		components[i] = components[i] / norm;
	}
	return *this;
}

Vector3D Vector3D::operator-(const Vector3D &another)
{
	Vector3D result = Vector3D((*this)[0] - another[0], (*this)[1] - another[1], (*this)[2] - another[2]);
	return result;
}

Vector3D Vector3D::cross_product(const Vector3D &another) const
{
	Vector3D result;
	result[0] = (*this)[1] * another[2] - (*this)[2] * another[1];
	result[1] = (*this)[2] * another[0] - (*this)[0] * another[2];
	result[2] = (*this)[0] * another[1] - (*this)[1] * another[0];
	return result;
}

double Vector3D::BoxProduct(const Vector3D &a, const Vector3D &b)
{
	double result = (a[0] * b[1] - a[1] * b[0]) * this->components[2] + (a[1] * b[2] - a[2] * b[1]) * this->components[0] + (a[2] * b[0] - a[0] * b[2]) * this->components[1];
	return result;
}


//2D
Vector2D::Vector2D()
{
}

Vector2D::Vector2D(double *vin)
{
	for (int i = 0; i < d_dim; i++)
	{
		components[i] = vin[i];
	}
}

Vector2D::Vector2D(double c1, double c2)
{
	components[0] = c1;
	components[1] = c2;
}
Vector2D::Vector2D(const Vector2D &another)
{
	(*this)[0] = another[0];
	(*this)[1] = another[1];
}
double &Vector2D::operator[](int i)
{
	return components[i];
}
const double &Vector2D::operator[](int i) const
{
	return components[i];
}
Vector2D Vector2D::operator-(const Vector2D &another) const
{
	Vector2D result = Vector2D(components[0] - another.components[0],
							   components[1] - another.components[1]);
	return result;
}
double Vector2D::dot_product(Vector2D &another) const
{
	double result = 0;
	for (int i = 0; i < d_dim; i++)
	{
		result += this->components[i] * another.components[i];
	}
	return result;
}
double Vector2D::dot_product(const Vector2D &another) const
{
	double result = 0;
	for (int i = 0; i < d_dim; i++)
	{
		result += this->components[i] * another.components[i];
	}
	return result;
}
Vector2D &Vector2D::operator=(const Vector2D &another)
{
	(*this)[0] = another[0];
	(*this)[1] = another[1];
	return *this;
}

Vector2D &Vector2D::normalized()
{
	double norm = 0;
	for (int i = 0; i < 2; i++)
	{
		norm += std::pow(components[i], 2);
	}
	norm = sqrt(norm);
	for (int i = 0; i < 2; i++)
	{
		components[i] = components[i] / norm;
	}
	return *this;
}

Vector2D Vector2D::operator-(const Vector2D &another)
{
	Vector2D result = Vector2D((*this)[0] - another[0], (*this)[1] - another[1]);
	return result;
}
