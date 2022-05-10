#pragma once

class Vector3D
{
private:
	int d_dim = 3;
	double components[3] = { 0,0,0 };//3d at most;
public:
	//constructors
	Vector3D();
	Vector3D(double* vin);
	Vector3D(double c1, double c2, double c3);
	Vector3D(const Vector3D& another);
	Vector3D& normalized();
	//operators
	double& operator[](int i);
	const double& operator[](int i)const;
	Vector3D operator-(const Vector3D& another) const;
	Vector3D& operator = (const Vector3D& another);
	//operation
	double dot_product(const Vector3D& another) const;
	double dot_product(Vector3D& another) const;
	Vector3D cross_product(const Vector3D& another) const;
	double BoxProduct(const Vector3D& a, const Vector3D& b);
	Vector3D operator - (const Vector3D&another);
};

class Vector2D
{
private:
	int d_dim = 2;
	double components[3] = { 0,0};//2d at most;
public:
	//constructors
	Vector2D();
	Vector2D(double* vin);
	Vector2D(double c1, double c2);
	Vector2D(const Vector2D& another);
	Vector2D& normalized();
	//operators
	double& operator[](int i);
	const double& operator[](int i)const;
	Vector2D operator-(const Vector2D& another) const;
	Vector2D& operator = (const Vector2D& another);
	//operation
	double dot_product(const Vector2D& another) const;
	double dot_product(Vector2D& another) const;
	//Vector2D cross_product(const Vector2D& another) const;
	//double BoxProduct(const Vector2D& a, const Vector2D& b);
	Vector2D operator - (const Vector2D&another);
};

