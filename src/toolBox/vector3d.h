#pragma once
#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
namespace GeomElements
{
	//A vector3d is a vector living in a 3d world.
	template<int ndim, typename T>
	class vector3d
	{
	private:
		T components[ndim];
	public:
		vector3d();
		vector3d(const T& c0, const T& c1, const T& c2 = 0);
		vector3d(const vector3d<ndim, T>& another);

		T& operator[](int ind);
		vector3d<ndim, T> operator+(const vector3d<ndim, T>& another) const;
		vector3d<ndim, T>& operator=(const vector3d<ndim, T>& another);
		vector3d<ndim, T> operator-(const vector3d<ndim, T>& another) const;
		vector3d<ndim, T> operator/(const T& divisor) const;
		vector3d<ndim, T> operator*(const T& multiple) const;
		T dot_product(const vector3d<ndim, T>& another);
		std::string to_string() const
		{
			std::string vec_string = "vector" + std::to_string(ndim) + "d :(";
			for (int i = 0; i < ndim; i++)
			{
				vec_string += std::to_string(components[i]) + " ";
			}
			vec_string += ")\n";
			return vec_string;
		}
		vector3d<ndim, T> cross_product(const vector3d<ndim, T>& another) const;
		T normalize();
		T L2Square()
		{
			double VMSquare = 0;
			for(int i =0;i<ndim;i++)
			{
				VMSquare+= components[i]*components[i];
			}
			return VMSquare;
		}
		void reset()
		{
			for(int i =0;i<ndim;i++)
			{
				components[i] = 0;
			}
		}

		void show();
	};



	template<int ndim, typename T>
	T vector3d<ndim, T>::dot_product(const vector3d<ndim, T>& another)
	{
		T result = 0;
		for (int i = 0; i < ndim; i++)
		{
			result += this->components[i] * another.components[i];
		}
		return result;
	}

	template<int ndim, typename T>
	vector3d<ndim, T> vector3d<ndim, T>::cross_product(const vector3d<ndim, T>& another) const
	{
		if (ndim == 2)
		{
			throw std::runtime_error("2D vector making a cross_product\n");
		}
		T c0 = this->components[1] * another.components[2] - this->components[2] * another.components[1];
		T c1 = this->components[2] * another.components[0] - this->components[0] * another.components[2];
		T c2 = this->components[0] * another.components[1] - this->components[1] * another.components[0];

		return vector3d<ndim, T>(c0, c1, c2);
	}
	
	template<int ndim, typename T>
	T& vector3d<ndim, T>::operator[](int ind)
	{
		if (ind < 0 or ind >= ndim)
		{
			throw std::runtime_error("Vector3d index out of range");
		}
		return components[ind];

	}
	
	template<int ndim, typename T>
	vector3d<ndim, T>& vector3d<ndim, T>::operator=(const vector3d<ndim, T>& another)
	{
		for (int i = 0; i < ndim; i++)
		{
			this->components[i] = another.components[i];
		}
		return *this;
	}
	
	template<int ndim, typename T>
	vector3d<ndim, T> vector3d<ndim, T>::operator+(const vector3d<ndim, T>& another) const
	{
		vector3d<ndim, T> result;
		for (int i = 0; i < ndim; i++)
		{
			result[i] = this->components[i] + another.components[i];
		}
		return result;
	}
	
	template<int ndim, typename T>
	vector3d<ndim, T> vector3d<ndim, T>::operator-(const vector3d<ndim, T>& another) const
	{
		vector3d<ndim, T> result;
		for (int i = 0; i < ndim; i++)
		{
			result[i] = this->components[i] - another.components[i];
		}
		return result;
	}
	
	template<int ndim, typename T>
	vector3d<ndim, T> vector3d<ndim, T>::operator/(const T& divisor) const
	{
		vector3d<ndim,T> result;
		for (int i = 0; i < ndim; i++)
		{
			result.components[i] = this->components[i] / divisor;
		}
		return result;
	}
	
	template<int ndim, typename T>
	vector3d<ndim, T> vector3d<ndim, T>::operator*(const T& multiple) const
	{
		vector3d<ndim,T> result;
		for (int i = 0; i < ndim; i++)
		{
			result.components[i] = this->components[i] *multiple;
		}
		return result;
	}
	
	template<int ndim, typename T>
	T vector3d<ndim, T>::normalize()
	{
		double L2 = 0;
		for (int i = 0; i < ndim; i++)
		{
			L2 += (components[i] * components[i]);
		}
		L2 = sqrt(L2);
		for (int i = 0; i < ndim; i++)
		{
			components[i] = components[i] / L2;
		}
		return L2;
	}

	template<int ndim, typename T>
	vector3d<ndim, T>::vector3d()
	{
		for(int i = 0;i<ndim;i++)
		{
			components[i] = 0;
		}
	}
	
	template<int ndim, typename T>
	vector3d<ndim, T>::vector3d(const T& c0, const T& c1, const T& c2)
	{
		if (ndim == 2)
		{
			components[0] = c0;
			components[1] = c1;
		}
		else if (ndim == 3)
		{
			components[0] = c0;
			components[1] = c1;
			components[2] = c2;
		}
	}
	
	template<int ndim, typename T>
	vector3d<ndim, T>::vector3d(const vector3d<ndim, T>& another)
	{
		for (int i = 0; i < ndim; i++)
		{
			components[i] = another.components[i];
		}
	}

	template<int ndim, typename T>
	void vector3d<ndim, T>::show()
	{
		for (int i = 0; i < ndim; i++)
		{
			std::cout << components[i] << '\t';
		}
		std::cout << '\n';
	}

};
