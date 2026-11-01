#pragma once

//std
#include <cstdint>
#include <initializer_list>

namespace math
{
	namespace linear
	{
		class Vec2;
	}
}

namespace math
{
	namespace linear
	{
		class Mat2
		{
		public:
			//constructors
			Mat2(void);
			Mat2(const Mat2&) = default;
			explicit Mat2(const double*);
			explicit Mat2(const Vec2&, const Vec2&);
			Mat2(std::initializer_list<std::initializer_list<double>>);

			//destructor
			~Mat2(void) = default;

			//data
			double* data(void);
			const double* data(void) const;

			//operators
			Mat2 operator+(void) const;
			Mat2 operator-(void) const;

			Mat2 operator/(double) const;
			Mat2 operator+(const Mat2&) const;
			Mat2 operator-(const Mat2&) const;
			Mat2 operator*(const Mat2&) const;
			Vec2 operator*(const Vec2&) const;

			Mat2& operator+=(double);
			Mat2& operator-=(double);
			Mat2& operator*=(double);
			Mat2& operator/=(double);

			Mat2& operator+=(const Mat2&);
			Mat2& operator-=(const Mat2&);

			double& operator[](uint32_t);
			double& operator()(uint32_t, uint32_t);

			const double& operator[](uint32_t) const;
			const double& operator()(uint32_t, uint32_t) const;

			//linear
			static Mat2 eye(void);
			Mat2 inverse(void) const;
			Mat2 transpose(void) const;
			Vec2 eigenvalues(void) const;

			double trace(void) const;
			double determinant(void) const;
			double invariant_I1(void) const;
			double invariant_I2(void) const;
			double deviatoric_J2(void) const;

			//friends
			friend Mat2 operator*(double, const Mat2&);

		private:
			//data
			double m_data[4];

			//friends
			friend class math::linear::Vec2;
		};
	}
}