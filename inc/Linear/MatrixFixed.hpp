#pragma once

//std
#include <cstdint>
#include <cstring>

namespace math
{
	template<uint32_t r, uint32_t c> class MatrixFixed
	{
	public:
		//constructor
		MatrixFixed(void);
		MatrixFixed(const MatrixFixed&);

		explicit MatrixFixed(double*);
		explicit MatrixFixed(const double*);

		//destructor
		~MatrixFixed(void) = default;

		//data
		double* data(void);
		const double* data(void) const;

		//operators
		MatrixFixed operator+(void) const;
		MatrixFixed operator-(void) const;
		
		MatrixFixed& operator+=(double);
		MatrixFixed& operator-=(double);
		MatrixFixed& operator*=(double);
		MatrixFixed& operator/=(double);
		MatrixFixed& operator+=(const MatrixFixed&);
		MatrixFixed& operator-=(const MatrixFixed&);
		
		MatrixFixed& operator=(const MatrixFixed&);
		MatrixFixed operator+(const MatrixFixed&) const;
		MatrixFixed operator-(const MatrixFixed&) const;

		double& operator[](uint32_t);
		double& operator()(uint32_t, uint32_t);
		const double& operator[](uint32_t) const;
		const double& operator()(uint32_t, uint32_t) const;

		template<uint32_t d> MatrixFixed<r, d> operator*(const MatrixFixed<c, d>&) const;

		//linear
		MatrixFixed<c, r> transpose(void) const;

		//friends
		template<uint32_t a, uint32_t b> friend MatrixFixed<a, b> operator*(double, const MatrixFixed<a, b>&);

	protected:
		//data
		double* m_data_ptr;
		const double* m_data_ref;
		double m_data_mem[r * c];

		//friends
		template<uint32_t, uint32_t> friend class MatrixFixed;
	};

	//constructor
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>::MatrixFixed(void) : m_data_ptr{m_data_mem}, m_data_ref{m_data_mem}
	{
		return;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>::MatrixFixed(const MatrixFixed& M) : m_data_ptr{m_data_mem}, m_data_ref{m_data_mem}
	{
		memcpy(m_data_mem, M.m_data_ref, r * c * sizeof(double));
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>::MatrixFixed(double* data_ptr) : m_data_ptr{data_ptr}, m_data_ref{data_ptr}
	{
		return;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>::MatrixFixed(const double* data_ref) : m_data_ptr{nullptr}, m_data_ref{data_ref}
	{
		return;
	}

	//data
	template<uint32_t r, uint32_t c>  double* MatrixFixed<r, c>::data(void)
	{
		return m_data_ptr;
	}
	template<uint32_t r, uint32_t c> const double* MatrixFixed<r, c>::data(void) const
	{
		return m_data_ref;
	}

	//operators
	template<uint32_t r, uint32_t c> MatrixFixed<r, c> MatrixFixed<r, c>::operator+(void) const
	{
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c> MatrixFixed<r, c>::operator-(void) const
	{
		MatrixFixed R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data_ptr[i] = -m_data_ref[i];
		return R;
	}

	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator+=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data_ptr[i] += s;
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator-=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data_ptr[i] -= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator*=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data_ptr[i] *= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator/=(double s)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data_ptr[i] /= s;
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator+=(const MatrixFixed& M)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data_ptr[i] += M.m_data_ref[i];
		return *this;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c>& MatrixFixed<r, c>::operator-=(const MatrixFixed& M)
	{
		for(uint32_t i = 0; i < r * c; i++) m_data_ptr[i] -= M.m_data_ref[i];
		return *this;
	}

	template<uint32_t r, uint32_t c>  MatrixFixed<r, c>& MatrixFixed<r, c>::operator=(const MatrixFixed& M)
	{
		memcpy(m_data_ptr, M.m_data_ref, r * c * sizeof(double));
		return *this;
	}
	template<uint32_t r, uint32_t c>  MatrixFixed<r, c> MatrixFixed<r, c>::operator+(const MatrixFixed& M) const
	{
		MatrixFixed R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data_ptr[i] = m_data_ref[i] + M.m_data_ref[i];
		return R;
	}
	template<uint32_t r, uint32_t c> MatrixFixed<r, c> MatrixFixed<r, c>::operator-(const MatrixFixed& M) const
	{
		MatrixFixed R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data_ptr[i] = m_data_ref[i] - M.m_data_ref[i];
		return R;
	}

	template<uint32_t r, uint32_t c> double& MatrixFixed<r, c>::operator[](uint32_t index)
	{
		return m_data_ptr[index];
	}
	template<uint32_t r, uint32_t c> double& MatrixFixed<r, c>::operator()(uint32_t index_1, uint32_t index_2)
	{
		return m_data_ptr[index_1 + r * index_2];
	}
	template<uint32_t r, uint32_t c> const double& MatrixFixed<r, c>::operator[](uint32_t index) const
	{
		return m_data_ref[index];
	}
	template<uint32_t r, uint32_t c> const double& MatrixFixed<r, c>::operator()(uint32_t index_1, uint32_t index_2) const
	{
		return m_data_ref[index_1 + r * index_2];
	}

	template<uint32_t r, uint32_t c> template<uint32_t d> MatrixFixed<r, d> MatrixFixed<r, c>::operator*(const MatrixFixed<c, d>& M) const
	{
		double v[r];
		MatrixFixed<r, d> R;
		for(uint32_t j = 0; j < d; j++)
		{
			for(uint32_t i = 0; i < r; i++) v[i] = 0;
			for(uint32_t k = 0; k < c; k++)
			{
				const double s = M.m_data_ref[k + c * j];
				for(uint32_t i = 0; i < r; i++) v[i] += m_data_ref[i + r * k] * s;
			}
			for(uint32_t i = 0; i < r; i++) R.m_data_ptr[i + r * j] = v[i];
		}
		return R;
	}

	//linear
	template<uint32_t r, uint32_t c> MatrixFixed<c, r> MatrixFixed<r, c>::transpose(void) const
	{
		MatrixFixed<c, r> R;
		for(uint32_t i = 0; i < c; i++)
		{
			for(uint32_t j = 0; j < r; i++)
			{
				R.m_data_ptr[i + c * j] = m_data_ref[j + r * i];
			}
		}
		return R;
	}

	//friends
	template<uint32_t r, uint32_t c> MatrixFixed<r, c> operator*(double s, const MatrixFixed<r, c>& M)
	{
		MatrixFixed<r, c> R;
		for(uint32_t i = 0; i < r * c; i++) R.m_data_ptr[i] = s * M.m_data_ref[i];
		return R;
		
	}
}