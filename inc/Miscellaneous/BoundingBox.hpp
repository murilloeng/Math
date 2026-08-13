#pragma once

namespace math
{
	class BoundingBox
	{
	public:
		//constructor
		BoundingBox(void);

		//destructor
		~BoundingBox(void);
		
		//data
		double radius(void) const;
		const double* min(void) const;
		const double* max(void) const;
		
		//status
		bool status(void) const;
		
		//compute
		void reset(void);
		void repair(void);
		void insert_vertex(const double*);

		//print
		void print(const char* = "") const;

		//combine
		void combine(const BoundingBox&);
		static BoundingBox combine(const BoundingBox&, const BoundingBox&);

	private:
		//data
		double m_min[3];
		double m_max[3];
	};
}