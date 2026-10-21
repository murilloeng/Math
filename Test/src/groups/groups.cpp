//std
#include <cmath>
#include <ctime>

//Math
#include "Math/inc/Groups/ASO2.hpp"
#include "Math/inc/Groups/GSO2.hpp"

#include "Math/inc/Groups/ASE2.hpp"
#include "Math/inc/Groups/GSE2.hpp"

#include "Math/inc/Groups/ASO3.hpp"
#include "Math/inc/Groups/GSO3.hpp"

#include "Math/inc/Groups/ASE3.hpp"
#include "Math/inc/Groups/GSE3.hpp"

//Test
#include "Math/Test/inc/groups.hpp"

void tests::groups::groups(void)
{
	srand(time(nullptr));
	math::ASO2::test_adjoint(1000, 1e-12);
	math::ASO2::test_bracket(1000, 1e-12);
	math::ASO2::test_exp_log(1000, 1e-12);
	math::ASO2::test_jacobian(1000, 1e-08);
	math::ASO2::test_jacobian_inverse(1000, 1e-12);

	math::GSO2::test_matrix(1000, 1e-12);
	math::GSO2::test_inverse(1000, 1e-12);
	math::GSO2::test_log_exp(1000, 1e-12);

	math::ASE2::test_adjoint(1000, 1e-12);
	math::ASE2::test_bracket(1000, 1e-12);
	math::ASE2::test_exp_log(1000, 1e-12);
	math::ASE2::test_jacobian(1000, 1e-08);
	math::ASE2::test_jacobian_inverse(1000, 1e-12);

	math::GSE2::test_matrix(1000, 1e-12);
	math::GSE2::test_inverse(1000, 1e-12);
	math::GSE2::test_log_exp(1000, 1e-12);

	math::ASO3::test_adjoint(1000, 1e-12);
	math::ASO3::test_bracket(1000, 1e-12);
	math::ASO3::test_exp_log(1000, 1e-12);
	math::ASO3::test_jacobian(1000, 1e-08);
	math::ASO3::test_jacobian_inverse(1000, 1e-12);

	math::GSO3::test_matrix(1000, 1e-12);
	math::GSO3::test_inverse(1000, 1e-12);
	math::GSO3::test_log_exp(1000, 1e-12);

	math::ASE3::test_adjoint(1000, 1e-12);
	math::ASE3::test_bracket(1000, 1e-12);
	math::ASE3::test_exp_log(1000, 1e-12);
	math::ASE3::test_jacobian(1000, 1e-08);
	math::ASE3::test_jacobian_inverse(1000, 1e-12);

	math::GSE3::test_matrix(1000, 1e-12);
	math::GSE3::test_inverse(1000, 1e-12);
	math::GSE3::test_log_exp(1000, 1e-12);
}