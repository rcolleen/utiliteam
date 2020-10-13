#ifndef MATH_HH
#define MATH_HH
#include "../../submodules/eigen-git-mirror/Eigen/Core"

template <typename MatrixType>
bool is_unitary(const MatrixType& input_matrix, double tol)
{
    MatrixType product = input_matrix.transpose() * input_matrix;
    return product.isIdentity(tol);
}

Eigen::Vector3d projector_operator(Eigen::Vector3d vect_a, Eigen::Vector3d vect_b)
{
    //returns projection of vect_a onto vect_b
    return (vect_b*vect_b.transpose()*vect_a);
}
#endif
