#include "./wyckoff.hpp"
#include "../eigen-git-mirror/Eigen/Eigenvalues"
#include <cmath>

std::vector<Eigen::Vector3d> Subspace::orthonormalize_basis(const Eigen::Vector3d& input_vec_0, const Eigen::Vector3d& input_vec_1, const Eigen::Vector3d& input_vec_2)
{
    Eigen::Vector3d output_vec_0=input_vec_0.normalized();
    Eigen::Vector3d output_vec_1=input_vec_1 - projector_operator(input_vec_1,output_vec_0);
    output_vec_1.normalize();
    Eigen::Vector3d output_vec_2= input_vec_2 - projector_operator(input_vec_2, output_vec_0) -projector_operator(input_vec_2, output_vec_1);
    output_vec_2.normalize();
    return {output_vec_0, output_vec_1, output_vec_2};
}

Subspace::Subspace(const Eigen::Vector3d& input_vec_0, const Eigen::Vector3d& input_vec_1,
                   const Eigen::Vector3d& input_vec_2, const Eigen::Vector3d& offset){
    std::vector<Eigen::Vector3d> basis_vectors=orthonormalize_basis( input_vec_0, input_vec_1, input_vec_2);

    m_offset = offset -projector_operator(offset, basis_vectors[0]) - projector_operator(offset,basis_vectors[1]) - projector_operator(offset, basis_vectors[2]);
    for (int i = 0; i < 3; ++i)
    {
        m_basis_col_matrix(i, 0) = basis_vectors[0](i);
        m_basis_col_matrix(i, 1) = basis_vectors[1](i);
        m_basis_col_matrix(i, 2) = basis_vectors[2](i);
    }

    //This may still be a different subspace though order if given axis, or plane could be define in different orders
    //DONE: make sure planes are defined by orthoganol vectors?
    //      ake sure offset is not in subspace
}

Subspace::Subspace(const Eigen::Matrix3d& input_basis_vectors, const Eigen::Vector3d& input_offset)
    :Subspace(input_basis_vectors.col(0), input_basis_vectors.col(1), input_basis_vectors.col(2), input_offset) {}

Subspace::Subspace(): m_offset(Eigen::Vector3d::Zero()) , m_basis_col_matrix(Eigen::Matrix3d::Zero()){}

Eigen::Matrix3d Subspace::basis_col_matrix() const { return this->m_basis_col_matrix; }

Eigen::Vector3d Subspace::offset() const { return this->m_offset; }

Subspace Subspace::operator*(const SymOp& lhs)
{
    Eigen::Matrix3d basis_vectors = lhs.get_cart_matrix() * this->m_basis_col_matrix;
    Eigen::Vector3d offset_vector = lhs.get_cart_matrix() * this->m_offset + lhs.get_translation();
    Subspace product_wycoff_position(basis_vectors, offset_vector);
    return product_wycoff_position;
}

std::string Subspace::formula() const
{
    double tol = 1e-5;
    std::string formula = "{";
    std::vector<std::string> xyz = {"x", "y", "z"};

    for (int i = 0; i < 3; ++i)
    {
        std::string temp_string;
        for (int j = 0; j < 3; ++j)
        {
            if (std::abs(this->m_basis_col_matrix(i,j)) > tol)
            {
                std::string sign = this->m_basis_col_matrix(i,j) > 0 ? "+" : "";
                temp_string += sign + std::to_string(this->m_basis_col_matrix(i,j)) + xyz[j];
            }
        }
        if (temp_string.size() == 0)
        {
            temp_string = "0";
        }
        formula += temp_string + ", ";
    }

    formula.pop_back();
    formula.pop_back();
    formula += "}";
    return formula;
}

std::vector<SymOp> find_coset(PeriodicGroup factor_group, PeriodicGroup subgroup)
{
    /*finds the coset to a subgroup in a factor group and returns it as
     * a std::vector of SymOps*/
    std::vector<SymOp> coset;
    auto binary_compare = factor_group.get_comparator();
    for (SymOp factor_group_symop : factor_group.operations())
    {
        auto compare = [factor_group_symop, binary_compare](SymOp subgroup_symop) {
            return binary_compare(factor_group_symop, subgroup_symop);
        };
        if (std::find_if(subgroup.operations().begin(), subgroup.operations().end(), compare) ==
            subgroup.operations().end())
        {
            coset.push_back(factor_group_symop);
        }
    }
    return coset;
}

Eigen::Matrix4d make_symop_4dmatrix(SymOp input_symop)
{
    Eigen::Matrix4d symop_4d;
    symop_4d.row(3) << 0, 0, 0, 1;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            symop_4d(i, j) = input_symop.get_cart_matrix()(i, j);
        }
    }
    for (int i = 0; i < 3; i++)
    {
        symop_4d(i, 3) = input_symop.get_translation()(i);
    }
    return symop_4d;
}

Eigen::Matrix4d make_4d_reynolds_operator(PeriodicGroup subgroup)
{
    Eigen::Matrix4d average_op = Eigen::Matrix4d::Zero();
    for (SymOp sym_op : subgroup.operations())
    {
        average_op += make_symop_4dmatrix(sym_op);
    }
    average_op = average_op / subgroup.operations().size();
    return average_op;
}

/*Subspace find_invariant_subspace(PeriodicGroup subgroup){
 return;
}*/

// Subspace operator*(const SymOp& lhs, const Subspace& rhs)
//{
//    Eigen::Matrix3d basis_vectors = lhs.get_cart_matrix() * rhs.basis_col_matrix();
//    Eigen::Vector3d offset_vector = lhs.get_cart_matrix() * rhs.offset() + lhs.get_translation();
//    Subspace product_wycoff_position(basis_vectors, offset_vector);
//    return product_wycoff_position;
//}

Eigen::Matrix3d make_3d_reynolds_operator(PeriodicGroup subgroup)
{
    Eigen::Matrix3d average_op = Eigen::Matrix3d::Zero();
    for (const auto& sym_op : subgroup.operations())
    {
        average_op += sym_op.get_cart_matrix();
    }

    average_op = average_op / subgroup.operations().size();
    return average_op;
}

Subspace find_invariant_subspace(PeriodicGroup subgroup)
{
    // TODO: How is tolerance is being dealt? Any global variable to use or should it be an input arg?
    double tol = 1e-5;
    Eigen::Matrix3d reynolds_operator = make_3d_reynolds_operator(subgroup);
    Eigen::EigenSolver<Eigen::Matrix3d> eigen_solver(reynolds_operator);
    Eigen::VectorXcd eigen_values = eigen_solver.eigenvalues();
    Eigen::MatrixXcd eigen_vectors = eigen_solver.eigenvectors();
    Eigen::Matrix3d basis_col_matrix = Eigen::Matrix3d::Zero();
    // TODO: Should compute the average translation in subgroup to get the offset
    Eigen::Vector3d offset(0, 0, 0);

    for (int i = 0; i < 3; ++i)
    {
        if (std::abs(eigen_values(i).real() - 1) < tol && std::abs(eigen_values(i).imag()) < tol)
        {
            for (int j = 0; j < 3; ++j)
            {
                basis_col_matrix(j, i) = eigen_vectors(j, i).real();
            }
        }
    }

    Subspace invariant_subspace(basis_col_matrix, offset);
    return invariant_subspace;
}

std::vector<Subspace> find_symmetrically_equivalent_wyckoff_positions(std::vector<SymOp> coset,
                                                                      Subspace wyckoff_position)
{
    std::vector<Subspace> equivalent_wyckoff_positions;
    equivalent_wyckoff_positions.push_back(wyckoff_position);
    for (SymOp symop : coset)
    {
        equivalent_wyckoff_positions.push_back(wyckoff_position * symop);
        //TODO: check if in list before push_back
    }

    return equivalent_wyckoff_positions;
}

std::vector<Eigen::Vector3d> convert_subspace_to_vector(Subspace input_subspace, double tol)
{
    std::vector<Eigen::Vector3d> output_vector;
    for( int i=0; i<3; i++)
    {       Eigen::Vector3d temp_vector;
        temp_vector<<input_subspace.basis_col_matrix()(0,i), input_subspace.basis_col_matrix()(1,i), input_subspace.basis_col_matrix()(2,i);
        if ( temp_vector.norm()>tol)
            {output_vector.push_back(temp_vector);}
    }
    return output_vector;
}

/*moved function to math.hh in avdv-factor-group
 * Eigen::Vector3d projector_operator(Eigen::Vector3d vect_a, Eigen::Vector3d vect_b)
{
    //project a onto b
    return ( vect_b*vect_b.transpose()*vect_a);
}
*/

bool subspaces_are_equal(Subspace lhs, Subspace rhs, double tol)
{
// different depending on the dimension of the subspace.
//
    std::vector<Eigen::Vector3d> subspace_1_vectors = convert_subspace_to_vector(lhs, tol);
    std::vector<Eigen::Vector3d> subspace_2_vectors = convert_subspace_to_vector(rhs, tol);

    //subspaces have different dimensions:
    if (subspace_1_vectors.size() !=subspace_2_vectors.size()){ return false;}

    //subspaces are both dimension 3:
    //here they shouls have exactly the same matrix and offset
    if (subspace_1_vectors.size()==4 || subspace_1_vectors.size()==1)
    { return (lhs.basis_col_matrix().isApprox(rhs.basis_col_matrix(), tol) && lhs.offset().isApprox(rhs.offset(), tol));}

    //subspaces are both planes:
    //here they should have the same normal vector and their offsets should have the components that are normal to the plane
    if (subspace_1_vectors.size() ==3)
    {
        Eigen::Vector3d cross_vector_1=subspace_1_vectors[0].cross(subspace_1_vectors[1]).normalized();
        Eigen::Vector3d cross_vector_2=subspace_2_vectors[0].cross(subspace_2_vectors[1]).normalized();
        if (compare_vectors(cross_vector_1, cross_vector_2, tol)){
            return compare_vectors(projector_operator(subspace_1_vectors[2], cross_vector_1), projector_operator(subspace_2_vectors[2], cross_vector_2), tol);
                    }
        return false;
     }

    //subspaces are both vectors:
    //here the dot of both vectors hsould be 1 when normalized, and the offset should have a component perpendicular to the vector which is identical.
    if (subspace_1_vectors.size() ==2)
    {
        if (abs(subspace_1_vectors[0].normalized().dot(subspace_2_vectors[0].normalized()))-1<tol)
        {
           Eigen::Vector3d offset_1= subspace_1_vectors[1]-projector_operator(subspace_1_vectors[1], subspace_1_vectors[0]);
           Eigen::Vector3d offset_2= subspace_2_vectors[1]-projector_operator(subspace_2_vectors[1], subspace_2_vectors[0]);
           return compare_vectors(offset_1, offset_2, tol);
           
        }
        return false;
    }


    //TODO: insufficient comparison, multiple vectors can span the same subspace
    // is A is contained in B and B is contained in A then they are the same subspace (projections?)
}

bool wyckoff_positions_are_equal(std::vector<Subspace> lhs, std::vector<Subspace> rhs, double tol)
{ //TODO: 
    return false;
}
