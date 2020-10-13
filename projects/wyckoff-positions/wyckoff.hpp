#ifndef WYCKOFF_HH
#define WYCKOFF_HH

#include "../../submodules/eigen-git-mirror/Eigen/Core"
#include "../avdv-factor-group/include.hpp"
#include <vector>

typedef SymGroup<SymOp, BinarySymOpPeriodicCompare_f, BinarySymOpPeriodicMultiplier_f> PeriodicGroup;

class Subspace
{
public:
    /* inherent from symop??*/
    std::vector<Eigen::Vector3d> orthonormalize_basis(const Eigen::Vector3d& basis_vec_0, const Eigen::Vector3d& basis_vec_1, const Eigen::Vector3d& basis_vec_2);
    Subspace(const Eigen::Vector3d& basis_vec_0, const Eigen::Vector3d& basis_vec_1, const Eigen::Vector3d& basis_vec_2, const Eigen::Vector3d& offset);
    Subspace(const Eigen::Matrix3d& basis_col_matrix, const Eigen::Vector3d& offset);
    Subspace();
    int get_dimension();
    std::string formula() const; // Sesha suggested we overload the << operator
    Eigen::Vector3d offset() const;
    Eigen::Matrix3d basis_col_matrix() const; /* [100    [x
                                           010 dot y  + [offset]
                                               001]    z]*/
    Subspace operator*(const SymOp& lhs);
    Eigen::Vector3d m_offset;
    Eigen::Matrix3d m_basis_col_matrix;
    
private:
    int dimension;
};


typedef std::vector<std::vector<Subspace>> WyckoffList;
// Subspace operator*(const SymOp& lhs, const Subspace& rhs);
std::vector<SymOp> find_coset(PeriodicGroup factor_group, PeriodicGroup subgroup);

Eigen::Matrix4d make_symop_4dmatrix(SymOp input_op);

// Reynolds?
Eigen::Matrix4d make_4d_reynolds_operator(PeriodicGroup subgroup);
Eigen::Matrix3d make_3d_reynolds_operator(PeriodicGroup subgroup);
/*template this-->(need begin and end iterator)? vector of symops?*/
Subspace find_invariant_subspace(PeriodicGroup subgroup);

std::vector<Subspace> find_symmetrically_equivalent_wyckoff_positions(std::vector<SymOp> coset,
                                                                      Subspace wyckoff_position);
std::vector<Eigen::Vector3d> convert_subspace_to_vector(Subspace input_subspace, double tol);

Eigen::Vector3d projector_operator(Eigen::Vector3d vect_a, Eigen::Vector3d vect_b);

bool subspaces_are_equal(Subspace lhs, Subspace rhs, double tol);

bool wyckoff_positions_are_equal(std::vector<Subspace> lhs, std::vector<Subspace> rhs, double tol);

////////////////////////

/* ~~~Ouline of functions needed for new wyckoff alorithm~~~~
 *
 * convert symop to fractional coordinates (maybe should be part of symop class, tkae inputs of symop and lattice)
 * combine all symop in group into M-I matrix
 * calc the smith normal form (there are some funcitons to this end in LinearAgebra in CASM, not part of eigen
 * use smith normal form method to calculate wyckoff positions (there are some logic if statements and loops to make sure each one is captured)
 * compare size of wyckoff orbit to sze of coset to determine if the orbit is kept
 * WyckoffList wyckoff_positions generate_space_group_wyckoff_positions(PeriodicGroup point_group);
 */

//////////////////////////////////////////////////////////////

/* ~~THE following should not be necessary given the above algorithm:
 *
 * can make a cyclic subgroup option for generate_subgroups, default not just cyclic.
 * 
 * OtherGroupType center_subgroup(PeriodicGroup cyclic_subgroup);
 * (test that this will work)
 *
 typedef std::vector<std::vector<Subspace>> WyckoffList
 * bool is_wyckoff_unique(WyckoffList group_wyckoff);
 * 
 * WyckoffList extended_wyckoff_list generate_translation_equivalent_wyckoff_lists(WyckoffList group_wyckoff_list);
 * 
 * void intersect_wyckoff_positions(WyckoffList extended_wyckoff_list, WyckoffList group_wyckoff_list);
 *
 * WyckoffList wyckoff_positions generate_point_group_wyckoff_position(PeriodicGroup point_group);
 *

 */
#endif
