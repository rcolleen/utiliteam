#include "./categorize_symop.hh"
#include "../avdv-factor-group/tests.hpp"

bool test_categorize_identity(double tol)
{
    Lattice lattice({1,0,0}, {0,1,0}, {0,0,1});
    Eigen::Matrix3d identity_matrix=Eigen::Matrix3d::Identity();
    SymOp identity_op(identity_matrix);
    if (check_op_type(identity_op, lattice, tol)==SYMOP_TYPE::IDENTITY)
    {return true;}
    return false;
}

bool test_categorize_rotation(double tol)
{
    Lattice lattice({1,0,0}, {0,1,0}, {0,0,1});
    Eigen::Matrix3d rotation_matrix=make_z_rotation_matrix(60);
    SymOp rotation_op(rotation_matrix);
    if (check_op_type(rotation_op, lattice, tol)==SYMOP_TYPE::ROTATION)
    {return true;}
    return false;
}

bool test_categorize_screw(double tol)
{
    Lattice lattice({1,0,0}, {0,1,0}, {0,0,1});
    Eigen::Matrix3d rotation_matrix=make_z_rotation_matrix(60);
    Eigen::Vector3d translation;
    translation<<0.0, 0.0, 0.5;
    SymOp screw_op(rotation_matrix, translation);
    if (check_op_type(screw_op, lattice, tol)==SYMOP_TYPE::SCREW)
    {return true;}
    return false;
}

bool test_categorize_not_screw(double tol)
{
    Lattice lattice({1,0,0}, {0,1,0}, {0,0,1});
    Eigen::Matrix3d rotation_matrix=make_z_rotation_matrix(60);
    Eigen::Vector3d translation;
    translation<<0.0, 0.0, 1.0;
    SymOp screw_op(rotation_matrix, translation);
    if (check_op_type(screw_op, lattice, tol)==SYMOP_TYPE::SCREW)
    {return false;}
    if (check_op_type(screw_op, lattice, tol)==SYMOP_TYPE::ROTATION)
    {return true;}
    return false;
}

bool test_invariant_subspace_rotation(double tol)
{
    Lattice lattice({1,0,0}, {0,1,0}, {0,0,1});
    Eigen::Matrix3d rotation_matrix=make_z_rotation_matrix(60);
    SymOp rotation_op(rotation_matrix);
    
    Subspace expected_subspace({0,0,1}, {0,0,0},{0,0,0},{0,0,0});
    std::optional<Subspace> invariant_subspace=find_invariant_subspace(rotation_op,lattice,tol); 
    if (invariant_subspace.has_value()){
            std::cout<<"Calclulated_invariant_subspace"<<invariant_subspace.value().formula()<<std::endl;
            return subspaces_are_equal(expected_subspace, invariant_subspace.value(), tol);
    }
    return false;
}


bool test_invariant_subspace_screw(double tol)
{
    Lattice lattice({1,0,0}, {0,1,0}, {0,0,1});
    Eigen::Matrix3d rotation_matrix=make_z_rotation_matrix(60);
    Eigen::Vector3d translation;
    translation<<0.0, 0.0, 0.5;
    SymOp screw_op(rotation_matrix, translation);
    Subspace expected_subspace({0,0,0}, {0,0,0},{0,0,1},{0,0,0});
    std::optional<Subspace> invariant_subspace=find_invariant_subspace(screw_op,lattice,tol); 
    
    return (invariant_subspace.has_value()){
        return false;
    }
    return true;
}

int main()
{ 
    double tol=1e-6;
    EXPECT_TRUE(test_categorize_identity(tol), "test categorize identity");
    EXPECT_TRUE(test_categorize_rotation(tol), "test categorize rotation");
    EXPECT_TRUE(test_categorize_screw(tol), "test categorize screw");
    EXPECT_TRUE(test_categorize_not_screw(tol), "test categorize not screw");
    EXPECT_TRUE(test_invariant_subspace_rotation(tol), "test find invariant subspace of rotation");
    EXPECT_TRUE(test_invariant_subspace_screw(tol), "test find invariant subspace screw");

    return 0;

}
