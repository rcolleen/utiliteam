#include "symop.hpp"
#include "tests.hpp"
#include <type_traits>
#include <iostream>
int main()
{
	//make symop with translation and test
    std::cout<<"---- Running SymOp Tests ----"<<std::endl;
    std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl; 
    
    Eigen::Matrix3d test_matrix= make_z_rotation_matrix(60);
	Eigen::Vector3d translation;
	translation<<0,0,1;
	SymOp With_Rotations_60= SymOp(test_matrix);
	SymOp With_Rotations_screw60= SymOp(test_matrix, translation);	
	EXPECT_TRUE(With_Rotations_60.get_cart_matrix()== test_matrix, "test matrix from default construction without translation");
	EXPECT_TRUE(With_Rotations_screw60.get_cart_matrix()== test_matrix, "test matrix from default construction with translation");
	EXPECT_TRUE(With_Rotations_60.get_translation()== Eigen::Vector3d(0,0,0), "translation from default construction without translation");
	EXPECT_TRUE(With_Rotations_screw60.get_translation()== translation, "transaltion from default construction with translation");
	//test both the compare functions
	
	SymOpCompare_f compare_60rot(With_Rotations_60);
	EXPECT_TRUE(compare_60rot(With_Rotations_60), "60 rotations should be itself");
	EXPECT_TRUE(!compare_60rot(With_Rotations_screw60), "60 rotations should not be screw rotation");
	//test the multiplication function
	Eigen::Matrix3d rot_120=make_z_rotation_matrix(120);
	SymOp With_Rotations_120= SymOp(rot_120);
	auto With_Rotations_60_product=With_Rotations_60*With_Rotations_60;
	EXPECT_TRUE(With_Rotations_120.get_cart_matrix().isApprox(With_Rotations_60_product.get_cart_matrix(), 1e-6), "120 degrees is two 60 degree rotations");
	return 0;
}