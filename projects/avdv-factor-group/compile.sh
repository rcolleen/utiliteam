#g++ --std=c++17 lattice.cpp tests.cpp coordinate_test.cpp coordinate.cpp -o coordinate_test
g++ --std=c++17 lattice.cpp tests.cpp lattice_test.cpp -o lattice_test
#g++ --std=c++17 symop.cpp lattice.cpp coordinate.cpp tests.cpp site_test.cpp site.cpp -o site_test
#g++ --std=c++17 symop.cpp lattice.cpp coordinate.cpp site.cpp tests.cpp structure_test.cpp structure.cpp -o structure_test
<<<<<<< HEAD
#g++ --std=c++17 lattice.cpp tests.cpp symop.cpp symop_test.cpp -o test_symop
g++ --std=c++17 lattice.cpp coordinate.cpp site.cpp structure.cpp symop.cpp tests.cpp point_group.cxx test_pointgroup.cpp -o test_pointgroup
#g++ --std=c++17 lattice.cpp site.cpp coordinate.cpp structure.cpp symop.cpp tests.cpp point_group.cxx fastsymmetry.cpp test_fastsym.cpp -o test_fastsymmetry
#g++ --std=c++17 lattice.cpp coordinate.cpp symop.cpp structure.cpp site.cpp tests.cpp io_test.cpp io.hpp -o io_test
#g++ --std=c++17 lattice.cpp site.cpp coordinate.cpp structure.cpp symop.cpp tests.cpp point_group.cxx fastsymmetry.cpp factor_group.cpp factor_group_test.cpp -o test_factorgroup
=======
g++ --std=c++17 lattice.cpp site.cpp coordinate.cpp tests.cpp symop.cpp symop_test.cpp -o test_symop
#g++ --std=c++17 lattice.cpp symop.cpp tests.cpp point_group.cxx test_pointgroup.cpp -o test_pointgroup
#g++ --std=c++17 lattice.cpp site.cpp coordinate.cpp structure.cpp symop.cpp tests.cpp point_group.cxx fastsymmetry.cpp test_fastsym.cpp -o test_fastsymmetry
#g++ --std=c++17 lattice.cpp coordinate.cpp symop.cpp structure.cpp site.cpp tests.cpp io_test.cpp io.hpp -o io_test
g++ --std=c++17 lattice.cpp site.cpp coordinate.cpp structure.cpp symop.cpp tests.cpp point_group.cxx factor_group.cpp factor_group_test.cpp -o test_factorgroup
>>>>>>> upstream/covid_wyckoff
#g++ --std=c++17 lattice.cpp symop.cpp tests.cpp symgroup_test.cpp symgroup.hpp -o symgroup_test
