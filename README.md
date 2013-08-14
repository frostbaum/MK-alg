MK-alg
======

O(3) implementation transferred to Fortran from:

src/main/java/blogspot/software_and_algorithms/stern_library/optimization/HungarianAlgorithm.java

Use:
integer :: dim                            ! dimension of problem
double precision :: costmatrix(dim, dim)  ! cost matrix
integer :: solution(dim)                  ! assigned column index by row

call mk_ass(dim, costmatrix, solution)
