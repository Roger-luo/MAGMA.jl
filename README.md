# MAGMA

To use the MAGMA package, one should be cautious about the following:

1. The Installation of MAGMA and its location of shared library. With the constant 'libmagma' correctly set (by default const libmagma = "/usr/local/magma/lib/libmagma.so"), one should be able to use ccall to call the functions in MAGMA shared lib.

2. The LOAD_PATH variables. One should push the path to PROJECT/REPOSITORY MAGMA.jl into LOAD_PATH, using 'push!(LOAD_PATH, path/to/the/project)'. And then the test.jl should be able to run successfully
