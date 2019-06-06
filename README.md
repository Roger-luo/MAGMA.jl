#MAGMA

A Julia binding for MAGMA, don't drink it.

Disclaimer: This is an on-going JSoC 2019 project, it is in very early stage and still work in progress.
Installation

To use this package, use dev in Pkg mode.

(pkg)> dev https://github.com/Roger-luo/MAGMA.jl.git

The Installation of MAGMA and its location of shared library. With the constant 'libmagma' correctly set (by default const libmagma = "/usr/local/magma/lib/libmagma.so"), one should be able to use ccall to call the functions in MAGMA shared lib.
License

MIT
