# Automatically generated using Clang.jl


const MAGMA_API = 2

# Skipping MacroDefinition: magma_free ( ptr ) magma_free_internal ( ptr , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_free_pinned ( ptr ) magma_free_pinned_internal ( ptr , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_queue_create ( device , queue_ptr ) magma_queue_create_internal ( device , queue_ptr , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_queue_create_from_cuda ( device , cuda_stream , cublas_handle , cusparse_handle , queue_ptr ) magma_queue_create_from_cuda_internal ( device , cuda_stream , cublas_handle , cusparse_handle , queue_ptr , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_queue_create_from_opencl ( device , cl_queue , queue_ptr ) magma_queue_create_from_opencl_internal ( device , cl_queue , queue_ptr , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_queue_destroy ( queue ) magma_queue_destroy_internal ( queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_queue_sync ( queue ) magma_queue_sync_internal ( queue , __func__ , __FILE__ , __LINE__ )

const MAX_THREADS_BLG = 256
const magmaFloatComplex = Cint

@cenum magma_side_t::UInt32 begin
    MagmaLeft = 141
    MagmaRight = 142
    MagmaBothSides = 143
end


const real_Double_t = Cdouble

struct gbstrct_blg
    dQ1::PtrOrCuPtr{magmaFloatComplex}
    dT1::PtrOrCuPtr{magmaFloatComplex}
    dT2::PtrOrCuPtr{magmaFloatComplex}
    dV2::PtrOrCuPtr{magmaFloatComplex}
    dE::PtrOrCuPtr{magmaFloatComplex}
    T::PtrOrCuPtr{magmaFloatComplex}
    A::PtrOrCuPtr{magmaFloatComplex}
    V::PtrOrCuPtr{magmaFloatComplex}
    TAU::PtrOrCuPtr{magmaFloatComplex}
    E::PtrOrCuPtr{magmaFloatComplex}
    E_CPU::PtrOrCuPtr{magmaFloatComplex}
    cores_num::Cint
    locores_num::Cint
    overlapQ1::Cint
    usemulticpu::Cint
    NB::Cint
    NBTILES::Cint
    N::Cint
    NE::Cint
    N_CPU::Cint
    N_GPU::Cint
    LDA::Cint
    LDE::Cint
    BAND::Cint
    grsiz::Cint
    Vblksiz::Cint
    WANTZ::Cint
    SIDE::magma_side_t
    timeblg::PtrOrCuPtr{real_Double_t}
    timeaplQ::PtrOrCuPtr{real_Double_t}
    ss_prog::PtrOrCuPtr{Cint}
end

const magma_int_t = Cint
const magmaFloatComplex_ptr = PtrOrCuPtr{magmaFloatComplex}
const magma_queue = Cvoid
const magma_queue_t = PtrOrCuPtr{magma_queue}

struct cgehrd_data
    ngpu::magma_int_t
    ldda::magma_int_t
    ldv::magma_int_t
    ldvd::magma_int_t
    dA::NTuple{8, magmaFloatComplex_ptr}
    dV::NTuple{8, magmaFloatComplex_ptr}
    dVd::NTuple{8, magmaFloatComplex_ptr}
    dY::NTuple{8, magmaFloatComplex_ptr}
    dW::NTuple{8, magmaFloatComplex_ptr}
    dTi::NTuple{8, magmaFloatComplex_ptr}
    queues::NTuple{8, magma_queue_t}
end

const magmaDouble_ptr = PtrOrCuPtr{Cdouble}

struct dgehrd_data
    ngpu::magma_int_t
    ldda::magma_int_t
    ldv::magma_int_t
    ldvd::magma_int_t
    dA::NTuple{8, magmaDouble_ptr}
    dV::NTuple{8, magmaDouble_ptr}
    dVd::NTuple{8, magmaDouble_ptr}
    dY::NTuple{8, magmaDouble_ptr}
    dW::NTuple{8, magmaDouble_ptr}
    dTi::NTuple{8, magmaDouble_ptr}
    queues::NTuple{8, magma_queue_t}
end

const magmaFloat_ptr = PtrOrCuPtr{Cfloat}

struct sgehrd_data
    ngpu::magma_int_t
    ldda::magma_int_t
    ldv::magma_int_t
    ldvd::magma_int_t
    dA::NTuple{8, magmaFloat_ptr}
    dV::NTuple{8, magmaFloat_ptr}
    dVd::NTuple{8, magmaFloat_ptr}
    dY::NTuple{8, magmaFloat_ptr}
    dW::NTuple{8, magmaFloat_ptr}
    dTi::NTuple{8, magmaFloat_ptr}
    queues::NTuple{8, magma_queue_t}
end

const MAGMA_D_ZERO = 0.0
const MAGMA_D_ONE = 1.0
const MAGMA_D_HALF = 0.5
const MAGMA_D_NEG_ONE = -1.0
const MAGMA_D_NEG_HALF = -0.5
const MAGMA_S_ZERO = 0.0
const MAGMA_S_ONE = 1.0
const MAGMA_S_HALF = 0.5
const MAGMA_S_NEG_ONE = -1.0
const MAGMA_S_NEG_HALF = -0.5

const MAGMA_VERSION_MAJOR = 2
const MAGMA_VERSION_MINOR = 5
const MAGMA_VERSION_MICRO = 1
const MAGMA_VERSION_STAGE = "alpha1"
const MagmaMaxGPUs = 8
const MagmaMaxAccelerators = 8
const MagmaMaxSubs = 16
const MagmaBigTileSize = 1000000
const MAGMA_SUCCESS = 0
const MAGMA_ERR = -100
const MAGMA_ERR_NOT_INITIALIZED = -101
const MAGMA_ERR_REINITIALIZED = -102
const MAGMA_ERR_NOT_SUPPORTED = -103
const MAGMA_ERR_ILLEGAL_VALUE = -104
const MAGMA_ERR_NOT_FOUND = -105
const MAGMA_ERR_ALLOCATION = -106
const MAGMA_ERR_INTERNAL_LIMIT = -107
const MAGMA_ERR_UNALLOCATED = -108
const MAGMA_ERR_FILESYSTEM = -109
const MAGMA_ERR_UNEXPECTED = -110
const MAGMA_ERR_SEQUENCE_FLUSHED = -111
const MAGMA_ERR_HOST_ALLOC = -112
const MAGMA_ERR_DEVICE_ALLOC = -113
const MAGMA_ERR_CUDASTREAM = -114
const MAGMA_ERR_INVALID_PTR = -115
const MAGMA_ERR_UNKNOWN = -116
const MAGMA_ERR_NOT_IMPLEMENTED = -117
const MAGMA_ERR_NAN = -118
const MAGMA_SLOW_CONVERGENCE = -201
const MAGMA_DIVERGENCE = -202
const MAGMA_NONSPD = -203
const MAGMA_ERR_BADPRECOND = -204
const MAGMA_NOTCONVERGED = -205
const MAGMA_ERR_CUSPARSE = -3000
const MAGMA_ERR_CUSPARSE_NOT_INITIALIZED = -3001
const MAGMA_ERR_CUSPARSE_ALLOC_FAILED = -3002
const MAGMA_ERR_CUSPARSE_INVALID_VALUE = -3003
const MAGMA_ERR_CUSPARSE_ARCH_MISMATCH = -3004
const MAGMA_ERR_CUSPARSE_MAPPING_ERROR = -3005
const MAGMA_ERR_CUSPARSE_EXECUTION_FAILED = -3006
const MAGMA_ERR_CUSPARSE_INTERNAL_ERROR = -3007
const MAGMA_ERR_CUSPARSE_MATRIX_TYPE_NOT_SUPPORTED = -3008
const MAGMA_ERR_CUSPARSE_ZERO_PIVOT = -3009

@cenum magma_bool_t::UInt32 begin
    MagmaFalse = 0
    MagmaTrue = 1
end


const Magma2lapack_Min = MagmaFalse

@cenum magma_storev_t::UInt32 begin
    MagmaColumnwise = 401
    MagmaRowwise = 402
end


const Magma2lapack_Max = MagmaRowwise
const MagmaRowMajorStr = "Row"
const MagmaColMajorStr = "Col"
const MagmaNoTransStr = "NoTrans"
const MagmaTransStr = "Trans"
const MagmaConjTransStr = "ConjTrans"
const Magma_ConjTransStr = "ConjTrans"
const MagmaUpperStr = "Upper"
const MagmaLowerStr = "Lower"
const MagmaFullStr = "Full"
const MagmaNonUnitStr = "NonUnit"
const MagmaUnitStr = "Unit"
const MagmaLeftStr = "Left"
const MagmaRightStr = "Right"
const MagmaBothSidesStr = "Both"
const MagmaOneNormStr = "1"
const MagmaTwoNormStr = "2"
const MagmaFrobeniusNormStr = "Fro"
const MagmaInfNormStr = "Inf"
const MagmaMaxNormStr = "Max"
const MagmaForwardStr = "Forward"
const MagmaBackwardStr = "Backward"
const MagmaColumnwiseStr = "Columnwise"
const MagmaRowwiseStr = "Rowwise"
const MagmaNoVecStr = "NoVec"
const MagmaVecStr = "Vec"
const MagmaIVecStr = "IVec"
const MagmaAllVecStr = "All"
const MagmaSomeVecStr = "Some"
const MagmaOverwriteVecStr = "Overwrite"
const magma_index_t = Cint
const magma_uindex_t = UInt32
const magma_event_t = Cint
const magma_device_t = magma_int_t
const magmaHalf = Int16
const magmaDoubleComplex = Cint
const magma_ptr = PtrOrCuPtr{Cvoid}
const magmaInt_ptr = PtrOrCuPtr{magma_int_t}
const magmaIndex_ptr = PtrOrCuPtr{magma_index_t}
const magmaUIndex_ptr = PtrOrCuPtr{magma_uindex_t}
const magmaDoubleComplex_ptr = PtrOrCuPtr{magmaDoubleComplex}
const magmaHalf_ptr = PtrOrCuPtr{magmaHalf}
const magma_const_ptr = PtrOrCuPtr{Cvoid}
const magmaInt_const_ptr = PtrOrCuPtr{magma_int_t}
const magmaIndex_const_ptr = PtrOrCuPtr{magma_index_t}
const magmaUIndex_const_ptr = PtrOrCuPtr{magma_uindex_t}
const magmaFloat_const_ptr = PtrOrCuPtr{Cfloat}
const magmaDouble_const_ptr = PtrOrCuPtr{Cdouble}
const magmaFloatComplex_const_ptr = PtrOrCuPtr{magmaFloatComplex}
const magmaDoubleComplex_const_ptr = PtrOrCuPtr{magmaDoubleComplex}
const magmaHalf_const_ptr = PtrOrCuPtr{magmaHalf}

@cenum magma_order_t::UInt32 begin
    MagmaRowMajor = 101
    MagmaColMajor = 102
end

@cenum magma_trans_t::UInt32 begin
    MagmaNoTrans = 111
    MagmaTrans = 112
    MagmaConjTrans = 113
    Magma_ConjTrans = 113
end

@cenum magma_uplo_t::UInt32 begin
    MagmaUpper = 121
    MagmaLower = 122
    MagmaFull = 123
    MagmaHessenberg = 124
end


const magma_type_t = magma_uplo_t

@cenum magma_diag_t::UInt32 begin
    MagmaNonUnit = 131
    MagmaUnit = 132
end

@cenum magma_norm_t::UInt32 begin
    MagmaOneNorm = 171
    MagmaRealOneNorm = 172
    MagmaTwoNorm = 173
    MagmaFrobeniusNorm = 174
    MagmaInfNorm = 175
    MagmaRealInfNorm = 176
    MagmaMaxNorm = 177
    MagmaRealMaxNorm = 178
end

@cenum magma_dist_t::UInt32 begin
    MagmaDistUniform = 201
    MagmaDistSymmetric = 202
    MagmaDistNormal = 203
end

@cenum magma_sym_t::UInt32 begin
    MagmaHermGeev = 241
    MagmaHermPoev = 242
    MagmaNonsymPosv = 243
    MagmaSymPosv = 244
end

@cenum magma_pack_t::UInt32 begin
    MagmaNoPacking = 291
    MagmaPackSubdiag = 292
    MagmaPackSupdiag = 293
    MagmaPackColumn = 294
    MagmaPackRow = 295
    MagmaPackLowerBand = 296
    MagmaPackUpeprBand = 297
    MagmaPackAll = 298
end

@cenum magma_vec_t::UInt32 begin
    MagmaNoVec = 301
    MagmaVec = 302
    MagmaIVec = 303
    MagmaAllVec = 304
    MagmaSomeVec = 305
    MagmaOverwriteVec = 306
    MagmaBacktransVec = 307
end

@cenum magma_range_t::UInt32 begin
    MagmaRangeAll = 311
    MagmaRangeV = 312
    MagmaRangeI = 313
end

@cenum magma_vect_t::UInt32 begin
    MagmaQ = 322
    MagmaP = 323
end

@cenum magma_direct_t::UInt32 begin
    MagmaForward = 391
    MagmaBackward = 392
end

@cenum magma_mode_t::UInt32 begin
    MagmaHybrid = 701
    MagmaNative = 702
end

@cenum magma_storage_t::UInt32 begin
    Magma_CSR = 611
    Magma_ELLPACKT = 612
    Magma_ELL = 613
    Magma_DENSE = 614
    Magma_BCSR = 615
    Magma_CSC = 616
    Magma_HYB = 617
    Magma_COO = 618
    Magma_ELLRT = 619
    Magma_SPMVFUNCTION = 620
    Magma_SELLP = 621
    Magma_ELLD = 622
    Magma_CSRLIST = 623
    Magma_CSRD = 624
    Magma_CSRL = 627
    Magma_CSRU = 628
    Magma_CSRCOO = 629
    Magma_CUCSR = 630
    Magma_COOLIST = 631
    Magma_CSR5 = 632
end

@cenum magma_solver_type::UInt32 begin
    Magma_CG = 431
    Magma_CGMERGE = 432
    Magma_GMRES = 433
    Magma_BICGSTAB = 434
    Magma_BICGSTABMERGE = 435
    Magma_BICGSTABMERGE2 = 436
    Magma_JACOBI = 437
    Magma_GS = 438
    Magma_ITERREF = 439
    Magma_BCSRLU = 440
    Magma_PCG = 441
    Magma_PGMRES = 442
    Magma_PBICGSTAB = 443
    Magma_PASTIX = 444
    Magma_ILU = 445
    Magma_ICC = 446
    Magma_PARILU = 447
    Magma_PARIC = 448
    Magma_BAITER = 449
    Magma_LOBPCG = 450
    Magma_NONE = 451
    Magma_FUNCTION = 452
    Magma_IDR = 453
    Magma_PIDR = 454
    Magma_CGS = 455
    Magma_PCGS = 456
    Magma_CGSMERGE = 457
    Magma_PCGSMERGE = 458
    Magma_TFQMR = 459
    Magma_PTFQMR = 460
    Magma_TFQMRMERGE = 461
    Magma_PTFQMRMERGE = 462
    Magma_QMR = 463
    Magma_PQMR = 464
    Magma_QMRMERGE = 465
    Magma_PQMRMERGE = 466
    Magma_BOMBARD = 490
    Magma_BOMBARDMERGE = 491
    Magma_PCGMERGE = 492
    Magma_BAITERO = 493
    Magma_IDRMERGE = 494
    Magma_PBICGSTABMERGE = 495
    Magma_PARICT = 496
    Magma_CUSTOMIC = 497
    Magma_CUSTOMILU = 498
    Magma_PIDRMERGE = 499
    Magma_BICG = 500
    Magma_BICGMERGE = 501
    Magma_PBICG = 502
    Magma_PBICGMERGE = 503
    Magma_LSQR = 504
    Magma_PARILUT = 505
    Magma_ISAI = 506
    Magma_CUSOLVE = 507
    Magma_VBJACOBI = 508
    Magma_PARDISO = 509
    Magma_SYNCFREESOLVE = 510
    Magma_ILUT = 511
end

@cenum magma_ortho_t::UInt32 begin
    Magma_CGSO = 561
    Magma_FUSED_CGSO = 562
    Magma_MGSO = 563
end

@cenum magma_location_t::UInt32 begin
    Magma_CPU = 571
    Magma_DEV = 572
end

@cenum magma_symmetry_t::UInt32 begin
    Magma_GENERAL = 581
    Magma_SYMMETRIC = 582
end

@cenum magma_diagorder_t::UInt32 begin
    Magma_ORDERED = 591
    Magma_DIAGFIRST = 592
    Magma_UNITY = 593
    Magma_VALUE = 594
end

@cenum magma_precision::UInt32 begin
    Magma_DCOMPLEX = 501
    Magma_FCOMPLEX = 502
    Magma_DOUBLE = 503
    Magma_FLOAT = 504
end

@cenum magma_scale_t::UInt32 begin
    Magma_NOSCALE = 511
    Magma_UNITROW = 512
    Magma_UNITDIAG = 513
    Magma_UNITCOL = 514
    Magma_UNITROWCOL = 515
    Magma_UNITDIAGCOL = 516
end

@cenum magma_operation_t::UInt32 begin
    Magma_SOLVE = 801
    Magma_SETUPSOLVE = 802
    Magma_APPLYSOLVE = 803
    Magma_DESTROYSOLVE = 804
    Magma_INFOSOLVE = 805
    Magma_GENERATEPREC = 806
    Magma_PRECONDLEFT = 807
    Magma_PRECONDRIGHT = 808
    Magma_TRANSPOSE = 809
    Magma_SPMV = 810
end

@cenum magma_refinement_t::UInt32 begin
    Magma_PREC_SS = 900
    Magma_PREC_SST = 901
    Magma_PREC_HS = 902
    Magma_PREC_HST = 903
    Magma_PREC_SH = 904
    Magma_PREC_SHT = 905
    Magma_PREC_XHS_H = 910
    Magma_PREC_XHS_HTC = 911
    Magma_PREC_XHS_161616 = 912
    Magma_PREC_XHS_161616TC = 913
    Magma_PREC_XHS_161632TC = 914
    Magma_PREC_XSH_S = 915
    Magma_PREC_XSH_STC = 916
    Magma_PREC_XSH_163232TC = 917
    Magma_PREC_XSH_323232TC = 918
    Magma_REFINE_IRSTRS = 920
    Magma_REFINE_IRDTRS = 921
    Magma_REFINE_IRGMSTRS = 922
    Magma_REFINE_IRGMDTRS = 923
    Magma_REFINE_GMSTRS = 924
    Magma_REFINE_GMDTRS = 925
    Magma_REFINE_GMGMSTRS = 926
    Magma_REFINE_GMGMDTRS = 927
    Magma_PREC_HD = 930
end

@cenum magma_mp_type_t::UInt32 begin
    Magma_MP_BASE_SS = 950
    Magma_MP_BASE_DD = 951
    Magma_MP_BASE_XHS = 952
    Magma_MP_BASE_XSH = 953
    Magma_MP_BASE_XHD = 954
    Magma_MP_BASE_XDH = 955
    Magma_MP_ENABLE_DFLT_MATH = 960
    Magma_MP_ENABLE_TC_MATH = 961
    Magma_MP_SGEMM = 962
    Magma_MP_HGEMM = 963
    Magma_MP_GEMEX_I32_O32_C32 = 964
    Magma_MP_GEMEX_I16_O32_C32 = 965
    Magma_MP_GEMEX_I16_O16_C32 = 966
    Magma_MP_GEMEX_I16_O16_C16 = 967
    Magma_MP_TC_SGEMM = 968
    Magma_MP_TC_HGEMM = 969
    Magma_MP_TC_GEMEX_I32_O32_C32 = 970
    Magma_MP_TC_GEMEX_I16_O32_C32 = 971
    Magma_MP_TC_GEMEX_I16_O16_C32 = 972
    Magma_MP_TC_GEMEX_I16_O16_C16 = 973
end


struct zgehrd_data
    ngpu::magma_int_t
    ldda::magma_int_t
    ldv::magma_int_t
    ldvd::magma_int_t
    dA::NTuple{8, magmaDoubleComplex_ptr}
    dV::NTuple{8, magmaDoubleComplex_ptr}
    dVd::NTuple{8, magmaDoubleComplex_ptr}
    dY::NTuple{8, magmaDoubleComplex_ptr}
    dW::NTuple{8, magmaDoubleComplex_ptr}
    dTi::NTuple{8, magmaDoubleComplex_ptr}
    queues::NTuple{8, magma_queue_t}
end
#
# const magmablas_ctranspose_inplace = magmablas_ctranspose_inplace_v1
# const magmablas_ctranspose_conj_inplace = magmablas_ctranspose_conj_inplace_v1
# const magmablas_ctranspose = magmablas_ctranspose_v1
# const magmablas_ctranspose_conj = magmablas_ctranspose_conj_v1
# const magmablas_cgetmatrix_transpose = magmablas_cgetmatrix_transpose_v1
# const magmablas_csetmatrix_transpose = magmablas_csetmatrix_transpose_v1
# const magmablas_cprbt = magmablas_cprbt_v1
# const magmablas_cprbt_mv = magmablas_cprbt_mv_v1
# const magmablas_cprbt_mtv = magmablas_cprbt_mtv_v1
# const magma_cgetmatrix_1D_col_bcyclic = magma_cgetmatrix_1D_col_bcyclic_v1
# const magma_csetmatrix_1D_col_bcyclic = magma_csetmatrix_1D_col_bcyclic_v1
# const magma_cgetmatrix_1D_row_bcyclic = magma_cgetmatrix_1D_row_bcyclic_v1
# const magma_csetmatrix_1D_row_bcyclic = magma_csetmatrix_1D_row_bcyclic_v1
# const magmablas_cgeadd = magmablas_cgeadd_v1
# const magmablas_cgeadd2 = magmablas_cgeadd2_v1
# const magmablas_clacpy = magmablas_clacpy_v1
# const magmablas_clacpy_conj = magmablas_clacpy_conj_v1
# const magmablas_clacpy_sym_in = magmablas_clacpy_sym_in_v1
# const magmablas_clacpy_sym_out = magmablas_clacpy_sym_out_v1
# const magmablas_clange = magmablas_clange_v1
# const magmablas_clanhe = magmablas_clanhe_v1
# const magmablas_clansy = magmablas_clansy_v1
# const magmablas_clarfg = magmablas_clarfg_v1
# const magmablas_clascl = magmablas_clascl_v1
# const magmablas_clascl_2x2 = magmablas_clascl_2x2_v1
# const magmablas_clascl2 = magmablas_clascl2_v1
# const magmablas_clascl_diag = magmablas_clascl_diag_v1
# const magmablas_claset = magmablas_claset_v1
# const magmablas_claset_band = magmablas_claset_band_v1
# const magmablas_claswp = magmablas_claswp_v1
# const magmablas_claswp2 = magmablas_claswp2_v1
# const magmablas_claswp_sym = magmablas_claswp_sym_v1
# const magmablas_claswpx = magmablas_claswpx_v1
# const magmablas_csymmetrize = magmablas_csymmetrize_v1
# const magmablas_csymmetrize_tiles = magmablas_csymmetrize_tiles_v1
# const magmablas_ctrtri_diag = magmablas_ctrtri_diag_v1
# const magmablas_scnrm2_adjust = magmablas_scnrm2_adjust_v1
# const magmablas_scnrm2_check = magmablas_scnrm2_check_v1
# const magmablas_scnrm2_cols = magmablas_scnrm2_cols_v1
# const magmablas_scnrm2_row_check_adjust = magmablas_scnrm2_row_check_adjust_v1
# const magma_clarfb_gpu = magma_clarfb_gpu_v1
# const magma_clarfb_gpu_gemm = magma_clarfb_gpu_gemm_v1
# const magma_clarfbx_gpu = magma_clarfbx_gpu_v1
# const magma_clarfg_gpu = magma_clarfg_gpu_v1
# const magma_clarfgtx_gpu = magma_clarfgtx_gpu_v1
# const magma_clarfgx_gpu = magma_clarfgx_gpu_v1
# const magma_clarfx_gpu = magma_clarfx_gpu_v1
# const magmablas_caxpycp = magmablas_caxpycp_v1
# const magmablas_cswap = magmablas_cswap_v1
# const magmablas_cswapblk = magmablas_cswapblk_v1
# const magmablas_cswapdblk = magmablas_cswapdblk_v1
# const magmablas_cgemv = magmablas_cgemv_v1
# const magmablas_cgemv_conj = magmablas_cgemv_conj_v1
# const magmablas_chemv = magmablas_chemv_v1
# const magmablas_csymv = magmablas_csymv_v1
# const magmablas_cgemm = magmablas_cgemm_v1
# const magmablas_cgemm_reduce = magmablas_cgemm_reduce_v1
# const magmablas_chemm = magmablas_chemm_v1
# const magmablas_csymm = magmablas_csymm_v1
# const magmablas_csyr2k = magmablas_csyr2k_v1
# const magmablas_cher2k = magmablas_cher2k_v1
# const magmablas_csyrk = magmablas_csyrk_v1
# const magmablas_cherk = magmablas_cherk_v1
# const magmablas_ctrsm = magmablas_ctrsm_v1
# const magmablas_ctrsm_outofplace = magmablas_ctrsm_outofplace_v1
# const magmablas_ctrsm_work = magmablas_ctrsm_work_v1
# const magma_csetvector = magma_csetvector_v1
# const magma_cgetvector = magma_cgetvector_v1
# const magma_ccopyvector = magma_ccopyvector_v1
# const magma_csetmatrix = magma_csetmatrix_v1
# const magma_cgetmatrix = magma_cgetmatrix_v1
# const magma_ccopymatrix = magma_ccopymatrix_v1
# const magma_icamax = magma_icamax_v1
# const magma_icamin = magma_icamin_v1
# const magma_scasum = magma_scasum_v1
# const magma_caxpy = magma_caxpy_v1
# const magma_ccopy = magma_ccopy_v1
# const magma_cdotc = magma_cdotc_v1
# const magma_cdotu = magma_cdotu_v1
# const magma_scnrm2 = magma_scnrm2_v1
# const magma_crot = magma_crot_v1
# const magma_csrot = magma_csrot_v1
# const magma_crotm = magma_crotm_v1
# const magma_crotmg = magma_crotmg_v1
# const magma_cscal = magma_cscal_v1
# const magma_csscal = magma_csscal_v1
# const magma_cswap = magma_cswap_v1
# const magma_cgemv = magma_cgemv_v1
# const magma_cgerc = magma_cgerc_v1
# const magma_cgeru = magma_cgeru_v1
# const magma_chemv = magma_chemv_v1
# const magma_cher = magma_cher_v1
# const magma_cher2 = magma_cher2_v1
# const magma_ctrmv = magma_ctrmv_v1
# const magma_ctrsv = magma_ctrsv_v1
# const magma_cgemm = magma_cgemm_v1
# const magma_csymm = magma_csymm_v1
# const magma_chemm = magma_chemm_v1
# const magma_csyr2k = magma_csyr2k_v1
# const magma_cher2k = magma_cher2k_v1
# const magma_csyrk = magma_csyrk_v1
# const magma_cherk = magma_cherk_v1
# const magma_ctrmm = magma_ctrmm_v1
# const magma_ctrsm = magma_ctrsm_v1
#
# # Skipping MacroDefinition: magma_dsetvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_dsetvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dgetvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_dgetvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dcopyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_dcopyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dsetvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_dsetvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dgetvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_dgetvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dcopyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_dcopyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dsetmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_dsetmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dgetmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_dgetmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dcopymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_dcopymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dsetmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_dsetmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dgetmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_dgetmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dcopymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_dcopymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dsetvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_dsetvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dgetvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_dgetvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dcopyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_dcopyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dsetmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_dsetmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dgetmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_dgetmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_dcopymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_dcopymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
#
# const magmablas_dtranspose_inplace = magmablas_dtranspose_inplace_v1
# const magmablas_dtranspose = magmablas_dtranspose_v1
# const magmablas_dgetmatrix_transpose = magmablas_dgetmatrix_transpose_v1
# const magmablas_dsetmatrix_transpose = magmablas_dsetmatrix_transpose_v1
# const magmablas_dprbt = magmablas_dprbt_v1
# const magmablas_dprbt_mv = magmablas_dprbt_mv_v1
# const magmablas_dprbt_mtv = magmablas_dprbt_mtv_v1
# const magma_dgetmatrix_1D_col_bcyclic = magma_dgetmatrix_1D_col_bcyclic_v1
# const magma_dsetmatrix_1D_col_bcyclic = magma_dsetmatrix_1D_col_bcyclic_v1
# const magma_dgetmatrix_1D_row_bcyclic = magma_dgetmatrix_1D_row_bcyclic_v1
# const magma_dsetmatrix_1D_row_bcyclic = magma_dsetmatrix_1D_row_bcyclic_v1
# const magmablas_dgeadd = magmablas_dgeadd_v1
# const magmablas_dgeadd2 = magmablas_dgeadd2_v1
# const magmablas_dlacpy = magmablas_dlacpy_v1
# const magmablas_dlacpy_conj = magmablas_dlacpy_conj_v1
# const magmablas_dlacpy_sym_in = magmablas_dlacpy_sym_in_v1
# const magmablas_dlacpy_sym_out = magmablas_dlacpy_sym_out_v1
# const magmablas_dlange = magmablas_dlange_v1
# const magmablas_dlansy = magmablas_dlansy_v1
# const magmablas_dlarfg = magmablas_dlarfg_v1
# const magmablas_dlascl = magmablas_dlascl_v1
# const magmablas_dlascl_2x2 = magmablas_dlascl_2x2_v1
# const magmablas_dlascl2 = magmablas_dlascl2_v1
# const magmablas_dlascl_diag = magmablas_dlascl_diag_v1
# const magmablas_dlaset = magmablas_dlaset_v1
# const magmablas_dlaset_band = magmablas_dlaset_band_v1
# const magmablas_dlaswp = magmablas_dlaswp_v1
# const magmablas_dlaswp2 = magmablas_dlaswp2_v1
# const magmablas_dlaswp_sym = magmablas_dlaswp_sym_v1
# const magmablas_dlaswpx = magmablas_dlaswpx_v1
# const magmablas_dsymmetrize = magmablas_dsymmetrize_v1
# const magmablas_dsymmetrize_tiles = magmablas_dsymmetrize_tiles_v1
# const magmablas_dtrtri_diag = magmablas_dtrtri_diag_v1
# const magmablas_dnrm2_adjust = magmablas_dnrm2_adjust_v1
# const magmablas_dnrm2_check = magmablas_dnrm2_check_v1
# const magmablas_dnrm2_cols = magmablas_dnrm2_cols_v1
# const magmablas_dnrm2_row_check_adjust = magmablas_dnrm2_row_check_adjust_v1
# const magma_dlarfb_gpu = magma_dlarfb_gpu_v1
# const magma_dlarfb_gpu_gemm = magma_dlarfb_gpu_gemm_v1
# const magma_dlarfbx_gpu = magma_dlarfbx_gpu_v1
# const magma_dlarfg_gpu = magma_dlarfg_gpu_v1
# const magma_dlarfgtx_gpu = magma_dlarfgtx_gpu_v1
# const magma_dlarfgx_gpu = magma_dlarfgx_gpu_v1
# const magma_dlarfx_gpu = magma_dlarfx_gpu_v1
# const magmablas_daxpycp = magmablas_daxpycp_v1
# const magmablas_dswap = magmablas_dswap_v1
# const magmablas_dswapblk = magmablas_dswapblk_v1
# const magmablas_dswapdblk = magmablas_dswapdblk_v1
# const magmablas_dgemv = magmablas_dgemv_v1
# const magmablas_dgemv_conj = magmablas_dgemv_conj_v1
# const magmablas_dsymv = magmablas_dsymv_v1
# const magmablas_dgemm = magmablas_dgemm_v1
# const magmablas_dgemm_reduce = magmablas_dgemm_reduce_v1
# const magmablas_dsymm = magmablas_dsymm_v1
# const magmablas_dsyr2k = magmablas_dsyr2k_v1
# const magmablas_dsyrk = magmablas_dsyrk_v1
# const magmablas_dtrsm = magmablas_dtrsm_v1
# const magmablas_dtrsm_outofplace = magmablas_dtrsm_outofplace_v1
# const magmablas_dtrsm_work = magmablas_dtrsm_work_v1
# const magma_dsetvector = magma_dsetvector_v1
# const magma_dgetvector = magma_dgetvector_v1
# const magma_dcopyvector = magma_dcopyvector_v1
# const magma_dsetmatrix = magma_dsetmatrix_v1
# const magma_dgetmatrix = magma_dgetmatrix_v1
# const magma_dcopymatrix = magma_dcopymatrix_v1
# const magma_idamax = magma_idamax_v1
# const magma_idamin = magma_idamin_v1
# const magma_dasum = magma_dasum_v1
# const magma_daxpy = magma_daxpy_v1
# const magma_dcopy = magma_dcopy_v1
# const magma_ddot = magma_ddot_v1
# const magma_dnrm2 = magma_dnrm2_v1
# const magma_drot = magma_drot_v1
# const magma_drotm = magma_drotm_v1
# const magma_drotmg = magma_drotmg_v1
# const magma_dscal = magma_dscal_v1
# const magma_dswap = magma_dswap_v1
# const magma_dgemv = magma_dgemv_v1
# const magma_dger = magma_dger_v1
# const magma_dsymv = magma_dsymv_v1
# const magma_dsyr = magma_dsyr_v1
# const magma_dsyr2 = magma_dsyr2_v1
# const magma_dtrmv = magma_dtrmv_v1
# const magma_dtrsv = magma_dtrsv_v1
# const magma_dgemm = magma_dgemm_v1
# const magma_dsymm = magma_dsymm_v1
# const magma_dsyr2k = magma_dsyr2k_v1
# const magma_dsyrk = magma_dsyrk_v1
# const magma_dtrmm = magma_dtrmm_v1
# const magma_dtrsm = magma_dtrsm_v1
# const magmablas_dsaxpycp = magmablas_dsaxpycp_v1
# const magmablas_dslaswp = magmablas_dslaswp_v1
# const magmablas_dlag2s = magmablas_dlag2s_v1
# const magmablas_slag2d = magmablas_slag2d_v1
# const magmablas_dlat2s = magmablas_dlat2s_v1
# const magmablas_slat2d = magmablas_slat2d_v1
#
# # Skipping MacroDefinition: magma_ssetvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_ssetvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_sgetvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_sgetvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_scopyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_scopyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_ssetvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_ssetvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_sgetvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_sgetvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_scopyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_scopyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_ssetmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_ssetmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_sgetmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_sgetmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_scopymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_scopymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_ssetmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_ssetmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_sgetmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_sgetmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_scopymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_scopymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_ssetvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_ssetvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_sgetvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_sgetvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_scopyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_scopyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_ssetmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_ssetmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_sgetmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_sgetmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_scopymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_scopymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
#
# const magmablas_stranspose_inplace = magmablas_stranspose_inplace_v1
# const magmablas_stranspose = magmablas_stranspose_v1
# const magmablas_sgetmatrix_transpose = magmablas_sgetmatrix_transpose_v1
# const magmablas_ssetmatrix_transpose = magmablas_ssetmatrix_transpose_v1
# const magmablas_sprbt = magmablas_sprbt_v1
# const magmablas_sprbt_mv = magmablas_sprbt_mv_v1
# const magmablas_sprbt_mtv = magmablas_sprbt_mtv_v1
# const magma_sgetmatrix_1D_col_bcyclic = magma_sgetmatrix_1D_col_bcyclic_v1
# const magma_ssetmatrix_1D_col_bcyclic = magma_ssetmatrix_1D_col_bcyclic_v1
# const magma_sgetmatrix_1D_row_bcyclic = magma_sgetmatrix_1D_row_bcyclic_v1
# const magma_ssetmatrix_1D_row_bcyclic = magma_ssetmatrix_1D_row_bcyclic_v1
# const magmablas_sgeadd = magmablas_sgeadd_v1
# const magmablas_sgeadd2 = magmablas_sgeadd2_v1
# const magmablas_slacpy = magmablas_slacpy_v1
# const magmablas_slacpy_conj = magmablas_slacpy_conj_v1
# const magmablas_slacpy_sym_in = magmablas_slacpy_sym_in_v1
# const magmablas_slacpy_sym_out = magmablas_slacpy_sym_out_v1
# const magmablas_slange = magmablas_slange_v1
# const magmablas_slansy = magmablas_slansy_v1
# const magmablas_slarfg = magmablas_slarfg_v1
# const magmablas_slascl = magmablas_slascl_v1
# const magmablas_slascl_2x2 = magmablas_slascl_2x2_v1
# const magmablas_slascl2 = magmablas_slascl2_v1
# const magmablas_slascl_diag = magmablas_slascl_diag_v1
# const magmablas_slaset = magmablas_slaset_v1
# const magmablas_slaset_band = magmablas_slaset_band_v1
# const magmablas_slaswp = magmablas_slaswp_v1
# const magmablas_slaswp2 = magmablas_slaswp2_v1
# const magmablas_slaswp_sym = magmablas_slaswp_sym_v1
# const magmablas_slaswpx = magmablas_slaswpx_v1
# const magmablas_ssymmetrize = magmablas_ssymmetrize_v1
# const magmablas_ssymmetrize_tiles = magmablas_ssymmetrize_tiles_v1
# const magmablas_strtri_diag = magmablas_strtri_diag_v1
# const magmablas_snrm2_adjust = magmablas_snrm2_adjust_v1
# const magmablas_snrm2_check = magmablas_snrm2_check_v1
# const magmablas_snrm2_cols = magmablas_snrm2_cols_v1
# const magmablas_snrm2_row_check_adjust = magmablas_snrm2_row_check_adjust_v1
# const magma_slarfb_gpu = magma_slarfb_gpu_v1
# const magma_slarfb_gpu_gemm = magma_slarfb_gpu_gemm_v1
# const magma_slarfbx_gpu = magma_slarfbx_gpu_v1
# const magma_slarfg_gpu = magma_slarfg_gpu_v1
# const magma_slarfgtx_gpu = magma_slarfgtx_gpu_v1
# const magma_slarfgx_gpu = magma_slarfgx_gpu_v1
# const magma_slarfx_gpu = magma_slarfx_gpu_v1
# const magmablas_saxpycp = magmablas_saxpycp_v1
# const magmablas_sswap = magmablas_sswap_v1
# const magmablas_sswapblk = magmablas_sswapblk_v1
# const magmablas_sswapdblk = magmablas_sswapdblk_v1
# const magmablas_sgemv = magmablas_sgemv_v1
# const magmablas_sgemv_conj = magmablas_sgemv_conj_v1
# const magmablas_ssymv = magmablas_ssymv_v1
# const magmablas_sgemm = magmablas_sgemm_v1
# const magmablas_sgemm_reduce = magmablas_sgemm_reduce_v1
# const magmablas_ssymm = magmablas_ssymm_v1
# const magmablas_ssyr2k = magmablas_ssyr2k_v1
# const magmablas_ssyrk = magmablas_ssyrk_v1
# const magmablas_strsm = magmablas_strsm_v1
# const magmablas_strsm_outofplace = magmablas_strsm_outofplace_v1
# const magmablas_strsm_work = magmablas_strsm_work_v1
# const magma_ssetvector = magma_ssetvector_v1
# const magma_sgetvector = magma_sgetvector_v1
# const magma_scopyvector = magma_scopyvector_v1
# const magma_ssetmatrix = magma_ssetmatrix_v1
# const magma_sgetmatrix = magma_sgetmatrix_v1
# const magma_scopymatrix = magma_scopymatrix_v1
# const magma_isamax = magma_isamax_v1
# const magma_isamin = magma_isamin_v1
# const magma_sasum = magma_sasum_v1
# const magma_saxpy = magma_saxpy_v1
# const magma_scopy = magma_scopy_v1
# const magma_sdot = magma_sdot_v1
# const magma_snrm2 = magma_snrm2_v1
# const magma_srot = magma_srot_v1
# const magma_srotm = magma_srotm_v1
# const magma_srotmg = magma_srotmg_v1
# const magma_sscal = magma_sscal_v1
# const magma_sswap = magma_sswap_v1
# const magma_sgemv = magma_sgemv_v1
# const magma_sger = magma_sger_v1
# const magma_ssymv = magma_ssymv_v1
# const magma_ssyr = magma_ssyr_v1
# const magma_ssyr2 = magma_ssyr2_v1
# const magma_strmv = magma_strmv_v1
# const magma_strsv = magma_strsv_v1
# const magma_sgemm = magma_sgemm_v1
# const magma_ssymm = magma_ssymm_v1
# const magma_ssyr2k = magma_ssyr2k_v1
# const magma_ssyrk = magma_ssyrk_v1
# const magma_strmm = magma_strmm_v1
# const magma_strsm = magma_strsm_v1
#
# # Skipping MacroDefinition: magma_queue_create_v1 ( queue_ptr ) magma_queue_create_v1_internal ( queue_ptr , __func__ , __FILE__ , __LINE__ )
#
# const MagmaUpperLower = MagmaFull
# const MagmaUpperLowerStr = MagmaFullStr
#
# # Skipping MacroDefinition: MAGMA_Z_CNJG ( a ) MAGMA_Z_CONJ ( a )
# # Skipping MacroDefinition: MAGMA_C_CNJG ( a ) MAGMA_C_CONJ ( a )
# # Skipping MacroDefinition: MAGMA_D_CNJG ( a ) MAGMA_D_CONJ ( a )
# # Skipping MacroDefinition: MAGMA_S_CNJG ( a ) MAGMA_S_CONJ ( a )
#
# const magma_queue_create = magma_queue_create_v1
# const magma_setvector = magma_setvector_v1
# const magma_getvector = magma_getvector_v1
# const magma_copyvector = magma_copyvector_v1
# const magma_setmatrix = magma_setmatrix_v1
# const magma_getmatrix = magma_getmatrix_v1
# const magma_copymatrix = magma_copymatrix_v1
# const magma_isetvector = magma_isetvector_v1
# const magma_igetvector = magma_igetvector_v1
# const magma_icopyvector = magma_icopyvector_v1
# const magma_isetmatrix = magma_isetmatrix_v1
# const magma_igetmatrix = magma_igetmatrix_v1
# const magma_icopymatrix = magma_icopymatrix_v1
# const magma_index_setvector = magma_index_setvector_v1
# const magma_index_getvector = magma_index_getvector_v1
# const magma_index_copyvector = magma_index_copyvector_v1
# const magma_index_setmatrix = magma_index_setmatrix_v1
# const magma_index_getmatrix = magma_index_getmatrix_v1
# const magma_index_copymatrix = magma_index_copymatrix_v1
#
# # Skipping MacroDefinition: magma_zsetvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_zsetvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zgetvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_zgetvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zcopyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_zcopyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zsetvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_zsetvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zgetvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_zgetvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zcopyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_zcopyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zsetmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_zsetmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zgetmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_zgetmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zcopymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_zcopymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zsetmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_zsetmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zgetmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_zgetmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zcopymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_zcopymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zsetvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_zsetvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zgetvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_zgetvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zcopyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_zcopyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zsetmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_zsetmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zgetmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_zgetmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# # Skipping MacroDefinition: magma_zcopymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_zcopymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
#
# const magmablas_ztranspose_inplace = magmablas_ztranspose_inplace_v1
# const magmablas_ztranspose_conj_inplace = magmablas_ztranspose_conj_inplace_v1
# const magmablas_ztranspose = magmablas_ztranspose_v1
# const magmablas_ztranspose_conj = magmablas_ztranspose_conj_v1
# const magmablas_zgetmatrix_transpose = magmablas_zgetmatrix_transpose_v1
# const magmablas_zsetmatrix_transpose = magmablas_zsetmatrix_transpose_v1
# const magmablas_zprbt = magmablas_zprbt_v1
# const magmablas_zprbt_mv = magmablas_zprbt_mv_v1
# const magmablas_zprbt_mtv = magmablas_zprbt_mtv_v1
# const magma_zgetmatrix_1D_col_bcyclic = magma_zgetmatrix_1D_col_bcyclic_v1
# const magma_zsetmatrix_1D_col_bcyclic = magma_zsetmatrix_1D_col_bcyclic_v1
# const magma_zgetmatrix_1D_row_bcyclic = magma_zgetmatrix_1D_row_bcyclic_v1
# const magma_zsetmatrix_1D_row_bcyclic = magma_zsetmatrix_1D_row_bcyclic_v1
# const magmablas_zgeadd = magmablas_zgeadd_v1
# const magmablas_zgeadd2 = magmablas_zgeadd2_v1
# const magmablas_zlacpy = magmablas_zlacpy_v1
# const magmablas_zlacpy_conj = magmablas_zlacpy_conj_v1
# const magmablas_zlacpy_sym_in = magmablas_zlacpy_sym_in_v1
# const magmablas_zlacpy_sym_out = magmablas_zlacpy_sym_out_v1
# const magmablas_zlange = magmablas_zlange_v1
# const magmablas_zlanhe = magmablas_zlanhe_v1
# const magmablas_zlansy = magmablas_zlansy_v1
# const magmablas_zlarfg = magmablas_zlarfg_v1
# const magmablas_zlascl = magmablas_zlascl_v1
# const magmablas_zlascl_2x2 = magmablas_zlascl_2x2_v1
# const magmablas_zlascl2 = magmablas_zlascl2_v1
# const magmablas_zlascl_diag = magmablas_zlascl_diag_v1
# const magmablas_zlaset = magmablas_zlaset_v1
# const magmablas_zlaset_band = magmablas_zlaset_band_v1
# const magmablas_zlaswp = magmablas_zlaswp_v1
# const magmablas_zlaswp2 = magmablas_zlaswp2_v1
# const magmablas_zlaswp_sym = magmablas_zlaswp_sym_v1
# const magmablas_zlaswpx = magmablas_zlaswpx_v1
# const magmablas_zsymmetrize = magmablas_zsymmetrize_v1
# const magmablas_zsymmetrize_tiles = magmablas_zsymmetrize_tiles_v1
# const magmablas_ztrtri_diag = magmablas_ztrtri_diag_v1
# const magmablas_dznrm2_adjust = magmablas_dznrm2_adjust_v1
# const magmablas_dznrm2_check = magmablas_dznrm2_check_v1
# const magmablas_dznrm2_cols = magmablas_dznrm2_cols_v1
# const magmablas_dznrm2_row_check_adjust = magmablas_dznrm2_row_check_adjust_v1
# const magma_zlarfb_gpu = magma_zlarfb_gpu_v1
# const magma_zlarfb_gpu_gemm = magma_zlarfb_gpu_gemm_v1
# const magma_zlarfbx_gpu = magma_zlarfbx_gpu_v1
# const magma_zlarfg_gpu = magma_zlarfg_gpu_v1
# const magma_zlarfgtx_gpu = magma_zlarfgtx_gpu_v1
# const magma_zlarfgx_gpu = magma_zlarfgx_gpu_v1
# const magma_zlarfx_gpu = magma_zlarfx_gpu_v1
# const magmablas_zaxpycp = magmablas_zaxpycp_v1
# const magmablas_zswap = magmablas_zswap_v1
# const magmablas_zswapblk = magmablas_zswapblk_v1
# const magmablas_zswapdblk = magmablas_zswapdblk_v1
# const magmablas_zgemv = magmablas_zgemv_v1
# const magmablas_zgemv_conj = magmablas_zgemv_conj_v1
# const magmablas_zhemv = magmablas_zhemv_v1
# const magmablas_zsymv = magmablas_zsymv_v1
# const magmablas_zgemm = magmablas_zgemm_v1
# const magmablas_zgemm_reduce = magmablas_zgemm_reduce_v1
# const magmablas_zhemm = magmablas_zhemm_v1
# const magmablas_zsymm = magmablas_zsymm_v1
# const magmablas_zsyr2k = magmablas_zsyr2k_v1
# const magmablas_zher2k = magmablas_zher2k_v1
# const magmablas_zsyrk = magmablas_zsyrk_v1
# const magmablas_zherk = magmablas_zherk_v1
# const magmablas_ztrsm = magmablas_ztrsm_v1
# const magmablas_ztrsm_outofplace = magmablas_ztrsm_outofplace_v1
# const magmablas_ztrsm_work = magmablas_ztrsm_work_v1
# const magma_zsetvector = magma_zsetvector_v1
# const magma_zgetvector = magma_zgetvector_v1
# const magma_zcopyvector = magma_zcopyvector_v1
# const magma_zsetmatrix = magma_zsetmatrix_v1
# const magma_zgetmatrix = magma_zgetmatrix_v1
# const magma_zcopymatrix = magma_zcopymatrix_v1
# const magma_izamax = magma_izamax_v1
# const magma_izamin = magma_izamin_v1
# const magma_dzasum = magma_dzasum_v1
# const magma_zaxpy = magma_zaxpy_v1
# const magma_zcopy = magma_zcopy_v1
# const magma_zdotc = magma_zdotc_v1
# const magma_zdotu = magma_zdotu_v1
# const magma_dznrm2 = magma_dznrm2_v1
# const magma_zrot = magma_zrot_v1
# const magma_zdrot = magma_zdrot_v1
# const magma_zrotm = magma_zrotm_v1
# const magma_zrotmg = magma_zrotmg_v1
# const magma_zscal = magma_zscal_v1
# const magma_zdscal = magma_zdscal_v1
# const magma_zswap = magma_zswap_v1
# const magma_zgemv = magma_zgemv_v1
# const magma_zgerc = magma_zgerc_v1
# const magma_zgeru = magma_zgeru_v1
# const magma_zhemv = magma_zhemv_v1
# const magma_zher = magma_zher_v1
# const magma_zher2 = magma_zher2_v1
# const magma_ztrmv = magma_ztrmv_v1
# const magma_ztrsv = magma_ztrsv_v1
# const magma_zgemm = magma_zgemm_v1
# const magma_zsymm = magma_zsymm_v1
# const magma_zhemm = magma_zhemm_v1
# const magma_zsyr2k = magma_zsyr2k_v1
# const magma_zher2k = magma_zher2k_v1
# const magma_zsyrk = magma_zsyrk_v1
# const magma_zherk = magma_zherk_v1
# const magma_ztrmm = magma_ztrmm_v1
# const magma_ztrsm = magma_ztrsm_v1
# const magmablas_zcaxpycp = magmablas_zcaxpycp_v1
# const magmablas_zclaswp = magmablas_zclaswp_v1
# const magmablas_zlag2c = magmablas_zlag2c_v1
# const magmablas_clag2z = magmablas_clag2z_v1
# const magmablas_zlat2c = magmablas_zlat2c_v1
# const magmablas_clat2z = magmablas_clat2z_v1
