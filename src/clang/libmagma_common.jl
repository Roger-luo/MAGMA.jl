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
    dQ1::Ptr{magmaFloatComplex}
    dT1::Ptr{magmaFloatComplex}
    dT2::Ptr{magmaFloatComplex}
    dV2::Ptr{magmaFloatComplex}
    dE::Ptr{magmaFloatComplex}
    T::Ptr{magmaFloatComplex}
    A::Ptr{magmaFloatComplex}
    V::Ptr{magmaFloatComplex}
    TAU::Ptr{magmaFloatComplex}
    E::Ptr{magmaFloatComplex}
    E_CPU::Ptr{magmaFloatComplex}
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
    timeblg::Ptr{real_Double_t}
    timeaplQ::Ptr{real_Double_t}
    ss_prog::Ptr{Cint}
end

const magma_int_t = Cint
const magmaFloatComplex_ptr = Ptr{magmaFloatComplex}
const magma_queue = Cvoid
const magma_queue_t = Ptr{magma_queue}

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

# Skipping MacroDefinition: blasf77_icamax FORTRAN_NAME ( icamax , ICAMAX )
# Skipping MacroDefinition: blasf77_caxpy FORTRAN_NAME ( caxpy , CAXPY )
# Skipping MacroDefinition: blasf77_ccopy FORTRAN_NAME ( ccopy , CCOPY )
# Skipping MacroDefinition: blasf77_cgemm FORTRAN_NAME ( cgemm , CGEMM )
# Skipping MacroDefinition: blasf77_cgemv FORTRAN_NAME ( cgemv , CGEMV )
# Skipping MacroDefinition: blasf77_cgerc FORTRAN_NAME ( cgerc , CGERC )
# Skipping MacroDefinition: blasf77_cgeru FORTRAN_NAME ( cgeru , CGERU )
# Skipping MacroDefinition: blasf77_chemm FORTRAN_NAME ( chemm , CHEMM )
# Skipping MacroDefinition: blasf77_chemv FORTRAN_NAME ( chemv , CHEMV )
# Skipping MacroDefinition: blasf77_cher FORTRAN_NAME ( cher , CHER )
# Skipping MacroDefinition: blasf77_cher2 FORTRAN_NAME ( cher2 , CHER2 )
# Skipping MacroDefinition: blasf77_cher2k FORTRAN_NAME ( cher2k , CHER2K )
# Skipping MacroDefinition: blasf77_cherk FORTRAN_NAME ( cherk , CHERK )
# Skipping MacroDefinition: blasf77_cscal FORTRAN_NAME ( cscal , CSCAL )
# Skipping MacroDefinition: blasf77_csscal FORTRAN_NAME ( csscal , CSSCAL )
# Skipping MacroDefinition: blasf77_cswap FORTRAN_NAME ( cswap , CSWAP )
# Skipping MacroDefinition: blasf77_csymm FORTRAN_NAME ( csymm , CSYMM )
# Skipping MacroDefinition: blasf77_csyr2k FORTRAN_NAME ( csyr2k , CSYR2K )
# Skipping MacroDefinition: blasf77_csyrk FORTRAN_NAME ( csyrk , CSYRK )
# Skipping MacroDefinition: blasf77_crotg FORTRAN_NAME ( crotg , CROTG )
# Skipping MacroDefinition: blasf77_crot FORTRAN_NAME ( crot , CROT )
# Skipping MacroDefinition: blasf77_csrot FORTRAN_NAME ( csrot , CSROT )
# Skipping MacroDefinition: blasf77_ctrmm FORTRAN_NAME ( ctrmm , CTRMM )
# Skipping MacroDefinition: blasf77_ctrmv FORTRAN_NAME ( ctrmv , CTRMV )
# Skipping MacroDefinition: blasf77_ctrsm FORTRAN_NAME ( ctrsm , CTRSM )
# Skipping MacroDefinition: blasf77_ctrsv FORTRAN_NAME ( ctrsv , CTRSV )
# Skipping MacroDefinition: lapackf77_slaed2 FORTRAN_NAME ( slaed2 , SLAED2 )
# Skipping MacroDefinition: lapackf77_slaed4 FORTRAN_NAME ( slaed4 , SLAED4 )
# Skipping MacroDefinition: lapackf77_slaln2 FORTRAN_NAME ( slaln2 , SLALN2 )
# Skipping MacroDefinition: lapackf77_slamc3 FORTRAN_NAME ( slamc3 , SLAMC3 )
# Skipping MacroDefinition: lapackf77_slamrg FORTRAN_NAME ( slamrg , SLAMRG )
# Skipping MacroDefinition: lapackf77_slasrt FORTRAN_NAME ( slasrt , SLASRT )
# Skipping MacroDefinition: lapackf77_sstebz FORTRAN_NAME ( sstebz , SSTEBZ )
# Skipping MacroDefinition: lapackf77_sbdsdc FORTRAN_NAME ( sbdsdc , SBDSDC )
# Skipping MacroDefinition: lapackf77_cbdsqr FORTRAN_NAME ( cbdsqr , CBDSQR )
# Skipping MacroDefinition: lapackf77_cgebak FORTRAN_NAME ( cgebak , CGEBAK )
# Skipping MacroDefinition: lapackf77_cgebal FORTRAN_NAME ( cgebal , CGEBAL )
# Skipping MacroDefinition: lapackf77_cgebd2 FORTRAN_NAME ( cgebd2 , CGEBD2 )
# Skipping MacroDefinition: lapackf77_cgebrd FORTRAN_NAME ( cgebrd , CGEBRD )
# Skipping MacroDefinition: lapackf77_cgbbrd FORTRAN_NAME ( cgbbrd , CGBBRD )
# Skipping MacroDefinition: lapackf77_cgbsv FORTRAN_NAME ( cgbsv , CGBSV )
# Skipping MacroDefinition: lapackf77_cgeev FORTRAN_NAME ( cgeev , CGEEV )
# Skipping MacroDefinition: lapackf77_cgehd2 FORTRAN_NAME ( cgehd2 , CGEHD2 )
# Skipping MacroDefinition: lapackf77_cgehrd FORTRAN_NAME ( cgehrd , CGEHRD )
# Skipping MacroDefinition: lapackf77_cgelqf FORTRAN_NAME ( cgelqf , CGELQF )
# Skipping MacroDefinition: lapackf77_cgels FORTRAN_NAME ( cgels , CGELS )
# Skipping MacroDefinition: lapackf77_cgeqlf FORTRAN_NAME ( cgeqlf , CGEQLF )
# Skipping MacroDefinition: lapackf77_cgeqp3 FORTRAN_NAME ( cgeqp3 , CGEQP3 )
# Skipping MacroDefinition: lapackf77_cgeqrf FORTRAN_NAME ( cgeqrf , CGEQRF )
# Skipping MacroDefinition: lapackf77_cgerqf FORTRAN_NAME ( cgerqf , CGERQF )
# Skipping MacroDefinition: lapackf77_cgesdd FORTRAN_NAME ( cgesdd , CGESDD )
# Skipping MacroDefinition: lapackf77_cgesv FORTRAN_NAME ( cgesv , CGESV )
# Skipping MacroDefinition: lapackf77_cgesvd FORTRAN_NAME ( cgesvd , CGESVD )
# Skipping MacroDefinition: lapackf77_cgetrf FORTRAN_NAME ( cgetrf , CGETRF )
# Skipping MacroDefinition: lapackf77_cgetri FORTRAN_NAME ( cgetri , CGETRI )
# Skipping MacroDefinition: lapackf77_cgetrs FORTRAN_NAME ( cgetrs , CGETRS )
# Skipping MacroDefinition: lapackf77_cgglse FORTRAN_NAME ( cgglse , CGGLSE )
# Skipping MacroDefinition: lapackf77_cggrqf FORTRAN_NAME ( cggrqf , CGGRQF )
# Skipping MacroDefinition: lapackf77_chetf2 FORTRAN_NAME ( chetf2 , CHETF2 )
# Skipping MacroDefinition: lapackf77_chetrs FORTRAN_NAME ( chetrs , CHETRS )
# Skipping MacroDefinition: lapackf77_chbtrd FORTRAN_NAME ( chbtrd , CHBTRD )
# Skipping MacroDefinition: lapackf77_cheev FORTRAN_NAME ( cheev , CHEEV )
# Skipping MacroDefinition: lapackf77_cheevd FORTRAN_NAME ( cheevd , CHEEVD )
# Skipping MacroDefinition: lapackf77_cheevr FORTRAN_NAME ( cheevr , CHEEVR )
# Skipping MacroDefinition: lapackf77_cheevx FORTRAN_NAME ( cheevx , CHEEVX )
# Skipping MacroDefinition: lapackf77_chegs2 FORTRAN_NAME ( chegs2 , CHEGS2 )
# Skipping MacroDefinition: lapackf77_chegst FORTRAN_NAME ( chegst , CHEGST )
# Skipping MacroDefinition: lapackf77_chegvd FORTRAN_NAME ( chegvd , CHEGVD )
# Skipping MacroDefinition: lapackf77_chetd2 FORTRAN_NAME ( chetd2 , CHETD2 )
# Skipping MacroDefinition: lapackf77_chetrd FORTRAN_NAME ( chetrd , CHETRD )
# Skipping MacroDefinition: lapackf77_chetrf FORTRAN_NAME ( chetrf , CHETRF )
# Skipping MacroDefinition: lapackf77_chesv FORTRAN_NAME ( chesv , CHESV )
# Skipping MacroDefinition: lapackf77_chseqr FORTRAN_NAME ( chseqr , CHSEQR )
# Skipping MacroDefinition: lapackf77_clabrd FORTRAN_NAME ( clabrd , CLABRD )
# Skipping MacroDefinition: lapackf77_clacgv FORTRAN_NAME ( clacgv , CLACGV )
# Skipping MacroDefinition: lapackf77_clacp2 FORTRAN_NAME ( clacp2 , CLACP2 )
# Skipping MacroDefinition: lapackf77_clacpy FORTRAN_NAME ( clacpy , CLACPY )
# Skipping MacroDefinition: lapackf77_clacrm FORTRAN_NAME ( clacrm , CLACRM )
# Skipping MacroDefinition: lapackf77_cladiv FORTRAN_NAME ( cladiv , CLADIV )
# Skipping MacroDefinition: lapackf77_clahef FORTRAN_NAME ( clahef , CLAHEF )
# Skipping MacroDefinition: lapackf77_clange FORTRAN_NAME ( clange , CLANGE )
# Skipping MacroDefinition: lapackf77_clanhe FORTRAN_NAME ( clanhe , CLANHE )
# Skipping MacroDefinition: lapackf77_clanht FORTRAN_NAME ( clanht , CLANHT )
# Skipping MacroDefinition: lapackf77_clansy FORTRAN_NAME ( clansy , CLANSY )
# Skipping MacroDefinition: lapackf77_clantr FORTRAN_NAME ( clantr , CLANTR )
# Skipping MacroDefinition: lapackf77_slapy3 FORTRAN_NAME ( slapy3 , SLAPY3 )
# Skipping MacroDefinition: lapackf77_claqp2 FORTRAN_NAME ( claqp2 , CLAQP2 )
# Skipping MacroDefinition: lapackf77_clarcm FORTRAN_NAME ( clarcm , CLARCM )
# Skipping MacroDefinition: lapackf77_clarf FORTRAN_NAME ( clarf , CLARF )
# Skipping MacroDefinition: lapackf77_clarfb FORTRAN_NAME ( clarfb , CLARFB )
# Skipping MacroDefinition: lapackf77_clarfg FORTRAN_NAME ( clarfg , CLARFG )
# Skipping MacroDefinition: lapackf77_clarft FORTRAN_NAME ( clarft , CLARFT )
# Skipping MacroDefinition: lapackf77_clarfx FORTRAN_NAME ( clarfx , CLARFX )
# Skipping MacroDefinition: lapackf77_clarnv FORTRAN_NAME ( clarnv , CLARNV )
# Skipping MacroDefinition: lapackf77_clartg FORTRAN_NAME ( clartg , CLARTG )
# Skipping MacroDefinition: lapackf77_clascl FORTRAN_NAME ( clascl , CLASCL )
# Skipping MacroDefinition: lapackf77_claset FORTRAN_NAME ( claset , CLASET )
# Skipping MacroDefinition: lapackf77_claswp FORTRAN_NAME ( claswp , CLASWP )
# Skipping MacroDefinition: lapackf77_clatrd FORTRAN_NAME ( clatrd , CLATRD )
# Skipping MacroDefinition: lapackf77_clatrs FORTRAN_NAME ( clatrs , CLATRS )
# Skipping MacroDefinition: lapackf77_clauum FORTRAN_NAME ( clauum , CLAUUM )
# Skipping MacroDefinition: lapackf77_clavhe FORTRAN_NAME ( clavhe , CLAVHE )
# Skipping MacroDefinition: lapackf77_cposv FORTRAN_NAME ( cposv , CPOSV )
# Skipping MacroDefinition: lapackf77_cpotrf FORTRAN_NAME ( cpotrf , CPOTRF )
# Skipping MacroDefinition: lapackf77_cpotri FORTRAN_NAME ( cpotri , CPOTRI )
# Skipping MacroDefinition: lapackf77_cpotrs FORTRAN_NAME ( cpotrs , CPOTRS )
# Skipping MacroDefinition: lapackf77_cstedc FORTRAN_NAME ( cstedc , CSTEDC )
# Skipping MacroDefinition: lapackf77_cstein FORTRAN_NAME ( cstein , CSTEIN )
# Skipping MacroDefinition: lapackf77_cstemr FORTRAN_NAME ( cstemr , CSTEMR )
# Skipping MacroDefinition: lapackf77_csteqr FORTRAN_NAME ( csteqr , CSTEQR )
# Skipping MacroDefinition: lapackf77_csymv FORTRAN_NAME ( csymv , CSYMV )
# Skipping MacroDefinition: lapackf77_csyr FORTRAN_NAME ( csyr , CSYR )
# Skipping MacroDefinition: lapackf77_csysv FORTRAN_NAME ( csysv , CSYSV )
# Skipping MacroDefinition: lapackf77_ctrevc FORTRAN_NAME ( ctrevc , CTREVC )
# Skipping MacroDefinition: lapackf77_ctrevc3 FORTRAN_NAME ( ctrevc3 , CTREVC3 )
# Skipping MacroDefinition: lapackf77_ctrtri FORTRAN_NAME ( ctrtri , CTRTRI )
# Skipping MacroDefinition: lapackf77_cung2r FORTRAN_NAME ( cung2r , CUNG2R )
# Skipping MacroDefinition: lapackf77_cungbr FORTRAN_NAME ( cungbr , CUNGBR )
# Skipping MacroDefinition: lapackf77_cunghr FORTRAN_NAME ( cunghr , CUNGHR )
# Skipping MacroDefinition: lapackf77_cunglq FORTRAN_NAME ( cunglq , CUNGLQ )
# Skipping MacroDefinition: lapackf77_cungql FORTRAN_NAME ( cungql , CUNGQL )
# Skipping MacroDefinition: lapackf77_cungqr FORTRAN_NAME ( cungqr , CUNGQR )
# Skipping MacroDefinition: lapackf77_cungtr FORTRAN_NAME ( cungtr , CUNGTR )
# Skipping MacroDefinition: lapackf77_cunm2r FORTRAN_NAME ( cunm2r , CUNM2R )
# Skipping MacroDefinition: lapackf77_cunmbr FORTRAN_NAME ( cunmbr , CUNMBR )
# Skipping MacroDefinition: lapackf77_cunmlq FORTRAN_NAME ( cunmlq , CUNMLQ )
# Skipping MacroDefinition: lapackf77_cunmql FORTRAN_NAME ( cunmql , CUNMQL )
# Skipping MacroDefinition: lapackf77_cunmqr FORTRAN_NAME ( cunmqr , CUNMQR )
# Skipping MacroDefinition: lapackf77_cunmrq FORTRAN_NAME ( cunmrq , CUNMRQ )
# Skipping MacroDefinition: lapackf77_cunmtr FORTRAN_NAME ( cunmtr , CUNMTR )
# Skipping MacroDefinition: lapackf77_cbdt01 FORTRAN_NAME ( cbdt01 , CBDT01 )
# Skipping MacroDefinition: lapackf77_cget22 FORTRAN_NAME ( cget22 , CGET22 )
# Skipping MacroDefinition: lapackf77_chet21 FORTRAN_NAME ( chet21 , CHET21 )
# Skipping MacroDefinition: lapackf77_chst01 FORTRAN_NAME ( chst01 , CHST01 )
# Skipping MacroDefinition: lapackf77_clarfy FORTRAN_NAME ( clarfy , CLARFY )
# Skipping MacroDefinition: lapackf77_clatms FORTRAN_NAME ( clatms , CLATMS )
# Skipping MacroDefinition: lapackf77_cqpt01 FORTRAN_NAME ( cqpt01 , CQPT01 )
# Skipping MacroDefinition: lapackf77_cqrt02 FORTRAN_NAME ( cqrt02 , CQRT02 )
# Skipping MacroDefinition: lapackf77_cstt21 FORTRAN_NAME ( cstt21 , CSTT21 )
# Skipping MacroDefinition: lapackf77_cunt01 FORTRAN_NAME ( cunt01 , CUNT01 )
# Skipping MacroDefinition: magma_setvector ( n , elemSize , hx_src , incx , dy_dst , incy , queue ) magma_setvector_internal ( n , elemSize , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_getvector ( n , elemSize , dx_src , incx , hy_dst , incy , queue ) magma_getvector_internal ( n , elemSize , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_copyvector ( n , elemSize , dx_src , incx , dy_dst , incy , queue ) magma_copyvector_internal ( n , elemSize , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_setvector_async ( n , elemSize , hx_src , incx , dy_dst , incy , queue ) magma_setvector_async_internal ( n , elemSize , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_getvector_async ( n , elemSize , dx_src , incx , hy_dst , incy , queue ) magma_getvector_async_internal ( n , elemSize , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_copyvector_async ( n , elemSize , dx_src , incx , dy_dst , incy , queue ) magma_copyvector_async_internal ( n , elemSize , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_setmatrix ( m , n , elemSize , hA_src , lda , dB_dst , lddb , queue ) magma_setmatrix_internal ( m , n , elemSize , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_getmatrix ( m , n , elemSize , dA_src , ldda , hB_dst , ldb , queue ) magma_getmatrix_internal ( m , n , elemSize , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_copymatrix ( m , n , elemSize , dA_src , ldda , dB_dst , lddb , queue ) magma_copymatrix_internal ( m , n , elemSize , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_setmatrix_async ( m , n , elemSize , hA_src , lda , dB_dst , lddb , queue ) magma_setmatrix_async_internal ( m , n , elemSize , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_getmatrix_async ( m , n , elemSize , dA_src , ldda , hB_dst , ldb , queue ) magma_getmatrix_async_internal ( m , n , elemSize , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_copymatrix_async ( m , n , elemSize , dA_src , ldda , dB_dst , lddb , queue ) magma_copymatrix_async_internal ( m , n , elemSize , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_isetvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_isetvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_igetvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_igetvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_icopyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_icopyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_isetvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_isetvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_igetvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_igetvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_icopyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_icopyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_isetmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_isetmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_igetmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_igetmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_icopymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_icopymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_isetmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_isetmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_igetmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_igetmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_icopymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_icopymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_setvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_index_setvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_getvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_index_getvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_copyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_index_copyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_setvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_index_setvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_getvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_index_getvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_copyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_index_copyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_uindex_setvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_uindex_setvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_uindex_getvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_uindex_getvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_uindex_copyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_uindex_copyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_uindex_setvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_uindex_setvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_uindex_getvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_uindex_getvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_uindex_copyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_uindex_copyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_setmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_index_setmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_getmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_index_getmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_copymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_index_copymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_setmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_index_setmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_getmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_index_getmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_copymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_index_copymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_setvector_v1 ( n , elemSize , hx_src , incx , dy_dst , incy ) magma_setvector_v1_internal ( n , elemSize , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_getvector_v1 ( n , elemSize , dx_src , incx , hy_dst , incy ) magma_getvector_v1_internal ( n , elemSize , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_copyvector_v1 ( n , elemSize , dx_src , incx , dy_dst , incy ) magma_copyvector_v1_internal ( n , elemSize , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_setmatrix_v1 ( m , n , elemSize , hA_src , lda , dB_dst , lddb ) magma_setmatrix_v1_internal ( m , n , elemSize , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_getmatrix_v1 ( m , n , elemSize , dA_src , ldda , hB_dst , ldb ) magma_getmatrix_v1_internal ( m , n , elemSize , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_copymatrix_v1 ( m , n , elemSize , dA_src , ldda , dB_dst , lddb ) magma_copymatrix_v1_internal ( m , n , elemSize , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_isetvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_isetvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_igetvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_igetvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_icopyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_icopyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_isetmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_isetmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_igetmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_igetmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_icopymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_icopymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_setvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_index_setvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_getvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_index_getvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_copyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_index_copyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_setmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_index_setmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_getmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_index_getmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_index_copymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_index_copymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )

const magmaDouble_ptr = Ptr{Cdouble}

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

# Skipping MacroDefinition: blasf77_idamax FORTRAN_NAME ( idamax , IDAMAX )
# Skipping MacroDefinition: blasf77_daxpy FORTRAN_NAME ( daxpy , DAXPY )
# Skipping MacroDefinition: blasf77_dcopy FORTRAN_NAME ( dcopy , DCOPY )
# Skipping MacroDefinition: blasf77_dgemm FORTRAN_NAME ( dgemm , DGEMM )
# Skipping MacroDefinition: blasf77_dgemv FORTRAN_NAME ( dgemv , DGEMV )
# Skipping MacroDefinition: blasf77_dger FORTRAN_NAME ( dger , DGER )
# Skipping MacroDefinition: blasf77_dsymm FORTRAN_NAME ( dsymm , DSYMM )
# Skipping MacroDefinition: blasf77_dsymv FORTRAN_NAME ( dsymv , DSYMV )
# Skipping MacroDefinition: blasf77_dsyr FORTRAN_NAME ( dsyr , DSYR )
# Skipping MacroDefinition: blasf77_dsyr2 FORTRAN_NAME ( dsyr2 , DSYR2 )
# Skipping MacroDefinition: blasf77_dsyr2k FORTRAN_NAME ( dsyr2k , DSYR2K )
# Skipping MacroDefinition: blasf77_dsyrk FORTRAN_NAME ( dsyrk , DSYRK )
# Skipping MacroDefinition: blasf77_dscal FORTRAN_NAME ( dscal , DSCAL )
# Skipping MacroDefinition: blasf77_dswap FORTRAN_NAME ( dswap , DSWAP )
# Skipping MacroDefinition: blasf77_drotg FORTRAN_NAME ( drotg , DROTG )
# Skipping MacroDefinition: blasf77_drot FORTRAN_NAME ( drot , DROT )
# Skipping MacroDefinition: blasf77_dtrmm FORTRAN_NAME ( dtrmm , DTRMM )
# Skipping MacroDefinition: blasf77_dtrmv FORTRAN_NAME ( dtrmv , DTRMV )
# Skipping MacroDefinition: blasf77_dtrsm FORTRAN_NAME ( dtrsm , DTRSM )
# Skipping MacroDefinition: blasf77_dtrsv FORTRAN_NAME ( dtrsv , DTRSV )
# Skipping MacroDefinition: lapackf77_dlaed2 FORTRAN_NAME ( dlaed2 , DLAED2 )
# Skipping MacroDefinition: lapackf77_dlaed4 FORTRAN_NAME ( dlaed4 , DLAED4 )
# Skipping MacroDefinition: lapackf77_dlaln2 FORTRAN_NAME ( dlaln2 , DLALN2 )
# Skipping MacroDefinition: lapackf77_dlamc3 FORTRAN_NAME ( dlamc3 , DLAMC3 )
# Skipping MacroDefinition: lapackf77_dlamrg FORTRAN_NAME ( dlamrg , DLAMRG )
# Skipping MacroDefinition: lapackf77_dlasrt FORTRAN_NAME ( dlasrt , DLASRT )
# Skipping MacroDefinition: lapackf77_dstebz FORTRAN_NAME ( dstebz , DSTEBZ )
# Skipping MacroDefinition: lapackf77_dbdsdc FORTRAN_NAME ( dbdsdc , DBDSDC )
# Skipping MacroDefinition: lapackf77_dbdsqr FORTRAN_NAME ( dbdsqr , DBDSQR )
# Skipping MacroDefinition: lapackf77_dgebak FORTRAN_NAME ( dgebak , DGEBAK )
# Skipping MacroDefinition: lapackf77_dgebal FORTRAN_NAME ( dgebal , DGEBAL )
# Skipping MacroDefinition: lapackf77_dgebd2 FORTRAN_NAME ( dgebd2 , DGEBD2 )
# Skipping MacroDefinition: lapackf77_dgebrd FORTRAN_NAME ( dgebrd , DGEBRD )
# Skipping MacroDefinition: lapackf77_dgbbrd FORTRAN_NAME ( dgbbrd , DGBBRD )
# Skipping MacroDefinition: lapackf77_dgbsv FORTRAN_NAME ( dgbsv , DGBSV )
# Skipping MacroDefinition: lapackf77_dgeev FORTRAN_NAME ( dgeev , DGEEV )
# Skipping MacroDefinition: lapackf77_dgehd2 FORTRAN_NAME ( dgehd2 , DGEHD2 )
# Skipping MacroDefinition: lapackf77_dgehrd FORTRAN_NAME ( dgehrd , DGEHRD )
# Skipping MacroDefinition: lapackf77_dgelqf FORTRAN_NAME ( dgelqf , DGELQF )
# Skipping MacroDefinition: lapackf77_dgels FORTRAN_NAME ( dgels , DGELS )
# Skipping MacroDefinition: lapackf77_dgeqlf FORTRAN_NAME ( dgeqlf , DGEQLF )
# Skipping MacroDefinition: lapackf77_dgeqp3 FORTRAN_NAME ( dgeqp3 , DGEQP3 )
# Skipping MacroDefinition: lapackf77_dgeqrf FORTRAN_NAME ( dgeqrf , DGEQRF )
# Skipping MacroDefinition: lapackf77_dgerqf FORTRAN_NAME ( dgerqf , DGERQF )
# Skipping MacroDefinition: lapackf77_dgesdd FORTRAN_NAME ( dgesdd , DGESDD )
# Skipping MacroDefinition: lapackf77_dgesv FORTRAN_NAME ( dgesv , DGESV )
# Skipping MacroDefinition: lapackf77_dgesvd FORTRAN_NAME ( dgesvd , DGESVD )
# Skipping MacroDefinition: lapackf77_dgetrf FORTRAN_NAME ( dgetrf , DGETRF )
# Skipping MacroDefinition: lapackf77_dgetri FORTRAN_NAME ( dgetri , DGETRI )
# Skipping MacroDefinition: lapackf77_dgetrs FORTRAN_NAME ( dgetrs , DGETRS )
# Skipping MacroDefinition: lapackf77_dgglse FORTRAN_NAME ( dgglse , DGGLSE )
# Skipping MacroDefinition: lapackf77_dggrqf FORTRAN_NAME ( dggrqf , DGGRQF )
# Skipping MacroDefinition: lapackf77_dsytf2 FORTRAN_NAME ( dsytf2 , DSYTF2 )
# Skipping MacroDefinition: lapackf77_dsytrs FORTRAN_NAME ( dsytrs , DSYTRS )
# Skipping MacroDefinition: lapackf77_dsbtrd FORTRAN_NAME ( dsbtrd , DSBTRD )
# Skipping MacroDefinition: lapackf77_dsyev FORTRAN_NAME ( dsyev , DSYEV )
# Skipping MacroDefinition: lapackf77_dsyevd FORTRAN_NAME ( dsyevd , DSYEVD )
# Skipping MacroDefinition: lapackf77_dsyevr FORTRAN_NAME ( dsyevr , DSYEVR )
# Skipping MacroDefinition: lapackf77_dsyevx FORTRAN_NAME ( dsyevx , DSYEVX )
# Skipping MacroDefinition: lapackf77_dsygs2 FORTRAN_NAME ( dsygs2 , DSYGS2 )
# Skipping MacroDefinition: lapackf77_dsygst FORTRAN_NAME ( dsygst , DSYGST )
# Skipping MacroDefinition: lapackf77_dsygvd FORTRAN_NAME ( dsygvd , DSYGVD )
# Skipping MacroDefinition: lapackf77_dsytd2 FORTRAN_NAME ( dsytd2 , DSYTD2 )
# Skipping MacroDefinition: lapackf77_dsytrd FORTRAN_NAME ( dsytrd , DSYTRD )
# Skipping MacroDefinition: lapackf77_dsytrf FORTRAN_NAME ( dsytrf , DSYTRF )
# Skipping MacroDefinition: lapackf77_dsysv FORTRAN_NAME ( dsysv , DSYSV )
# Skipping MacroDefinition: lapackf77_dhseqr FORTRAN_NAME ( dhseqr , DHSEQR )
# Skipping MacroDefinition: lapackf77_dlabrd FORTRAN_NAME ( dlabrd , DLABRD )
# Skipping MacroDefinition: lapackf77_dlacgv FORTRAN_NAME ( dlacgv , DLACGV )
# Skipping MacroDefinition: lapackf77_dlacp2 FORTRAN_NAME ( dlacp2 , DLACP2 )
# Skipping MacroDefinition: lapackf77_dlacpy FORTRAN_NAME ( dlacpy , DLACPY )
# Skipping MacroDefinition: lapackf77_dlacrm FORTRAN_NAME ( dlacrm , DLACRM )
# Skipping MacroDefinition: lapackf77_dladiv FORTRAN_NAME ( dladiv , DLADIV )
# Skipping MacroDefinition: lapackf77_dlasyf FORTRAN_NAME ( dlasyf , DLASYF )
# Skipping MacroDefinition: lapackf77_dlange FORTRAN_NAME ( dlange , DLANGE )
# Skipping MacroDefinition: lapackf77_dlansy FORTRAN_NAME ( dlansy , DLANSY )
# Skipping MacroDefinition: lapackf77_dlanst FORTRAN_NAME ( dlanst , DLANST )
# Skipping MacroDefinition: lapackf77_dlantr FORTRAN_NAME ( dlantr , DLANTR )
# Skipping MacroDefinition: lapackf77_dlapy3 FORTRAN_NAME ( dlapy3 , DLAPY3 )
# Skipping MacroDefinition: lapackf77_dlaqp2 FORTRAN_NAME ( dlaqp2 , DLAQP2 )
# Skipping MacroDefinition: lapackf77_dlarcm FORTRAN_NAME ( dlarcm , DLARCM )
# Skipping MacroDefinition: lapackf77_dlarf FORTRAN_NAME ( dlarf , DLARF )
# Skipping MacroDefinition: lapackf77_dlarfb FORTRAN_NAME ( dlarfb , DLARFB )
# Skipping MacroDefinition: lapackf77_dlarfg FORTRAN_NAME ( dlarfg , DLARFG )
# Skipping MacroDefinition: lapackf77_dlarft FORTRAN_NAME ( dlarft , DLARFT )
# Skipping MacroDefinition: lapackf77_dlarfx FORTRAN_NAME ( dlarfx , DLARFX )
# Skipping MacroDefinition: lapackf77_dlarnv FORTRAN_NAME ( dlarnv , DLARNV )
# Skipping MacroDefinition: lapackf77_dlartg FORTRAN_NAME ( dlartg , DLARTG )
# Skipping MacroDefinition: lapackf77_dlascl FORTRAN_NAME ( dlascl , DLASCL )
# Skipping MacroDefinition: lapackf77_dlaset FORTRAN_NAME ( dlaset , DLASET )
# Skipping MacroDefinition: lapackf77_dlaswp FORTRAN_NAME ( dlaswp , DLASWP )
# Skipping MacroDefinition: lapackf77_dlatrd FORTRAN_NAME ( dlatrd , DLATRD )
# Skipping MacroDefinition: lapackf77_dlatrs FORTRAN_NAME ( dlatrs , DLATRS )
# Skipping MacroDefinition: lapackf77_dlauum FORTRAN_NAME ( dlauum , DLAUUM )
# Skipping MacroDefinition: lapackf77_dlavsy FORTRAN_NAME ( dlavsy , DLAVSY )
# Skipping MacroDefinition: lapackf77_dposv FORTRAN_NAME ( dposv , DPOSV )
# Skipping MacroDefinition: lapackf77_dpotrf FORTRAN_NAME ( dpotrf , DPOTRF )
# Skipping MacroDefinition: lapackf77_dpotri FORTRAN_NAME ( dpotri , DPOTRI )
# Skipping MacroDefinition: lapackf77_dpotrs FORTRAN_NAME ( dpotrs , DPOTRS )
# Skipping MacroDefinition: lapackf77_dstedc FORTRAN_NAME ( dstedc , DSTEDC )
# Skipping MacroDefinition: lapackf77_dstein FORTRAN_NAME ( dstein , DSTEIN )
# Skipping MacroDefinition: lapackf77_dstemr FORTRAN_NAME ( dstemr , DSTEMR )
# Skipping MacroDefinition: lapackf77_dsteqr FORTRAN_NAME ( dsteqr , DSTEQR )
# Skipping MacroDefinition: lapackf77_dsymv FORTRAN_NAME ( dsymv , DSYMV )
# Skipping MacroDefinition: lapackf77_dsyr FORTRAN_NAME ( dsyr , DSYR )
# Skipping MacroDefinition: lapackf77_dtrevc FORTRAN_NAME ( dtrevc , DTREVC )
# Skipping MacroDefinition: lapackf77_dtrevc3 FORTRAN_NAME ( dtrevc3 , DTREVC3 )
# Skipping MacroDefinition: lapackf77_dtrtri FORTRAN_NAME ( dtrtri , DTRTRI )
# Skipping MacroDefinition: lapackf77_dorg2r FORTRAN_NAME ( dorg2r , DORG2R )
# Skipping MacroDefinition: lapackf77_dorgbr FORTRAN_NAME ( dorgbr , DORGBR )
# Skipping MacroDefinition: lapackf77_dorghr FORTRAN_NAME ( dorghr , DORGHR )
# Skipping MacroDefinition: lapackf77_dorglq FORTRAN_NAME ( dorglq , DORGLQ )
# Skipping MacroDefinition: lapackf77_dorgql FORTRAN_NAME ( dorgql , DORGQL )
# Skipping MacroDefinition: lapackf77_dorgqr FORTRAN_NAME ( dorgqr , DORGQR )
# Skipping MacroDefinition: lapackf77_dorgtr FORTRAN_NAME ( dorgtr , DORGTR )
# Skipping MacroDefinition: lapackf77_dorm2r FORTRAN_NAME ( dorm2r , DORM2R )
# Skipping MacroDefinition: lapackf77_dormbr FORTRAN_NAME ( dormbr , DORMBR )
# Skipping MacroDefinition: lapackf77_dormlq FORTRAN_NAME ( dormlq , DORMLQ )
# Skipping MacroDefinition: lapackf77_dormql FORTRAN_NAME ( dormql , DORMQL )
# Skipping MacroDefinition: lapackf77_dormqr FORTRAN_NAME ( dormqr , DORMQR )
# Skipping MacroDefinition: lapackf77_dormrq FORTRAN_NAME ( dormrq , DORMRQ )
# Skipping MacroDefinition: lapackf77_dormtr FORTRAN_NAME ( dormtr , DORMTR )
# Skipping MacroDefinition: lapackf77_dbdt01 FORTRAN_NAME ( dbdt01 , DBDT01 )
# Skipping MacroDefinition: lapackf77_dget22 FORTRAN_NAME ( dget22 , DGET22 )
# Skipping MacroDefinition: lapackf77_dsyt21 FORTRAN_NAME ( dsyt21 , DSYT21 )
# Skipping MacroDefinition: lapackf77_dhst01 FORTRAN_NAME ( dhst01 , DHST01 )
# Skipping MacroDefinition: lapackf77_dlarfy FORTRAN_NAME ( dlarfy , DLARFY )
# Skipping MacroDefinition: lapackf77_dlatms FORTRAN_NAME ( dlatms , DLATMS )
# Skipping MacroDefinition: lapackf77_dqpt01 FORTRAN_NAME ( dqpt01 , DQPT01 )
# Skipping MacroDefinition: lapackf77_dqrt02 FORTRAN_NAME ( dqrt02 , DQRT02 )
# Skipping MacroDefinition: lapackf77_dstt21 FORTRAN_NAME ( dstt21 , DSTT21 )
# Skipping MacroDefinition: lapackf77_dort01 FORTRAN_NAME ( dort01 , DORT01 )
# Skipping MacroDefinition: lapackf77_ieeeck FORTRAN_NAME ( ieeeck , IEEECK )
# Skipping MacroDefinition: lapackf77_lsame FORTRAN_NAME ( lsame , LSAME )
# Skipping MacroDefinition: lapackf77_slamch FORTRAN_NAME ( slamch , SLAMCH )
# Skipping MacroDefinition: lapackf77_dlamch FORTRAN_NAME ( dlamch , DLAMCH )
# Skipping MacroDefinition: lapackf77_slabad FORTRAN_NAME ( slabad , SLABAD )
# Skipping MacroDefinition: lapackf77_dlabad FORTRAN_NAME ( dlabad , DLABAD )
# Skipping MacroDefinition: lapackf77_zcgesv FORTRAN_NAME ( zcgesv , ZCGESV )
# Skipping MacroDefinition: lapackf77_dsgesv FORTRAN_NAME ( dsgesv , DSGESV )
# Skipping MacroDefinition: lapackf77_dsterf FORTRAN_NAME ( dsterf , DSTERF )
# Skipping MacroDefinition: lapackf77_ssterf FORTRAN_NAME ( ssterf , SSTERF )
# Skipping MacroDefinition: lapackf77_zlag2c FORTRAN_NAME ( zlag2c , ZLAG2C )
# Skipping MacroDefinition: lapackf77_clag2z FORTRAN_NAME ( clag2z , CLAG2Z )
# Skipping MacroDefinition: lapackf77_dlag2s FORTRAN_NAME ( dlag2s , DLAG2S )
# Skipping MacroDefinition: lapackf77_slag2d FORTRAN_NAME ( slag2d , SLAG2D )
# Skipping MacroDefinition: lapackf77_zlat2c FORTRAN_NAME ( zlat2c , ZLAT2C )
# Skipping MacroDefinition: lapackf77_dlat2s FORTRAN_NAME ( dlat2s , DLAT2S )
# Skipping MacroDefinition: lapackf77_dlapy2 FORTRAN_NAME ( dlapy2 , DLAPY2 )
# Skipping MacroDefinition: lapackf77_slapy2 FORTRAN_NAME ( slapy2 , SLAPY2 )
# Skipping MacroDefinition: FORTRAN_NAME ( lcname , UCNAME ) lcname ## _error

const magmaFloat_ptr = Ptr{Cfloat}

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

# Skipping MacroDefinition: blasf77_isamax FORTRAN_NAME ( isamax , ISAMAX )
# Skipping MacroDefinition: blasf77_saxpy FORTRAN_NAME ( saxpy , SAXPY )
# Skipping MacroDefinition: blasf77_scopy FORTRAN_NAME ( scopy , SCOPY )
# Skipping MacroDefinition: blasf77_sgemm FORTRAN_NAME ( sgemm , SGEMM )
# Skipping MacroDefinition: blasf77_sgemv FORTRAN_NAME ( sgemv , SGEMV )
# Skipping MacroDefinition: blasf77_sger FORTRAN_NAME ( sger , SGER )
# Skipping MacroDefinition: blasf77_ssymm FORTRAN_NAME ( ssymm , SSYMM )
# Skipping MacroDefinition: blasf77_ssymv FORTRAN_NAME ( ssymv , SSYMV )
# Skipping MacroDefinition: blasf77_ssyr FORTRAN_NAME ( ssyr , SSYR )
# Skipping MacroDefinition: blasf77_ssyr2 FORTRAN_NAME ( ssyr2 , SSYR2 )
# Skipping MacroDefinition: blasf77_ssyr2k FORTRAN_NAME ( ssyr2k , SSYR2K )
# Skipping MacroDefinition: blasf77_ssyrk FORTRAN_NAME ( ssyrk , SSYRK )
# Skipping MacroDefinition: blasf77_sscal FORTRAN_NAME ( sscal , SSCAL )
# Skipping MacroDefinition: blasf77_sswap FORTRAN_NAME ( sswap , SSWAP )
# Skipping MacroDefinition: blasf77_srotg FORTRAN_NAME ( srotg , SROTG )
# Skipping MacroDefinition: blasf77_srot FORTRAN_NAME ( srot , SROT )
# Skipping MacroDefinition: blasf77_strmm FORTRAN_NAME ( strmm , STRMM )
# Skipping MacroDefinition: blasf77_strmv FORTRAN_NAME ( strmv , STRMV )
# Skipping MacroDefinition: blasf77_strsm FORTRAN_NAME ( strsm , STRSM )
# Skipping MacroDefinition: blasf77_strsv FORTRAN_NAME ( strsv , STRSV )
# Skipping MacroDefinition: lapackf77_sbdsqr FORTRAN_NAME ( sbdsqr , SBDSQR )
# Skipping MacroDefinition: lapackf77_sgebak FORTRAN_NAME ( sgebak , SGEBAK )
# Skipping MacroDefinition: lapackf77_sgebal FORTRAN_NAME ( sgebal , SGEBAL )
# Skipping MacroDefinition: lapackf77_sgebd2 FORTRAN_NAME ( sgebd2 , SGEBD2 )
# Skipping MacroDefinition: lapackf77_sgebrd FORTRAN_NAME ( sgebrd , SGEBRD )
# Skipping MacroDefinition: lapackf77_sgbbrd FORTRAN_NAME ( sgbbrd , SGBBRD )
# Skipping MacroDefinition: lapackf77_sgbsv FORTRAN_NAME ( sgbsv , SGBSV )
# Skipping MacroDefinition: lapackf77_sgeev FORTRAN_NAME ( sgeev , SGEEV )
# Skipping MacroDefinition: lapackf77_sgehd2 FORTRAN_NAME ( sgehd2 , SGEHD2 )
# Skipping MacroDefinition: lapackf77_sgehrd FORTRAN_NAME ( sgehrd , SGEHRD )
# Skipping MacroDefinition: lapackf77_sgelqf FORTRAN_NAME ( sgelqf , SGELQF )
# Skipping MacroDefinition: lapackf77_sgels FORTRAN_NAME ( sgels , SGELS )
# Skipping MacroDefinition: lapackf77_sgeqlf FORTRAN_NAME ( sgeqlf , SGEQLF )
# Skipping MacroDefinition: lapackf77_sgeqp3 FORTRAN_NAME ( sgeqp3 , SGEQP3 )
# Skipping MacroDefinition: lapackf77_sgeqrf FORTRAN_NAME ( sgeqrf , SGEQRF )
# Skipping MacroDefinition: lapackf77_sgerqf FORTRAN_NAME ( sgerqf , SGERQF )
# Skipping MacroDefinition: lapackf77_sgesdd FORTRAN_NAME ( sgesdd , SGESDD )
# Skipping MacroDefinition: lapackf77_sgesv FORTRAN_NAME ( sgesv , SGESV )
# Skipping MacroDefinition: lapackf77_sgesvd FORTRAN_NAME ( sgesvd , SGESVD )
# Skipping MacroDefinition: lapackf77_sgetrf FORTRAN_NAME ( sgetrf , SGETRF )
# Skipping MacroDefinition: lapackf77_sgetri FORTRAN_NAME ( sgetri , SGETRI )
# Skipping MacroDefinition: lapackf77_sgetrs FORTRAN_NAME ( sgetrs , SGETRS )
# Skipping MacroDefinition: lapackf77_sgglse FORTRAN_NAME ( sgglse , SGGLSE )
# Skipping MacroDefinition: lapackf77_sggrqf FORTRAN_NAME ( sggrqf , SGGRQF )
# Skipping MacroDefinition: lapackf77_ssytf2 FORTRAN_NAME ( ssytf2 , SSYTF2 )
# Skipping MacroDefinition: lapackf77_ssytrs FORTRAN_NAME ( ssytrs , SSYTRS )
# Skipping MacroDefinition: lapackf77_ssbtrd FORTRAN_NAME ( ssbtrd , SSBTRD )
# Skipping MacroDefinition: lapackf77_ssyev FORTRAN_NAME ( ssyev , SSYEV )
# Skipping MacroDefinition: lapackf77_ssyevd FORTRAN_NAME ( ssyevd , SSYEVD )
# Skipping MacroDefinition: lapackf77_ssyevr FORTRAN_NAME ( ssyevr , SSYEVR )
# Skipping MacroDefinition: lapackf77_ssyevx FORTRAN_NAME ( ssyevx , SSYEVX )
# Skipping MacroDefinition: lapackf77_ssygs2 FORTRAN_NAME ( ssygs2 , SSYGS2 )
# Skipping MacroDefinition: lapackf77_ssygst FORTRAN_NAME ( ssygst , SSYGST )
# Skipping MacroDefinition: lapackf77_ssygvd FORTRAN_NAME ( ssygvd , SSYGVD )
# Skipping MacroDefinition: lapackf77_ssytd2 FORTRAN_NAME ( ssytd2 , SSYTD2 )
# Skipping MacroDefinition: lapackf77_ssytrd FORTRAN_NAME ( ssytrd , SSYTRD )
# Skipping MacroDefinition: lapackf77_ssytrf FORTRAN_NAME ( ssytrf , SSYTRF )
# Skipping MacroDefinition: lapackf77_ssysv FORTRAN_NAME ( ssysv , SSYSV )
# Skipping MacroDefinition: lapackf77_shseqr FORTRAN_NAME ( shseqr , SHSEQR )
# Skipping MacroDefinition: lapackf77_slabrd FORTRAN_NAME ( slabrd , SLABRD )
# Skipping MacroDefinition: lapackf77_slacgv FORTRAN_NAME ( slacgv , SLACGV )
# Skipping MacroDefinition: lapackf77_slacp2 FORTRAN_NAME ( slacp2 , SLACP2 )
# Skipping MacroDefinition: lapackf77_slacpy FORTRAN_NAME ( slacpy , SLACPY )
# Skipping MacroDefinition: lapackf77_slacrm FORTRAN_NAME ( slacrm , SLACRM )
# Skipping MacroDefinition: lapackf77_sladiv FORTRAN_NAME ( sladiv , SLADIV )
# Skipping MacroDefinition: lapackf77_slasyf FORTRAN_NAME ( slasyf , SLASYF )
# Skipping MacroDefinition: lapackf77_slange FORTRAN_NAME ( slange , SLANGE )
# Skipping MacroDefinition: lapackf77_slansy FORTRAN_NAME ( slansy , SLANSY )
# Skipping MacroDefinition: lapackf77_slanst FORTRAN_NAME ( slanst , SLANST )
# Skipping MacroDefinition: lapackf77_slantr FORTRAN_NAME ( slantr , SLANTR )
# Skipping MacroDefinition: lapackf77_slaqp2 FORTRAN_NAME ( slaqp2 , SLAQP2 )
# Skipping MacroDefinition: lapackf77_slarcm FORTRAN_NAME ( slarcm , SLARCM )
# Skipping MacroDefinition: lapackf77_slarf FORTRAN_NAME ( slarf , SLARF )
# Skipping MacroDefinition: lapackf77_slarfb FORTRAN_NAME ( slarfb , SLARFB )
# Skipping MacroDefinition: lapackf77_slarfg FORTRAN_NAME ( slarfg , SLARFG )
# Skipping MacroDefinition: lapackf77_slarft FORTRAN_NAME ( slarft , SLARFT )
# Skipping MacroDefinition: lapackf77_slarfx FORTRAN_NAME ( slarfx , SLARFX )
# Skipping MacroDefinition: lapackf77_slarnv FORTRAN_NAME ( slarnv , SLARNV )
# Skipping MacroDefinition: lapackf77_slartg FORTRAN_NAME ( slartg , SLARTG )
# Skipping MacroDefinition: lapackf77_slascl FORTRAN_NAME ( slascl , SLASCL )
# Skipping MacroDefinition: lapackf77_slaset FORTRAN_NAME ( slaset , SLASET )
# Skipping MacroDefinition: lapackf77_slaswp FORTRAN_NAME ( slaswp , SLASWP )
# Skipping MacroDefinition: lapackf77_slatrd FORTRAN_NAME ( slatrd , SLATRD )
# Skipping MacroDefinition: lapackf77_slatrs FORTRAN_NAME ( slatrs , SLATRS )
# Skipping MacroDefinition: lapackf77_slauum FORTRAN_NAME ( slauum , SLAUUM )
# Skipping MacroDefinition: lapackf77_slavsy FORTRAN_NAME ( slavsy , SLAVSY )
# Skipping MacroDefinition: lapackf77_sposv FORTRAN_NAME ( sposv , SPOSV )
# Skipping MacroDefinition: lapackf77_spotrf FORTRAN_NAME ( spotrf , SPOTRF )
# Skipping MacroDefinition: lapackf77_spotri FORTRAN_NAME ( spotri , SPOTRI )
# Skipping MacroDefinition: lapackf77_spotrs FORTRAN_NAME ( spotrs , SPOTRS )
# Skipping MacroDefinition: lapackf77_sstedc FORTRAN_NAME ( sstedc , SSTEDC )
# Skipping MacroDefinition: lapackf77_sstein FORTRAN_NAME ( sstein , SSTEIN )
# Skipping MacroDefinition: lapackf77_sstemr FORTRAN_NAME ( sstemr , SSTEMR )
# Skipping MacroDefinition: lapackf77_ssteqr FORTRAN_NAME ( ssteqr , SSTEQR )
# Skipping MacroDefinition: lapackf77_ssymv FORTRAN_NAME ( ssymv , SSYMV )
# Skipping MacroDefinition: lapackf77_ssyr FORTRAN_NAME ( ssyr , SSYR )
# Skipping MacroDefinition: lapackf77_strevc FORTRAN_NAME ( strevc , STREVC )
# Skipping MacroDefinition: lapackf77_strevc3 FORTRAN_NAME ( strevc3 , STREVC3 )
# Skipping MacroDefinition: lapackf77_strtri FORTRAN_NAME ( strtri , STRTRI )
# Skipping MacroDefinition: lapackf77_sorg2r FORTRAN_NAME ( sorg2r , SORG2R )
# Skipping MacroDefinition: lapackf77_sorgbr FORTRAN_NAME ( sorgbr , SORGBR )
# Skipping MacroDefinition: lapackf77_sorghr FORTRAN_NAME ( sorghr , SORGHR )
# Skipping MacroDefinition: lapackf77_sorglq FORTRAN_NAME ( sorglq , SORGLQ )
# Skipping MacroDefinition: lapackf77_sorgql FORTRAN_NAME ( sorgql , SORGQL )
# Skipping MacroDefinition: lapackf77_sorgqr FORTRAN_NAME ( sorgqr , SORGQR )
# Skipping MacroDefinition: lapackf77_sorgtr FORTRAN_NAME ( sorgtr , SORGTR )
# Skipping MacroDefinition: lapackf77_sorm2r FORTRAN_NAME ( sorm2r , SORM2R )
# Skipping MacroDefinition: lapackf77_sormbr FORTRAN_NAME ( sormbr , SORMBR )
# Skipping MacroDefinition: lapackf77_sormlq FORTRAN_NAME ( sormlq , SORMLQ )
# Skipping MacroDefinition: lapackf77_sormql FORTRAN_NAME ( sormql , SORMQL )
# Skipping MacroDefinition: lapackf77_sormqr FORTRAN_NAME ( sormqr , SORMQR )
# Skipping MacroDefinition: lapackf77_sormrq FORTRAN_NAME ( sormrq , SORMRQ )
# Skipping MacroDefinition: lapackf77_sormtr FORTRAN_NAME ( sormtr , SORMTR )
# Skipping MacroDefinition: lapackf77_sbdt01 FORTRAN_NAME ( sbdt01 , SBDT01 )
# Skipping MacroDefinition: lapackf77_sget22 FORTRAN_NAME ( sget22 , SGET22 )
# Skipping MacroDefinition: lapackf77_ssyt21 FORTRAN_NAME ( ssyt21 , SSYT21 )
# Skipping MacroDefinition: lapackf77_shst01 FORTRAN_NAME ( shst01 , SHST01 )
# Skipping MacroDefinition: lapackf77_slarfy FORTRAN_NAME ( slarfy , SLARFY )
# Skipping MacroDefinition: lapackf77_slatms FORTRAN_NAME ( slatms , SLATMS )
# Skipping MacroDefinition: lapackf77_sqpt01 FORTRAN_NAME ( sqpt01 , SQPT01 )
# Skipping MacroDefinition: lapackf77_sqrt02 FORTRAN_NAME ( sqrt02 , SQRT02 )
# Skipping MacroDefinition: lapackf77_sstt21 FORTRAN_NAME ( sstt21 , SSTT21 )
# Skipping MacroDefinition: lapackf77_sort01 FORTRAN_NAME ( sort01 , SORT01 )
# Skipping MacroDefinition: MAGMA_Z_MAKE ( r , i ) make_cuDoubleComplex ( r , i )
# Skipping MacroDefinition: MAGMA_Z_REAL ( a ) ( a ) . x
# Skipping MacroDefinition: MAGMA_Z_IMAG ( a ) ( a ) . y
# Skipping MacroDefinition: MAGMA_Z_ADD ( a , b ) cuCadd ( a , b )
# Skipping MacroDefinition: MAGMA_Z_SUB ( a , b ) cuCsub ( a , b )
# Skipping MacroDefinition: MAGMA_Z_MUL ( a , b ) cuCmul ( a , b )
# Skipping MacroDefinition: MAGMA_Z_DIV ( a , b ) cuCdiv ( a , b )
# Skipping MacroDefinition: MAGMA_Z_ABS ( a ) cuCabs ( a )
# Skipping MacroDefinition: MAGMA_Z_ABS1 ( a ) ( fabs ( ( a ) . x ) + fabs ( ( a ) . y ) )
# Skipping MacroDefinition: MAGMA_Z_CONJ ( a ) cuConj ( a )
# Skipping MacroDefinition: MAGMA_C_MAKE ( r , i ) make_cuFloatComplex ( r , i )
# Skipping MacroDefinition: MAGMA_C_REAL ( a ) ( a ) . x
# Skipping MacroDefinition: MAGMA_C_IMAG ( a ) ( a ) . y
# Skipping MacroDefinition: MAGMA_C_ADD ( a , b ) cuCaddf ( a , b )
# Skipping MacroDefinition: MAGMA_C_SUB ( a , b ) cuCsubf ( a , b )
# Skipping MacroDefinition: MAGMA_C_MUL ( a , b ) cuCmulf ( a , b )
# Skipping MacroDefinition: MAGMA_C_DIV ( a , b ) cuCdivf ( a , b )
# Skipping MacroDefinition: MAGMA_C_ABS ( a ) cuCabsf ( a )
# Skipping MacroDefinition: MAGMA_C_ABS1 ( a ) ( fabsf ( ( a ) . x ) + fabsf ( ( a ) . y ) )
# Skipping MacroDefinition: MAGMA_C_CONJ ( a ) cuConjf ( a )
# Skipping MacroDefinition: MAGMA_Z_EQUAL ( a , b ) ( MAGMA_Z_REAL ( a ) == MAGMA_Z_REAL ( b ) && MAGMA_Z_IMAG ( a ) == MAGMA_Z_IMAG ( b ) )
# Skipping MacroDefinition: MAGMA_Z_NEGATE ( a ) MAGMA_Z_MAKE ( - MAGMA_Z_REAL ( a ) , - MAGMA_Z_IMAG ( a ) )
# Skipping MacroDefinition: MAGMA_C_EQUAL ( a , b ) ( MAGMA_C_REAL ( a ) == MAGMA_C_REAL ( b ) && MAGMA_C_IMAG ( a ) == MAGMA_C_IMAG ( b ) )
# Skipping MacroDefinition: MAGMA_C_NEGATE ( a ) MAGMA_C_MAKE ( - MAGMA_C_REAL ( a ) , - MAGMA_C_IMAG ( a ) )
# Skipping MacroDefinition: MAGMA_D_MAKE ( r , i ) ( r )
# Skipping MacroDefinition: MAGMA_D_REAL ( x ) ( x )
# Skipping MacroDefinition: MAGMA_D_IMAG ( x ) ( 0.0 )
# Skipping MacroDefinition: MAGMA_D_ADD ( a , b ) ( ( a ) + ( b ) )
# Skipping MacroDefinition: MAGMA_D_SUB ( a , b ) ( ( a ) - ( b ) )
# Skipping MacroDefinition: MAGMA_D_MUL ( a , b ) ( ( a ) * ( b ) )
# Skipping MacroDefinition: MAGMA_D_DIV ( a , b ) ( ( a ) / ( b ) )
# Skipping MacroDefinition: MAGMA_D_ABS ( a ) ( ( a ) > 0 ? ( a ) : - ( a ) )
# Skipping MacroDefinition: MAGMA_D_ABS1 ( a ) ( ( a ) > 0 ? ( a ) : - ( a ) )
# Skipping MacroDefinition: MAGMA_D_CONJ ( a ) ( a )
# Skipping MacroDefinition: MAGMA_D_EQUAL ( a , b ) ( ( a ) == ( b ) )
# Skipping MacroDefinition: MAGMA_D_NEGATE ( a ) ( - a )
# Skipping MacroDefinition: MAGMA_S_MAKE ( r , i ) ( r )
# Skipping MacroDefinition: MAGMA_S_REAL ( x ) ( x )
# Skipping MacroDefinition: MAGMA_S_IMAG ( x ) ( 0.0 )
# Skipping MacroDefinition: MAGMA_S_ADD ( a , b ) ( ( a ) + ( b ) )
# Skipping MacroDefinition: MAGMA_S_SUB ( a , b ) ( ( a ) - ( b ) )
# Skipping MacroDefinition: MAGMA_S_MUL ( a , b ) ( ( a ) * ( b ) )
# Skipping MacroDefinition: MAGMA_S_DIV ( a , b ) ( ( a ) / ( b ) )
# Skipping MacroDefinition: MAGMA_S_ABS ( a ) ( ( a ) > 0 ? ( a ) : - ( a ) )
# Skipping MacroDefinition: MAGMA_S_ABS1 ( a ) ( ( a ) > 0 ? ( a ) : - ( a ) )
# Skipping MacroDefinition: MAGMA_S_CONJ ( a ) ( a )
# Skipping MacroDefinition: MAGMA_S_EQUAL ( a , b ) ( ( a ) == ( b ) )
# Skipping MacroDefinition: MAGMA_S_NEGATE ( a ) ( - a )
# Skipping MacroDefinition: MAGMA_Z_ZERO MAGMA_Z_MAKE ( 0.0 , 0.0 )
# Skipping MacroDefinition: MAGMA_Z_ONE MAGMA_Z_MAKE ( 1.0 , 0.0 )
# Skipping MacroDefinition: MAGMA_Z_HALF MAGMA_Z_MAKE ( 0.5 , 0.0 )
# Skipping MacroDefinition: MAGMA_Z_NEG_ONE MAGMA_Z_MAKE ( - 1.0 , 0.0 )
# Skipping MacroDefinition: MAGMA_Z_NEG_HALF MAGMA_Z_MAKE ( - 0.5 , 0.0 )
# Skipping MacroDefinition: MAGMA_C_ZERO MAGMA_C_MAKE ( 0.0 , 0.0 )
# Skipping MacroDefinition: MAGMA_C_ONE MAGMA_C_MAKE ( 1.0 , 0.0 )
# Skipping MacroDefinition: MAGMA_C_HALF MAGMA_C_MAKE ( 0.5 , 0.0 )
# Skipping MacroDefinition: MAGMA_C_NEG_ONE MAGMA_C_MAKE ( - 1.0 , 0.0 )
# Skipping MacroDefinition: MAGMA_C_NEG_HALF MAGMA_C_MAKE ( - 0.5 , 0.0 )

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

# Skipping MacroDefinition: CBLAS_SADDR ( a ) & ( a )

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
const magma_ptr = Ptr{Cvoid}
const magmaInt_ptr = Ptr{magma_int_t}
const magmaIndex_ptr = Ptr{magma_index_t}
const magmaUIndex_ptr = Ptr{magma_uindex_t}
const magmaDoubleComplex_ptr = Ptr{magmaDoubleComplex}
const magmaHalf_ptr = Ptr{magmaHalf}
const magma_const_ptr = Ptr{Cvoid}
const magmaInt_const_ptr = Ptr{magma_int_t}
const magmaIndex_const_ptr = Ptr{magma_index_t}
const magmaUIndex_const_ptr = Ptr{magma_uindex_t}
const magmaFloat_const_ptr = Ptr{Cfloat}
const magmaDouble_const_ptr = Ptr{Cdouble}
const magmaFloatComplex_const_ptr = Ptr{magmaFloatComplex}
const magmaDoubleComplex_const_ptr = Ptr{magmaDoubleComplex}
const magmaHalf_const_ptr = Ptr{magmaHalf}

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

# Skipping MacroDefinition: blasf77_izamax FORTRAN_NAME ( izamax , IZAMAX )
# Skipping MacroDefinition: blasf77_zaxpy FORTRAN_NAME ( zaxpy , ZAXPY )
# Skipping MacroDefinition: blasf77_zcopy FORTRAN_NAME ( zcopy , ZCOPY )
# Skipping MacroDefinition: blasf77_zgemm FORTRAN_NAME ( zgemm , ZGEMM )
# Skipping MacroDefinition: blasf77_zgemv FORTRAN_NAME ( zgemv , ZGEMV )
# Skipping MacroDefinition: blasf77_zgerc FORTRAN_NAME ( zgerc , ZGERC )
# Skipping MacroDefinition: blasf77_zgeru FORTRAN_NAME ( zgeru , ZGERU )
# Skipping MacroDefinition: blasf77_zhemm FORTRAN_NAME ( zhemm , ZHEMM )
# Skipping MacroDefinition: blasf77_zhemv FORTRAN_NAME ( zhemv , ZHEMV )
# Skipping MacroDefinition: blasf77_zher FORTRAN_NAME ( zher , ZHER )
# Skipping MacroDefinition: blasf77_zher2 FORTRAN_NAME ( zher2 , ZHER2 )
# Skipping MacroDefinition: blasf77_zher2k FORTRAN_NAME ( zher2k , ZHER2K )
# Skipping MacroDefinition: blasf77_zherk FORTRAN_NAME ( zherk , ZHERK )
# Skipping MacroDefinition: blasf77_zscal FORTRAN_NAME ( zscal , ZSCAL )
# Skipping MacroDefinition: blasf77_zdscal FORTRAN_NAME ( zdscal , ZDSCAL )
# Skipping MacroDefinition: blasf77_zswap FORTRAN_NAME ( zswap , ZSWAP )
# Skipping MacroDefinition: blasf77_zsymm FORTRAN_NAME ( zsymm , ZSYMM )
# Skipping MacroDefinition: blasf77_zsyr2k FORTRAN_NAME ( zsyr2k , ZSYR2K )
# Skipping MacroDefinition: blasf77_zsyrk FORTRAN_NAME ( zsyrk , ZSYRK )
# Skipping MacroDefinition: blasf77_zrotg FORTRAN_NAME ( zrotg , ZROTG )
# Skipping MacroDefinition: blasf77_zrot FORTRAN_NAME ( zrot , ZROT )
# Skipping MacroDefinition: blasf77_zdrot FORTRAN_NAME ( zdrot , ZDROT )
# Skipping MacroDefinition: blasf77_ztrmm FORTRAN_NAME ( ztrmm , ZTRMM )
# Skipping MacroDefinition: blasf77_ztrmv FORTRAN_NAME ( ztrmv , ZTRMV )
# Skipping MacroDefinition: blasf77_ztrsm FORTRAN_NAME ( ztrsm , ZTRSM )
# Skipping MacroDefinition: blasf77_ztrsv FORTRAN_NAME ( ztrsv , ZTRSV )
# Skipping MacroDefinition: lapackf77_zbdsqr FORTRAN_NAME ( zbdsqr , ZBDSQR )
# Skipping MacroDefinition: lapackf77_zgebak FORTRAN_NAME ( zgebak , ZGEBAK )
# Skipping MacroDefinition: lapackf77_zgebal FORTRAN_NAME ( zgebal , ZGEBAL )
# Skipping MacroDefinition: lapackf77_zgebd2 FORTRAN_NAME ( zgebd2 , ZGEBD2 )
# Skipping MacroDefinition: lapackf77_zgebrd FORTRAN_NAME ( zgebrd , ZGEBRD )
# Skipping MacroDefinition: lapackf77_zgbbrd FORTRAN_NAME ( zgbbrd , ZGBBRD )
# Skipping MacroDefinition: lapackf77_zgbsv FORTRAN_NAME ( zgbsv , ZGBSV )
# Skipping MacroDefinition: lapackf77_zgeev FORTRAN_NAME ( zgeev , ZGEEV )
# Skipping MacroDefinition: lapackf77_zgehd2 FORTRAN_NAME ( zgehd2 , ZGEHD2 )
# Skipping MacroDefinition: lapackf77_zgehrd FORTRAN_NAME ( zgehrd , ZGEHRD )
# Skipping MacroDefinition: lapackf77_zgelqf FORTRAN_NAME ( zgelqf , ZGELQF )
# Skipping MacroDefinition: lapackf77_zgels FORTRAN_NAME ( zgels , ZGELS )
# Skipping MacroDefinition: lapackf77_zgeqlf FORTRAN_NAME ( zgeqlf , ZGEQLF )
# Skipping MacroDefinition: lapackf77_zgeqp3 FORTRAN_NAME ( zgeqp3 , ZGEQP3 )
# Skipping MacroDefinition: lapackf77_zgeqrf FORTRAN_NAME ( zgeqrf , ZGEQRF )
# Skipping MacroDefinition: lapackf77_zgerqf FORTRAN_NAME ( zgerqf , ZGERQF )
# Skipping MacroDefinition: lapackf77_zgesdd FORTRAN_NAME ( zgesdd , ZGESDD )
# Skipping MacroDefinition: lapackf77_zgesv FORTRAN_NAME ( zgesv , ZGESV )
# Skipping MacroDefinition: lapackf77_zgesvd FORTRAN_NAME ( zgesvd , ZGESVD )
# Skipping MacroDefinition: lapackf77_zgetrf FORTRAN_NAME ( zgetrf , ZGETRF )
# Skipping MacroDefinition: lapackf77_zgetri FORTRAN_NAME ( zgetri , ZGETRI )
# Skipping MacroDefinition: lapackf77_zgetrs FORTRAN_NAME ( zgetrs , ZGETRS )
# Skipping MacroDefinition: lapackf77_zgglse FORTRAN_NAME ( zgglse , ZGGLSE )
# Skipping MacroDefinition: lapackf77_zggrqf FORTRAN_NAME ( zggrqf , ZGGRQF )
# Skipping MacroDefinition: lapackf77_zhetf2 FORTRAN_NAME ( zhetf2 , ZHETF2 )
# Skipping MacroDefinition: lapackf77_zhetrs FORTRAN_NAME ( zhetrs , ZHETRS )
# Skipping MacroDefinition: lapackf77_zhbtrd FORTRAN_NAME ( zhbtrd , ZHBTRD )
# Skipping MacroDefinition: lapackf77_zheev FORTRAN_NAME ( zheev , ZHEEV )
# Skipping MacroDefinition: lapackf77_zheevd FORTRAN_NAME ( zheevd , ZHEEVD )
# Skipping MacroDefinition: lapackf77_zheevr FORTRAN_NAME ( zheevr , ZHEEVR )
# Skipping MacroDefinition: lapackf77_zheevx FORTRAN_NAME ( zheevx , ZHEEVX )
# Skipping MacroDefinition: lapackf77_zhegs2 FORTRAN_NAME ( zhegs2 , ZHEGS2 )
# Skipping MacroDefinition: lapackf77_zhegst FORTRAN_NAME ( zhegst , ZHEGST )
# Skipping MacroDefinition: lapackf77_zhegvd FORTRAN_NAME ( zhegvd , ZHEGVD )
# Skipping MacroDefinition: lapackf77_zhetd2 FORTRAN_NAME ( zhetd2 , ZHETD2 )
# Skipping MacroDefinition: lapackf77_zhetrd FORTRAN_NAME ( zhetrd , ZHETRD )
# Skipping MacroDefinition: lapackf77_zhetrf FORTRAN_NAME ( zhetrf , ZHETRF )
# Skipping MacroDefinition: lapackf77_zhesv FORTRAN_NAME ( zhesv , ZHESV )
# Skipping MacroDefinition: lapackf77_zhseqr FORTRAN_NAME ( zhseqr , ZHSEQR )
# Skipping MacroDefinition: lapackf77_zlabrd FORTRAN_NAME ( zlabrd , ZLABRD )
# Skipping MacroDefinition: lapackf77_zlacgv FORTRAN_NAME ( zlacgv , ZLACGV )
# Skipping MacroDefinition: lapackf77_zlacp2 FORTRAN_NAME ( zlacp2 , ZLACP2 )
# Skipping MacroDefinition: lapackf77_zlacpy FORTRAN_NAME ( zlacpy , ZLACPY )
# Skipping MacroDefinition: lapackf77_zlacrm FORTRAN_NAME ( zlacrm , ZLACRM )
# Skipping MacroDefinition: lapackf77_zladiv FORTRAN_NAME ( zladiv , ZLADIV )
# Skipping MacroDefinition: lapackf77_zlahef FORTRAN_NAME ( zlahef , ZLAHEF )
# Skipping MacroDefinition: lapackf77_zlange FORTRAN_NAME ( zlange , ZLANGE )
# Skipping MacroDefinition: lapackf77_zlanhe FORTRAN_NAME ( zlanhe , ZLANHE )
# Skipping MacroDefinition: lapackf77_zlanht FORTRAN_NAME ( zlanht , ZLANHT )
# Skipping MacroDefinition: lapackf77_zlansy FORTRAN_NAME ( zlansy , ZLANSY )
# Skipping MacroDefinition: lapackf77_zlantr FORTRAN_NAME ( zlantr , ZLANTR )
# Skipping MacroDefinition: lapackf77_zlaqp2 FORTRAN_NAME ( zlaqp2 , ZLAQP2 )
# Skipping MacroDefinition: lapackf77_zlarcm FORTRAN_NAME ( zlarcm , ZLARCM )
# Skipping MacroDefinition: lapackf77_zlarf FORTRAN_NAME ( zlarf , ZLARF )
# Skipping MacroDefinition: lapackf77_zlarfb FORTRAN_NAME ( zlarfb , ZLARFB )
# Skipping MacroDefinition: lapackf77_zlarfg FORTRAN_NAME ( zlarfg , ZLARFG )
# Skipping MacroDefinition: lapackf77_zlarft FORTRAN_NAME ( zlarft , ZLARFT )
# Skipping MacroDefinition: lapackf77_zlarfx FORTRAN_NAME ( zlarfx , ZLARFX )
# Skipping MacroDefinition: lapackf77_zlarnv FORTRAN_NAME ( zlarnv , ZLARNV )
# Skipping MacroDefinition: lapackf77_zlartg FORTRAN_NAME ( zlartg , ZLARTG )
# Skipping MacroDefinition: lapackf77_zlascl FORTRAN_NAME ( zlascl , ZLASCL )
# Skipping MacroDefinition: lapackf77_zlaset FORTRAN_NAME ( zlaset , ZLASET )
# Skipping MacroDefinition: lapackf77_zlaswp FORTRAN_NAME ( zlaswp , ZLASWP )
# Skipping MacroDefinition: lapackf77_zlatrd FORTRAN_NAME ( zlatrd , ZLATRD )
# Skipping MacroDefinition: lapackf77_zlatrs FORTRAN_NAME ( zlatrs , ZLATRS )
# Skipping MacroDefinition: lapackf77_zlauum FORTRAN_NAME ( zlauum , ZLAUUM )
# Skipping MacroDefinition: lapackf77_zlavhe FORTRAN_NAME ( zlavhe , ZLAVHE )
# Skipping MacroDefinition: lapackf77_zposv FORTRAN_NAME ( zposv , ZPOSV )
# Skipping MacroDefinition: lapackf77_zpotrf FORTRAN_NAME ( zpotrf , ZPOTRF )
# Skipping MacroDefinition: lapackf77_zpotri FORTRAN_NAME ( zpotri , ZPOTRI )
# Skipping MacroDefinition: lapackf77_zpotrs FORTRAN_NAME ( zpotrs , ZPOTRS )
# Skipping MacroDefinition: lapackf77_zstedc FORTRAN_NAME ( zstedc , ZSTEDC )
# Skipping MacroDefinition: lapackf77_zstein FORTRAN_NAME ( zstein , ZSTEIN )
# Skipping MacroDefinition: lapackf77_zstemr FORTRAN_NAME ( zstemr , ZSTEMR )
# Skipping MacroDefinition: lapackf77_zsteqr FORTRAN_NAME ( zsteqr , ZSTEQR )
# Skipping MacroDefinition: lapackf77_zsymv FORTRAN_NAME ( zsymv , ZSYMV )
# Skipping MacroDefinition: lapackf77_zsyr FORTRAN_NAME ( zsyr , ZSYR )
# Skipping MacroDefinition: lapackf77_zsysv FORTRAN_NAME ( zsysv , ZSYSV )
# Skipping MacroDefinition: lapackf77_ztrevc FORTRAN_NAME ( ztrevc , ZTREVC )
# Skipping MacroDefinition: lapackf77_ztrevc3 FORTRAN_NAME ( ztrevc3 , ZTREVC3 )
# Skipping MacroDefinition: lapackf77_ztrtri FORTRAN_NAME ( ztrtri , ZTRTRI )
# Skipping MacroDefinition: lapackf77_zung2r FORTRAN_NAME ( zung2r , ZUNG2R )
# Skipping MacroDefinition: lapackf77_zungbr FORTRAN_NAME ( zungbr , ZUNGBR )
# Skipping MacroDefinition: lapackf77_zunghr FORTRAN_NAME ( zunghr , ZUNGHR )
# Skipping MacroDefinition: lapackf77_zunglq FORTRAN_NAME ( zunglq , ZUNGLQ )
# Skipping MacroDefinition: lapackf77_zungql FORTRAN_NAME ( zungql , ZUNGQL )
# Skipping MacroDefinition: lapackf77_zungqr FORTRAN_NAME ( zungqr , ZUNGQR )
# Skipping MacroDefinition: lapackf77_zungtr FORTRAN_NAME ( zungtr , ZUNGTR )
# Skipping MacroDefinition: lapackf77_zunm2r FORTRAN_NAME ( zunm2r , ZUNM2R )
# Skipping MacroDefinition: lapackf77_zunmbr FORTRAN_NAME ( zunmbr , ZUNMBR )
# Skipping MacroDefinition: lapackf77_zunmlq FORTRAN_NAME ( zunmlq , ZUNMLQ )
# Skipping MacroDefinition: lapackf77_zunmql FORTRAN_NAME ( zunmql , ZUNMQL )
# Skipping MacroDefinition: lapackf77_zunmqr FORTRAN_NAME ( zunmqr , ZUNMQR )
# Skipping MacroDefinition: lapackf77_zunmrq FORTRAN_NAME ( zunmrq , ZUNMRQ )
# Skipping MacroDefinition: lapackf77_zunmtr FORTRAN_NAME ( zunmtr , ZUNMTR )
# Skipping MacroDefinition: lapackf77_zbdt01 FORTRAN_NAME ( zbdt01 , ZBDT01 )
# Skipping MacroDefinition: lapackf77_zget22 FORTRAN_NAME ( zget22 , ZGET22 )
# Skipping MacroDefinition: lapackf77_zhet21 FORTRAN_NAME ( zhet21 , ZHET21 )
# Skipping MacroDefinition: lapackf77_zhst01 FORTRAN_NAME ( zhst01 , ZHST01 )
# Skipping MacroDefinition: lapackf77_zlarfy FORTRAN_NAME ( zlarfy , ZLARFY )
# Skipping MacroDefinition: lapackf77_zlatms FORTRAN_NAME ( zlatms , ZLATMS )
# Skipping MacroDefinition: lapackf77_zqpt01 FORTRAN_NAME ( zqpt01 , ZQPT01 )
# Skipping MacroDefinition: lapackf77_zqrt02 FORTRAN_NAME ( zqrt02 , ZQRT02 )
# Skipping MacroDefinition: lapackf77_zstt21 FORTRAN_NAME ( zstt21 , ZSTT21 )
# Skipping MacroDefinition: lapackf77_zunt01 FORTRAN_NAME ( zunt01 , ZUNT01 )
# Skipping MacroDefinition: magma_csetvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_csetvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_cgetvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_cgetvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ccopyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_ccopyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_csetvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_csetvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_cgetvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_cgetvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ccopyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_ccopyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_csetmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_csetmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_cgetmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_cgetmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ccopymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_ccopymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_csetmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_csetmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_cgetmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_cgetmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ccopymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_ccopymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_csetvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_csetvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_cgetvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_cgetvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ccopyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_ccopyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_csetmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_csetmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_cgetmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_cgetmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ccopymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_ccopymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )

const magmablas_ctranspose_inplace = magmablas_ctranspose_inplace_v1
const magmablas_ctranspose_conj_inplace = magmablas_ctranspose_conj_inplace_v1
const magmablas_ctranspose = magmablas_ctranspose_v1
const magmablas_ctranspose_conj = magmablas_ctranspose_conj_v1
const magmablas_cgetmatrix_transpose = magmablas_cgetmatrix_transpose_v1
const magmablas_csetmatrix_transpose = magmablas_csetmatrix_transpose_v1
const magmablas_cprbt = magmablas_cprbt_v1
const magmablas_cprbt_mv = magmablas_cprbt_mv_v1
const magmablas_cprbt_mtv = magmablas_cprbt_mtv_v1
const magma_cgetmatrix_1D_col_bcyclic = magma_cgetmatrix_1D_col_bcyclic_v1
const magma_csetmatrix_1D_col_bcyclic = magma_csetmatrix_1D_col_bcyclic_v1
const magma_cgetmatrix_1D_row_bcyclic = magma_cgetmatrix_1D_row_bcyclic_v1
const magma_csetmatrix_1D_row_bcyclic = magma_csetmatrix_1D_row_bcyclic_v1
const magmablas_cgeadd = magmablas_cgeadd_v1
const magmablas_cgeadd2 = magmablas_cgeadd2_v1
const magmablas_clacpy = magmablas_clacpy_v1
const magmablas_clacpy_conj = magmablas_clacpy_conj_v1
const magmablas_clacpy_sym_in = magmablas_clacpy_sym_in_v1
const magmablas_clacpy_sym_out = magmablas_clacpy_sym_out_v1
const magmablas_clange = magmablas_clange_v1
const magmablas_clanhe = magmablas_clanhe_v1
const magmablas_clansy = magmablas_clansy_v1
const magmablas_clarfg = magmablas_clarfg_v1
const magmablas_clascl = magmablas_clascl_v1
const magmablas_clascl_2x2 = magmablas_clascl_2x2_v1
const magmablas_clascl2 = magmablas_clascl2_v1
const magmablas_clascl_diag = magmablas_clascl_diag_v1
const magmablas_claset = magmablas_claset_v1
const magmablas_claset_band = magmablas_claset_band_v1
const magmablas_claswp = magmablas_claswp_v1
const magmablas_claswp2 = magmablas_claswp2_v1
const magmablas_claswp_sym = magmablas_claswp_sym_v1
const magmablas_claswpx = magmablas_claswpx_v1
const magmablas_csymmetrize = magmablas_csymmetrize_v1
const magmablas_csymmetrize_tiles = magmablas_csymmetrize_tiles_v1
const magmablas_ctrtri_diag = magmablas_ctrtri_diag_v1
const magmablas_scnrm2_adjust = magmablas_scnrm2_adjust_v1
const magmablas_scnrm2_check = magmablas_scnrm2_check_v1
const magmablas_scnrm2_cols = magmablas_scnrm2_cols_v1
const magmablas_scnrm2_row_check_adjust = magmablas_scnrm2_row_check_adjust_v1
const magma_clarfb_gpu = magma_clarfb_gpu_v1
const magma_clarfb_gpu_gemm = magma_clarfb_gpu_gemm_v1
const magma_clarfbx_gpu = magma_clarfbx_gpu_v1
const magma_clarfg_gpu = magma_clarfg_gpu_v1
const magma_clarfgtx_gpu = magma_clarfgtx_gpu_v1
const magma_clarfgx_gpu = magma_clarfgx_gpu_v1
const magma_clarfx_gpu = magma_clarfx_gpu_v1
const magmablas_caxpycp = magmablas_caxpycp_v1
const magmablas_cswap = magmablas_cswap_v1
const magmablas_cswapblk = magmablas_cswapblk_v1
const magmablas_cswapdblk = magmablas_cswapdblk_v1
const magmablas_cgemv = magmablas_cgemv_v1
const magmablas_cgemv_conj = magmablas_cgemv_conj_v1
const magmablas_chemv = magmablas_chemv_v1
const magmablas_csymv = magmablas_csymv_v1
const magmablas_cgemm = magmablas_cgemm_v1
const magmablas_cgemm_reduce = magmablas_cgemm_reduce_v1
const magmablas_chemm = magmablas_chemm_v1
const magmablas_csymm = magmablas_csymm_v1
const magmablas_csyr2k = magmablas_csyr2k_v1
const magmablas_cher2k = magmablas_cher2k_v1
const magmablas_csyrk = magmablas_csyrk_v1
const magmablas_cherk = magmablas_cherk_v1
const magmablas_ctrsm = magmablas_ctrsm_v1
const magmablas_ctrsm_outofplace = magmablas_ctrsm_outofplace_v1
const magmablas_ctrsm_work = magmablas_ctrsm_work_v1
const magma_csetvector = magma_csetvector_v1
const magma_cgetvector = magma_cgetvector_v1
const magma_ccopyvector = magma_ccopyvector_v1
const magma_csetmatrix = magma_csetmatrix_v1
const magma_cgetmatrix = magma_cgetmatrix_v1
const magma_ccopymatrix = magma_ccopymatrix_v1
const magma_icamax = magma_icamax_v1
const magma_icamin = magma_icamin_v1
const magma_scasum = magma_scasum_v1
const magma_caxpy = magma_caxpy_v1
const magma_ccopy = magma_ccopy_v1
const magma_cdotc = magma_cdotc_v1
const magma_cdotu = magma_cdotu_v1
const magma_scnrm2 = magma_scnrm2_v1
const magma_crot = magma_crot_v1
const magma_csrot = magma_csrot_v1
const magma_crotm = magma_crotm_v1
const magma_crotmg = magma_crotmg_v1
const magma_cscal = magma_cscal_v1
const magma_csscal = magma_csscal_v1
const magma_cswap = magma_cswap_v1
const magma_cgemv = magma_cgemv_v1
const magma_cgerc = magma_cgerc_v1
const magma_cgeru = magma_cgeru_v1
const magma_chemv = magma_chemv_v1
const magma_cher = magma_cher_v1
const magma_cher2 = magma_cher2_v1
const magma_ctrmv = magma_ctrmv_v1
const magma_ctrsv = magma_ctrsv_v1
const magma_cgemm = magma_cgemm_v1
const magma_csymm = magma_csymm_v1
const magma_chemm = magma_chemm_v1
const magma_csyr2k = magma_csyr2k_v1
const magma_cher2k = magma_cher2k_v1
const magma_csyrk = magma_csyrk_v1
const magma_cherk = magma_cherk_v1
const magma_ctrmm = magma_ctrmm_v1
const magma_ctrsm = magma_ctrsm_v1

# Skipping MacroDefinition: magma_dsetvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_dsetvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dgetvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_dgetvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dcopyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_dcopyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dsetvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_dsetvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dgetvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_dgetvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dcopyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_dcopyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dsetmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_dsetmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dgetmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_dgetmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dcopymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_dcopymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dsetmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_dsetmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dgetmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_dgetmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dcopymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_dcopymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dsetvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_dsetvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dgetvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_dgetvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dcopyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_dcopyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dsetmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_dsetmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dgetmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_dgetmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_dcopymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_dcopymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )

const magmablas_dtranspose_inplace = magmablas_dtranspose_inplace_v1
const magmablas_dtranspose = magmablas_dtranspose_v1
const magmablas_dgetmatrix_transpose = magmablas_dgetmatrix_transpose_v1
const magmablas_dsetmatrix_transpose = magmablas_dsetmatrix_transpose_v1
const magmablas_dprbt = magmablas_dprbt_v1
const magmablas_dprbt_mv = magmablas_dprbt_mv_v1
const magmablas_dprbt_mtv = magmablas_dprbt_mtv_v1
const magma_dgetmatrix_1D_col_bcyclic = magma_dgetmatrix_1D_col_bcyclic_v1
const magma_dsetmatrix_1D_col_bcyclic = magma_dsetmatrix_1D_col_bcyclic_v1
const magma_dgetmatrix_1D_row_bcyclic = magma_dgetmatrix_1D_row_bcyclic_v1
const magma_dsetmatrix_1D_row_bcyclic = magma_dsetmatrix_1D_row_bcyclic_v1
const magmablas_dgeadd = magmablas_dgeadd_v1
const magmablas_dgeadd2 = magmablas_dgeadd2_v1
const magmablas_dlacpy = magmablas_dlacpy_v1
const magmablas_dlacpy_conj = magmablas_dlacpy_conj_v1
const magmablas_dlacpy_sym_in = magmablas_dlacpy_sym_in_v1
const magmablas_dlacpy_sym_out = magmablas_dlacpy_sym_out_v1
const magmablas_dlange = magmablas_dlange_v1
const magmablas_dlansy = magmablas_dlansy_v1
const magmablas_dlarfg = magmablas_dlarfg_v1
const magmablas_dlascl = magmablas_dlascl_v1
const magmablas_dlascl_2x2 = magmablas_dlascl_2x2_v1
const magmablas_dlascl2 = magmablas_dlascl2_v1
const magmablas_dlascl_diag = magmablas_dlascl_diag_v1
const magmablas_dlaset = magmablas_dlaset_v1
const magmablas_dlaset_band = magmablas_dlaset_band_v1
const magmablas_dlaswp = magmablas_dlaswp_v1
const magmablas_dlaswp2 = magmablas_dlaswp2_v1
const magmablas_dlaswp_sym = magmablas_dlaswp_sym_v1
const magmablas_dlaswpx = magmablas_dlaswpx_v1
const magmablas_dsymmetrize = magmablas_dsymmetrize_v1
const magmablas_dsymmetrize_tiles = magmablas_dsymmetrize_tiles_v1
const magmablas_dtrtri_diag = magmablas_dtrtri_diag_v1
const magmablas_dnrm2_adjust = magmablas_dnrm2_adjust_v1
const magmablas_dnrm2_check = magmablas_dnrm2_check_v1
const magmablas_dnrm2_cols = magmablas_dnrm2_cols_v1
const magmablas_dnrm2_row_check_adjust = magmablas_dnrm2_row_check_adjust_v1
const magma_dlarfb_gpu = magma_dlarfb_gpu_v1
const magma_dlarfb_gpu_gemm = magma_dlarfb_gpu_gemm_v1
const magma_dlarfbx_gpu = magma_dlarfbx_gpu_v1
const magma_dlarfg_gpu = magma_dlarfg_gpu_v1
const magma_dlarfgtx_gpu = magma_dlarfgtx_gpu_v1
const magma_dlarfgx_gpu = magma_dlarfgx_gpu_v1
const magma_dlarfx_gpu = magma_dlarfx_gpu_v1
const magmablas_daxpycp = magmablas_daxpycp_v1
const magmablas_dswap = magmablas_dswap_v1
const magmablas_dswapblk = magmablas_dswapblk_v1
const magmablas_dswapdblk = magmablas_dswapdblk_v1
const magmablas_dgemv = magmablas_dgemv_v1
const magmablas_dgemv_conj = magmablas_dgemv_conj_v1
const magmablas_dsymv = magmablas_dsymv_v1
const magmablas_dgemm = magmablas_dgemm_v1
const magmablas_dgemm_reduce = magmablas_dgemm_reduce_v1
const magmablas_dsymm = magmablas_dsymm_v1
const magmablas_dsyr2k = magmablas_dsyr2k_v1
const magmablas_dsyrk = magmablas_dsyrk_v1
const magmablas_dtrsm = magmablas_dtrsm_v1
const magmablas_dtrsm_outofplace = magmablas_dtrsm_outofplace_v1
const magmablas_dtrsm_work = magmablas_dtrsm_work_v1
const magma_dsetvector = magma_dsetvector_v1
const magma_dgetvector = magma_dgetvector_v1
const magma_dcopyvector = magma_dcopyvector_v1
const magma_dsetmatrix = magma_dsetmatrix_v1
const magma_dgetmatrix = magma_dgetmatrix_v1
const magma_dcopymatrix = magma_dcopymatrix_v1
const magma_idamax = magma_idamax_v1
const magma_idamin = magma_idamin_v1
const magma_dasum = magma_dasum_v1
const magma_daxpy = magma_daxpy_v1
const magma_dcopy = magma_dcopy_v1
const magma_ddot = magma_ddot_v1
const magma_dnrm2 = magma_dnrm2_v1
const magma_drot = magma_drot_v1
const magma_drotm = magma_drotm_v1
const magma_drotmg = magma_drotmg_v1
const magma_dscal = magma_dscal_v1
const magma_dswap = magma_dswap_v1
const magma_dgemv = magma_dgemv_v1
const magma_dger = magma_dger_v1
const magma_dsymv = magma_dsymv_v1
const magma_dsyr = magma_dsyr_v1
const magma_dsyr2 = magma_dsyr2_v1
const magma_dtrmv = magma_dtrmv_v1
const magma_dtrsv = magma_dtrsv_v1
const magma_dgemm = magma_dgemm_v1
const magma_dsymm = magma_dsymm_v1
const magma_dsyr2k = magma_dsyr2k_v1
const magma_dsyrk = magma_dsyrk_v1
const magma_dtrmm = magma_dtrmm_v1
const magma_dtrsm = magma_dtrsm_v1
const magmablas_dsaxpycp = magmablas_dsaxpycp_v1
const magmablas_dslaswp = magmablas_dslaswp_v1
const magmablas_dlag2s = magmablas_dlag2s_v1
const magmablas_slag2d = magmablas_slag2d_v1
const magmablas_dlat2s = magmablas_dlat2s_v1
const magmablas_slat2d = magmablas_slat2d_v1

# Skipping MacroDefinition: magma_ssetvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_ssetvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_sgetvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_sgetvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_scopyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_scopyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ssetvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_ssetvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_sgetvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_sgetvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_scopyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_scopyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ssetmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_ssetmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_sgetmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_sgetmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_scopymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_scopymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ssetmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_ssetmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_sgetmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_sgetmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_scopymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_scopymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ssetvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_ssetvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_sgetvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_sgetvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_scopyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_scopyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_ssetmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_ssetmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_sgetmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_sgetmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_scopymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_scopymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )

const magmablas_stranspose_inplace = magmablas_stranspose_inplace_v1
const magmablas_stranspose = magmablas_stranspose_v1
const magmablas_sgetmatrix_transpose = magmablas_sgetmatrix_transpose_v1
const magmablas_ssetmatrix_transpose = magmablas_ssetmatrix_transpose_v1
const magmablas_sprbt = magmablas_sprbt_v1
const magmablas_sprbt_mv = magmablas_sprbt_mv_v1
const magmablas_sprbt_mtv = magmablas_sprbt_mtv_v1
const magma_sgetmatrix_1D_col_bcyclic = magma_sgetmatrix_1D_col_bcyclic_v1
const magma_ssetmatrix_1D_col_bcyclic = magma_ssetmatrix_1D_col_bcyclic_v1
const magma_sgetmatrix_1D_row_bcyclic = magma_sgetmatrix_1D_row_bcyclic_v1
const magma_ssetmatrix_1D_row_bcyclic = magma_ssetmatrix_1D_row_bcyclic_v1
const magmablas_sgeadd = magmablas_sgeadd_v1
const magmablas_sgeadd2 = magmablas_sgeadd2_v1
const magmablas_slacpy = magmablas_slacpy_v1
const magmablas_slacpy_conj = magmablas_slacpy_conj_v1
const magmablas_slacpy_sym_in = magmablas_slacpy_sym_in_v1
const magmablas_slacpy_sym_out = magmablas_slacpy_sym_out_v1
const magmablas_slange = magmablas_slange_v1
const magmablas_slansy = magmablas_slansy_v1
const magmablas_slarfg = magmablas_slarfg_v1
const magmablas_slascl = magmablas_slascl_v1
const magmablas_slascl_2x2 = magmablas_slascl_2x2_v1
const magmablas_slascl2 = magmablas_slascl2_v1
const magmablas_slascl_diag = magmablas_slascl_diag_v1
const magmablas_slaset = magmablas_slaset_v1
const magmablas_slaset_band = magmablas_slaset_band_v1
const magmablas_slaswp = magmablas_slaswp_v1
const magmablas_slaswp2 = magmablas_slaswp2_v1
const magmablas_slaswp_sym = magmablas_slaswp_sym_v1
const magmablas_slaswpx = magmablas_slaswpx_v1
const magmablas_ssymmetrize = magmablas_ssymmetrize_v1
const magmablas_ssymmetrize_tiles = magmablas_ssymmetrize_tiles_v1
const magmablas_strtri_diag = magmablas_strtri_diag_v1
const magmablas_snrm2_adjust = magmablas_snrm2_adjust_v1
const magmablas_snrm2_check = magmablas_snrm2_check_v1
const magmablas_snrm2_cols = magmablas_snrm2_cols_v1
const magmablas_snrm2_row_check_adjust = magmablas_snrm2_row_check_adjust_v1
const magma_slarfb_gpu = magma_slarfb_gpu_v1
const magma_slarfb_gpu_gemm = magma_slarfb_gpu_gemm_v1
const magma_slarfbx_gpu = magma_slarfbx_gpu_v1
const magma_slarfg_gpu = magma_slarfg_gpu_v1
const magma_slarfgtx_gpu = magma_slarfgtx_gpu_v1
const magma_slarfgx_gpu = magma_slarfgx_gpu_v1
const magma_slarfx_gpu = magma_slarfx_gpu_v1
const magmablas_saxpycp = magmablas_saxpycp_v1
const magmablas_sswap = magmablas_sswap_v1
const magmablas_sswapblk = magmablas_sswapblk_v1
const magmablas_sswapdblk = magmablas_sswapdblk_v1
const magmablas_sgemv = magmablas_sgemv_v1
const magmablas_sgemv_conj = magmablas_sgemv_conj_v1
const magmablas_ssymv = magmablas_ssymv_v1
const magmablas_sgemm = magmablas_sgemm_v1
const magmablas_sgemm_reduce = magmablas_sgemm_reduce_v1
const magmablas_ssymm = magmablas_ssymm_v1
const magmablas_ssyr2k = magmablas_ssyr2k_v1
const magmablas_ssyrk = magmablas_ssyrk_v1
const magmablas_strsm = magmablas_strsm_v1
const magmablas_strsm_outofplace = magmablas_strsm_outofplace_v1
const magmablas_strsm_work = magmablas_strsm_work_v1
const magma_ssetvector = magma_ssetvector_v1
const magma_sgetvector = magma_sgetvector_v1
const magma_scopyvector = magma_scopyvector_v1
const magma_ssetmatrix = magma_ssetmatrix_v1
const magma_sgetmatrix = magma_sgetmatrix_v1
const magma_scopymatrix = magma_scopymatrix_v1
const magma_isamax = magma_isamax_v1
const magma_isamin = magma_isamin_v1
const magma_sasum = magma_sasum_v1
const magma_saxpy = magma_saxpy_v1
const magma_scopy = magma_scopy_v1
const magma_sdot = magma_sdot_v1
const magma_snrm2 = magma_snrm2_v1
const magma_srot = magma_srot_v1
const magma_srotm = magma_srotm_v1
const magma_srotmg = magma_srotmg_v1
const magma_sscal = magma_sscal_v1
const magma_sswap = magma_sswap_v1
const magma_sgemv = magma_sgemv_v1
const magma_sger = magma_sger_v1
const magma_ssymv = magma_ssymv_v1
const magma_ssyr = magma_ssyr_v1
const magma_ssyr2 = magma_ssyr2_v1
const magma_strmv = magma_strmv_v1
const magma_strsv = magma_strsv_v1
const magma_sgemm = magma_sgemm_v1
const magma_ssymm = magma_ssymm_v1
const magma_ssyr2k = magma_ssyr2k_v1
const magma_ssyrk = magma_ssyrk_v1
const magma_strmm = magma_strmm_v1
const magma_strsm = magma_strsm_v1

# Skipping MacroDefinition: magma_queue_create_v1 ( queue_ptr ) magma_queue_create_v1_internal ( queue_ptr , __func__ , __FILE__ , __LINE__ )

const MagmaUpperLower = MagmaFull
const MagmaUpperLowerStr = MagmaFullStr

# Skipping MacroDefinition: MAGMA_Z_CNJG ( a ) MAGMA_Z_CONJ ( a )
# Skipping MacroDefinition: MAGMA_C_CNJG ( a ) MAGMA_C_CONJ ( a )
# Skipping MacroDefinition: MAGMA_D_CNJG ( a ) MAGMA_D_CONJ ( a )
# Skipping MacroDefinition: MAGMA_S_CNJG ( a ) MAGMA_S_CONJ ( a )

const magma_queue_create = magma_queue_create_v1
const magma_setvector = magma_setvector_v1
const magma_getvector = magma_getvector_v1
const magma_copyvector = magma_copyvector_v1
const magma_setmatrix = magma_setmatrix_v1
const magma_getmatrix = magma_getmatrix_v1
const magma_copymatrix = magma_copymatrix_v1
const magma_isetvector = magma_isetvector_v1
const magma_igetvector = magma_igetvector_v1
const magma_icopyvector = magma_icopyvector_v1
const magma_isetmatrix = magma_isetmatrix_v1
const magma_igetmatrix = magma_igetmatrix_v1
const magma_icopymatrix = magma_icopymatrix_v1
const magma_index_setvector = magma_index_setvector_v1
const magma_index_getvector = magma_index_getvector_v1
const magma_index_copyvector = magma_index_copyvector_v1
const magma_index_setmatrix = magma_index_setmatrix_v1
const magma_index_getmatrix = magma_index_getmatrix_v1
const magma_index_copymatrix = magma_index_copymatrix_v1

# Skipping MacroDefinition: magma_zsetvector ( n , hx_src , incx , dy_dst , incy , queue ) magma_zsetvector_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zgetvector ( n , dx_src , incx , hy_dst , incy , queue ) magma_zgetvector_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zcopyvector ( n , dx_src , incx , dy_dst , incy , queue ) magma_zcopyvector_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zsetvector_async ( n , hx_src , incx , dy_dst , incy , queue ) magma_zsetvector_async_internal ( n , hx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zgetvector_async ( n , dx_src , incx , hy_dst , incy , queue ) magma_zgetvector_async_internal ( n , dx_src , incx , hy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zcopyvector_async ( n , dx_src , incx , dy_dst , incy , queue ) magma_zcopyvector_async_internal ( n , dx_src , incx , dy_dst , incy , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zsetmatrix ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_zsetmatrix_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zgetmatrix ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_zgetmatrix_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zcopymatrix ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_zcopymatrix_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zsetmatrix_async ( m , n , hA_src , lda , dB_dst , lddb , queue ) magma_zsetmatrix_async_internal ( m , n , hA_src , lda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zgetmatrix_async ( m , n , dA_src , ldda , hB_dst , ldb , queue ) magma_zgetmatrix_async_internal ( m , n , dA_src , ldda , hB_dst , ldb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zcopymatrix_async ( m , n , dA_src , ldda , dB_dst , lddb , queue ) magma_zcopymatrix_async_internal ( m , n , dA_src , ldda , dB_dst , lddb , queue , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zsetvector_v1 ( n , hx_src , incx , dy_dst , incy ) magma_zsetvector_v1_internal ( n , hx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zgetvector_v1 ( n , dx_src , incx , hy_dst , incy ) magma_zgetvector_v1_internal ( n , dx_src , incx , hy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zcopyvector_v1 ( n , dx_src , incx , dy_dst , incy ) magma_zcopyvector_v1_internal ( n , dx_src , incx , dy_dst , incy , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zsetmatrix_v1 ( m , n , hA_src , lda , dB_dst , lddb ) magma_zsetmatrix_v1_internal ( m , n , hA_src , lda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zgetmatrix_v1 ( m , n , dA_src , ldda , hB_dst , ldb ) magma_zgetmatrix_v1_internal ( m , n , dA_src , ldda , hB_dst , ldb , __func__ , __FILE__ , __LINE__ )
# Skipping MacroDefinition: magma_zcopymatrix_v1 ( m , n , dA_src , ldda , dB_dst , lddb ) magma_zcopymatrix_v1_internal ( m , n , dA_src , ldda , dB_dst , lddb , __func__ , __FILE__ , __LINE__ )

const magmablas_ztranspose_inplace = magmablas_ztranspose_inplace_v1
const magmablas_ztranspose_conj_inplace = magmablas_ztranspose_conj_inplace_v1
const magmablas_ztranspose = magmablas_ztranspose_v1
const magmablas_ztranspose_conj = magmablas_ztranspose_conj_v1
const magmablas_zgetmatrix_transpose = magmablas_zgetmatrix_transpose_v1
const magmablas_zsetmatrix_transpose = magmablas_zsetmatrix_transpose_v1
const magmablas_zprbt = magmablas_zprbt_v1
const magmablas_zprbt_mv = magmablas_zprbt_mv_v1
const magmablas_zprbt_mtv = magmablas_zprbt_mtv_v1
const magma_zgetmatrix_1D_col_bcyclic = magma_zgetmatrix_1D_col_bcyclic_v1
const magma_zsetmatrix_1D_col_bcyclic = magma_zsetmatrix_1D_col_bcyclic_v1
const magma_zgetmatrix_1D_row_bcyclic = magma_zgetmatrix_1D_row_bcyclic_v1
const magma_zsetmatrix_1D_row_bcyclic = magma_zsetmatrix_1D_row_bcyclic_v1
const magmablas_zgeadd = magmablas_zgeadd_v1
const magmablas_zgeadd2 = magmablas_zgeadd2_v1
const magmablas_zlacpy = magmablas_zlacpy_v1
const magmablas_zlacpy_conj = magmablas_zlacpy_conj_v1
const magmablas_zlacpy_sym_in = magmablas_zlacpy_sym_in_v1
const magmablas_zlacpy_sym_out = magmablas_zlacpy_sym_out_v1
const magmablas_zlange = magmablas_zlange_v1
const magmablas_zlanhe = magmablas_zlanhe_v1
const magmablas_zlansy = magmablas_zlansy_v1
const magmablas_zlarfg = magmablas_zlarfg_v1
const magmablas_zlascl = magmablas_zlascl_v1
const magmablas_zlascl_2x2 = magmablas_zlascl_2x2_v1
const magmablas_zlascl2 = magmablas_zlascl2_v1
const magmablas_zlascl_diag = magmablas_zlascl_diag_v1
const magmablas_zlaset = magmablas_zlaset_v1
const magmablas_zlaset_band = magmablas_zlaset_band_v1
const magmablas_zlaswp = magmablas_zlaswp_v1
const magmablas_zlaswp2 = magmablas_zlaswp2_v1
const magmablas_zlaswp_sym = magmablas_zlaswp_sym_v1
const magmablas_zlaswpx = magmablas_zlaswpx_v1
const magmablas_zsymmetrize = magmablas_zsymmetrize_v1
const magmablas_zsymmetrize_tiles = magmablas_zsymmetrize_tiles_v1
const magmablas_ztrtri_diag = magmablas_ztrtri_diag_v1
const magmablas_dznrm2_adjust = magmablas_dznrm2_adjust_v1
const magmablas_dznrm2_check = magmablas_dznrm2_check_v1
const magmablas_dznrm2_cols = magmablas_dznrm2_cols_v1
const magmablas_dznrm2_row_check_adjust = magmablas_dznrm2_row_check_adjust_v1
const magma_zlarfb_gpu = magma_zlarfb_gpu_v1
const magma_zlarfb_gpu_gemm = magma_zlarfb_gpu_gemm_v1
const magma_zlarfbx_gpu = magma_zlarfbx_gpu_v1
const magma_zlarfg_gpu = magma_zlarfg_gpu_v1
const magma_zlarfgtx_gpu = magma_zlarfgtx_gpu_v1
const magma_zlarfgx_gpu = magma_zlarfgx_gpu_v1
const magma_zlarfx_gpu = magma_zlarfx_gpu_v1
const magmablas_zaxpycp = magmablas_zaxpycp_v1
const magmablas_zswap = magmablas_zswap_v1
const magmablas_zswapblk = magmablas_zswapblk_v1
const magmablas_zswapdblk = magmablas_zswapdblk_v1
const magmablas_zgemv = magmablas_zgemv_v1
const magmablas_zgemv_conj = magmablas_zgemv_conj_v1
const magmablas_zhemv = magmablas_zhemv_v1
const magmablas_zsymv = magmablas_zsymv_v1
const magmablas_zgemm = magmablas_zgemm_v1
const magmablas_zgemm_reduce = magmablas_zgemm_reduce_v1
const magmablas_zhemm = magmablas_zhemm_v1
const magmablas_zsymm = magmablas_zsymm_v1
const magmablas_zsyr2k = magmablas_zsyr2k_v1
const magmablas_zher2k = magmablas_zher2k_v1
const magmablas_zsyrk = magmablas_zsyrk_v1
const magmablas_zherk = magmablas_zherk_v1
const magmablas_ztrsm = magmablas_ztrsm_v1
const magmablas_ztrsm_outofplace = magmablas_ztrsm_outofplace_v1
const magmablas_ztrsm_work = magmablas_ztrsm_work_v1
const magma_zsetvector = magma_zsetvector_v1
const magma_zgetvector = magma_zgetvector_v1
const magma_zcopyvector = magma_zcopyvector_v1
const magma_zsetmatrix = magma_zsetmatrix_v1
const magma_zgetmatrix = magma_zgetmatrix_v1
const magma_zcopymatrix = magma_zcopymatrix_v1
const magma_izamax = magma_izamax_v1
const magma_izamin = magma_izamin_v1
const magma_dzasum = magma_dzasum_v1
const magma_zaxpy = magma_zaxpy_v1
const magma_zcopy = magma_zcopy_v1
const magma_zdotc = magma_zdotc_v1
const magma_zdotu = magma_zdotu_v1
const magma_dznrm2 = magma_dznrm2_v1
const magma_zrot = magma_zrot_v1
const magma_zdrot = magma_zdrot_v1
const magma_zrotm = magma_zrotm_v1
const magma_zrotmg = magma_zrotmg_v1
const magma_zscal = magma_zscal_v1
const magma_zdscal = magma_zdscal_v1
const magma_zswap = magma_zswap_v1
const magma_zgemv = magma_zgemv_v1
const magma_zgerc = magma_zgerc_v1
const magma_zgeru = magma_zgeru_v1
const magma_zhemv = magma_zhemv_v1
const magma_zher = magma_zher_v1
const magma_zher2 = magma_zher2_v1
const magma_ztrmv = magma_ztrmv_v1
const magma_ztrsv = magma_ztrsv_v1
const magma_zgemm = magma_zgemm_v1
const magma_zsymm = magma_zsymm_v1
const magma_zhemm = magma_zhemm_v1
const magma_zsyr2k = magma_zsyr2k_v1
const magma_zher2k = magma_zher2k_v1
const magma_zsyrk = magma_zsyrk_v1
const magma_zherk = magma_zherk_v1
const magma_ztrmm = magma_ztrmm_v1
const magma_ztrsm = magma_ztrsm_v1
const magmablas_zcaxpycp = magmablas_zcaxpycp_v1
const magmablas_zclaswp = magmablas_zclaswp_v1
const magmablas_zlag2c = magmablas_zlag2c_v1
const magmablas_clag2z = magmablas_clag2z_v1
const magmablas_zlat2c = magmablas_zlat2c_v1
const magmablas_clat2z = magmablas_clat2z_v1
