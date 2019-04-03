# These are transcribed from magma_types.h

const MagmaFalse         = 0
const MagmaTrue          = 1

const MagmaRowMajor      = 101
const MagmaColMajor      = 102

const MagmaNoTrans       = 111
const MagmaTrans         = 112
const MagmaConjTrans     = 113

const MagmaUpper         = 121
const MagmaLower         = 122
const MagmaUpperLower    = 123
const MagmaFull          = 123

const MagmaNonUnit       = 131
const MagmaUnit          = 132

const MagmaLeft          = 141
const MagmaRight         = 142
const MagmaBothSides     = 143

const MagmaOneNorm       = 171
const MagmaRealOneNorm   = 172
const MagmaTwoNorm       = 173
const MagmaFrobeniusNorm = 174
const MagmaInfNorm       = 175
const MagmaRealInfNorm   = 176
const MagmaMaxNorm       = 177
const MagmaRealMaxNorm   = 178

const MagmaDistUniform   = 201
const MagmaDistSymmetric = 202
const MagmaDistNormal    = 203

const MagmaHermGeev      = 241
const MagmaHermPoev      = 242
const MagmaNonsymPosv    = 243
const MagmaSymPosv       = 244

const MagmaNoPacking     = 291
const MagmaPackSubdiag   = 292
const MagmaPackSupdiag   = 293
const MagmaPackColumn    = 294
const MagmaPackRow       = 295
const MagmaPackLowerBand = 296
const MagmaPackUpeprBand = 297
const MagmaPackAll       = 298

# MAGMA constants indicating the vectors status as input/output for some functions
# For example, the gesvd functions will use MagmaNoVec, MagmaSomeVec, MagmaAllVec and MagmaOverwriteVec to indicate the
# strategies that will be applied to the SVD U matrix and VT matrix in A = U Σ V**T (for MagmaOverwriteVec it is going to overwrite A)
const MagmaNoVec         = 301
const MagmaVec           = 302
const MagmaIVec          = 303
const MagmaAllVec        = 304
const MagmaSomeVec       = 305
const MagmaOverwriteVec  = 306
const MagmaBacktransVec  = 307

const MagmaRangeAll      = 311
const MagmaRangeV        = 312
const MagmaRangeI        = 313

const MagmaQ             = 322
const MagmaP             = 323

const MagmaForward       = 391
const MagmaBackward      = 392

const MagmaColumnwise    = 401
const MagmaRowwise       = 402

const Magma_CSR          = 411
const Magma_ELLPACK      = 412
const Magma_ELL          = 413
const Magma_DENSE        = 414
const Magma_BCSR         = 415
const Magma_CSC          = 416
const Magma_HYB          = 417
const Magma_COO          = 418
const Magma_ELLRT        = 419
const Magma_SELLC        = 420
const Magma_SELLP        = 421
const Magma_ELLD         = 422
const Magma_ELLDD        = 423
const Magma_CSRD         = 424
const Magma_CSRL         = 427
const Magma_CSRU         = 428
const Magma_CSRCOO       = 429


const Magma_CG           = 431
const Magma_CGMERGE      = 432
const Magma_GMRES        = 433
const Magma_BICGSTAB     = 434
const Magma_BICGSTABMERGE  = 435
const Magma_BICGSTABMERGE2 = 436
const Magma_JACOBI       = 437
const Magma_GS           = 438
const Magma_ITERREF      = 439
const Magma_BCSRLU       = 440
const Magma_PCG          = 441
const Magma_PGMRES       = 442
const Magma_PBICGSTAB    = 443
const Magma_PASTIX       = 444
const Magma_ILU          = 445
const Magma_ICC          = 446
const Magma_AILU         = 447
const Magma_AICC         = 448
const Magma_BAITER       = 449
const Magma_LOBPCG       = 450
const Magma_NONE         = 451

const Magma_CGS          = 461
const Magma_FUSED_CGS    = 462
const Magma_MGS          = 463

const Magma_CPU          = 471
const Magma_DEV          = 472

const Magma_GENERAL      = 481
const Magma_SYMMETRIC    = 482

const Magma_ORDERED      = 491
const Magma_DIAGFIRST    = 492
const Magma_UNITY        = 493
const Magma_VALUE        = 494

const Magma_DCOMPLEX     = 501
const Magma_FCOMPLEX     = 502
const Magma_DOUBLE       = 503
const Magma_FLOAT        = 504

const Magma_NOSCALE      = 511
const Magma_UNITROW      = 512
const Magma_UNITDIAG     = 513

