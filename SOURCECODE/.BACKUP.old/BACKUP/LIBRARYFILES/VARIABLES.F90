        !*************************************************************!
        !*      SUM_HILLS-1.                                         *!
        !*                                                           *!
        !*      A PARALLEL IMPLEMENTATION OF METADYNAMICS FREE       *!
        !*      ENERGY SURFACE CONSTRUCTION.                         *!
        !*                                                           *!
        !*      DEVELOPED BY VENKAT KAPIL                            *!
        !*      CONTACT : venkat@iitk.ac.in                          *!
        !*                                                           *!
        !*      COPYRIGHT (C) 2015 VENKAT KAPIL                      *!
        !*      SOME RIGHTS RESERVED.                                *!
        !*                                                           *!
        !*************************************************************!

	!THE FOLLOWING IS THE LIST IF VARIABLES THAT HAVE BEEN USED.
	!*************************************************************!
	!IMPORTANT GRID VARIABLES.
	!*************************************************************!
	!*) NDIM		:	NUMBER OF DIMENSIONS.
	!*) GRIDMIN		:	MINIMUM OF GRID.
	!*) GRIDMAX		:	MAXIMUM OF GRID.
	!*) GRIDDIF		:	BIN WIDTHF OF THE GRID.
	!*) NBIN		:	NUMBER OF BINS.
	!*) S			:	AN ORDERED N - TUPLE REPRESENTING POSITION OF A POINT IN THE GRID.
	!*) I			:	AN ORDERED N - TUPLE REPRESENTING INDEX OF A POINT IN THE GRID.
	!*) NGRID		:	TOTAL NUMBER OF GRID POINTS.
	!*************************************************************!
	!IMPORTANT HILLS VARIABLES.
	!*************************************************************!
	!*) PACKAGE		:	CODE FOR PACKAGE USED FOR METADYNAMICS; 1 FOR PLUMED; 2 FOR CPMD.
	!*) NMETACOLVAR		:	NUMBER OF COLLECTIVE VARIABLE(S).
	!*) IMETACOLVAR		:	INDEX/INDICES OF COLLECTIVE VARIABLE(S).
	!*) NMETACOLVARGRID	:	NUMBER OF GRID POINTS IN COLLECTIVE VARIABLE SPACE.
	!*) NHILL		:	NUMBER OF COLL. COORDINATES USED.
	!*) CENTER		:	MEAN OF GAUSSIAN.
	!*) HEIGHT		:	MAX. VALUE OF GAUSSIAN.
	!*) WIDTH		:	STANDARD DEVIATION OF GAUSSIAN.
	!*) TIME		:	MD STEP AT WHICH GAUSSIAN WAS ADDED.
	!*) IMIN		:	STARTING INDEX FOR INTEGRATION OF BIAS.
	!*) IMAX		:	ENDING INDEX FOR INTEGRATION BIAS.
	!*) BIAS		:	TOTAL BIAS.
	!*) TEMP		:	TEMPERATURE OF ENSEMBLE.
	!*) WTAPLHA		:	MAGNUTUDE BY WHICH FREE ENERGY IS PROPORTIONAL TO BIAS.
	!*************************************************************!
	!IMPORTANT REWEIGHTING PARAMETERS.
	!*************************************************************!
	!*) C			:	REWEIGHTING FACTOR.
        !*) CNUMER		:	NUMERATOR OF C.
        !*) CDENOM		:	DENOMENATOR OF C.
        !*) CDUMP		:	FREQUENCY AT WHICH C IS PRINTED. MEASURED IN PER HILLS.
	!*) NHILLMIN		:	STARTING POINT FOR REWEIGHTING.
	!*) NHILLMAX		:	END POINT FOR REWEIGHTING.
	!*) GRIDPOINT		:	NEAREST GRID POINT TO COLLECTIVE VARIABLE IN TIME SERIES.
	!*) WEIGHT		:	STATISTICAL WEIGHT OF COLLECTIVE VARIABLE IN TIME SERIES.
	!*) PROB		:	PROBABILITY DISTRIBUTION FUNCTION.
	!*************************************************************!
	!IMPORTANT FILE PARAMETERS.
	!*************************************************************!
	!*) INPUTFILE		:	INPUTFILE.
	!*) LOGFILE		:	LOGFILE.
	!*) OUTPUTFILE1		:	CONTAINS FREE ENERGY SURFACE FROM  ADDING GAUSSIANS.
	!*) OUTPUTFILE2		:	CONTAINS FREE ENERGY SURFACE FROM REWEIGHTING.
	!*) CFILE		:	CONTAINS REWEIGHTING FACTOR AS A FUNCTION OF HILLS.
	!*) ERRORFILE		:	CONTAINS LIST OF ERRORS.
	!*) ERRORSTAT		:	FLAG FOR CHECKING ERRORS.
	!*************************************************************! 

	!*************************************************************! 
	!>>>>>ENVIR VARIABLES<<<<<!
	!*************************************************************! 
	MODULE GLOBAL_VAR
	!GRID VARIABLES.
	INTEGER NDIM
	INTEGER*8 NGRID
	REAL*8, DIMENSION(:), ALLOCATABLE ::  GRIDMIN, GRIDMAX, GRIDDIF,S
	INTEGER, DIMENSION(:), ALLOCATABLE :: NBIN, I
	!HILLS VARIABLES.
	INTEGER NMETACOLVAR, NHILL, PACKAGE, NMETACOLVARGRID
	REAL*8 TEMP, WTDT, WTALPHA
	LOGICAL WELLTEMP
	INTEGER, DIMENSION(:), ALLOCATABLE ::  TIME, IMETACOLVAR
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: IMIN, IMAX
	REAL*8, DIMENSION(:), ALLOCATABLE :: WIDTH, HEIGHT, BIAS
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: CENTER
	!REWEIGHTING.
	INTEGER NHILLMIN, NHILLMAX, CDUMP
	REAL*8 CNUM, CDEN
	LOGICAL REWEIGHT
	REAL*8, DIMENSION(:), ALLOCATABLE :: GRIDPOINT, PROB, C, CNUMER, CDENOM
	!FILES.
	CHARACTER (LEN = 100) INPUTFILE, LOGFILE, OUTPUTFILE1, OUTPUTFILE22, ERRORFILE, CFILE, CVFILE
	LOGICAL ERRORSTAT
	!#IF DEFINED (_PARALLEL)                                                                    
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	!DEFINES VARIABLES FOR MPI 
	INTEGER MPIERROR, MPIBIND, MPIRANK, MPITASKS
	INTEGER, DIMENSION(:), ALLOCATABLE :: MPIBIN
	REAL*8, DIMENSION(:), ALLOCATABLE :: TBIAS, TPROB, TCNUMER, TCDENOM
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	!#ENDIF
	END MODULE
	!*************************************************************! 
	!>>>>>ENVIR VARIABLES<<<<<!
	!*************************************************************! 
