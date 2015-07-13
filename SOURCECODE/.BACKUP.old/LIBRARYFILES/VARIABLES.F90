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

	!*************************************************************!
	!IMPORTANT VARIABLES.
	!*************************************************************!
	!THE FOLLOWING IS THE LIST IF VARIABLES THAT HAVE BEEN USED.
	!*) PACKAGE		:	CODE FOR PACKAGE USED FOR METADYNAMICS; 1 FOR PLUMED; 2 FOR CPMD. 
	!*) GRIDMIN		:	MINIMUM OF GRID.
	!*) GRIDMAX		:	MAXIMUM OF GRID.
	!*) GRIDDIF		:	BIN WIDTHF OF THE GRID.
	!*) NBIN		:	NUMBER OF BINS.
	!*) S			:	AN ORDERED N - TUPLE REPRESENTING POSITION OF A POINT IN THE GRID.
	!*) I			:	AN ORDERED N - TUPLE REPRESENTING INDEX OF A POINT IN THE GRID.
	!*) GRIDSIZE		:	TOTAL NUMBER OF GRID POINTS.
	!*) NMTD		:	NUMBER OF COLL. COORDINATES USED.
	!*) CENTER		:	MEAN OF GAUSSIAN.
	!*) HEIGHT		:	MAX. VALUE OF GAUSSIAN.
	!*) WIDTH		:	STANDARD DEVIATION OF GAUSSIAN.
	!*) IMIN		:	STARTING INDEX FOR INTEGRATION OF BIAS.
	!*) IMAX		:	ENDING INDEX FOR INTEGRATION BIAS.
	!*) BIAS		:	TOTAL BIAS.
	!*) C			:	REWEIGHTING FACTOR.
	!*) CNUMER		:	NUMERATOR OF C.
	!*) CDENOM		:	DENOMENATOR OF C.
	!*) CDUMP		:	FREQUENCY AT WHICH C IS PRINTED. MEASURED IN PER HILLS.
	!*) TEMP		:	TEMPERATURE.
	!*) WELLTEMP		:	FLAG FOR WELL TEMPERED METADYNAMICS.
	!*) WTDT		:	BIASING TEMPERATURE.
	!*) WTFAC		:	TEMPERING FACTOR.
	!*) GRIDPOINT		:	NEAREST GRID POINT TO COLLECTIVE VARIABLE IN TIME SERIES.
	!*) WEIGHT		:	STATISTICAL WEIGHT OF COLLECTIVE VARIABLE IN TIME SERIES.
	!*) INPUTFILE		:	INPUTFILE.
	!*) LOGFILE		:	LOGFILE.
	!*) OUTPUTFILE1		:	CONTAINS FREE ENERGY SURFACE FROM  ADDING GAUSSIANS.
	!*) OUTPUTFILE2		:	CONTAINS FREE ENERGY SURFACE FROM REWEIGHTING.
	!*) CFILE		:	CONTAINS REWEIGHTING FACTOR AS A FUNCTION OF HILLS.
	!*) ERRORFILE		:	CONTAINS LIST OF ERRORS.
	!*) ERRORSTAT		:	FLAG FOR CHECKING ERRORS.
	!*************************************************************! 

	!>>>>>ENVIR VARIABLES<<<<<!
	MODULE GLOBAL_VAR
	!GRID VARIABLES.
	INTEGER NDIM
	INTEGER*8 GRIDSIZE
	REAL*8, DIMENSION(:), ALLOCATABLE ::  GRIDMIN, GRIDMAX, GRIDDIF,S
	INTEGER, DIMENSION(:), ALLOCATABLE :: NBIN, I
	!METADYNAMICS VARIABLES.
	INTEGER NMTD, PACKAGE, CDUMP
	REAL*8 TEMP, WTDT, WTALPHA, CNUM, CDEN
	LOGICAL WELLTEMP
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: IMIN, IMAX
	REAL*8, DIMENSION(:), ALLOCATABLE :: WIDTH, HEIGHT, BIAS, C, CNUMER, CDENOM
	REAL*8, DIMENSION(:,:), ALLOCATABLE :: CENTER
	!TIME SERIES.
	REAL*8, DIMENSION(:), ALLOCATABLE :: GRIDPOINT
	REAL*8 WEIGHT
	!FILES.
	CHARACTER (LEN = 100) INPUTFILE, LOGFILE, OUTPUTFILE1, OUTPUTFILE2, ERRORFILE, CFILE
	LOGICAL ERRORSTAT
	!#IF DEFINED (_PARALLEL)                                                                    
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	!DEFINES VARIABLES FOR MPI 
	INTEGER MPIERROR, MPIBIND, MPIRANK, MPITASKS
	INTEGER, DIMENSION(:), ALLOCATABLE :: MPIBIN
	REAL*8, DIMENSION(:), ALLOCATABLE :: TBIAS, TCNUMER, TCDENOM
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	!#ENDIF
	END MODULE
	!>>>>>ENVIR VARIABLES<<<<<!
