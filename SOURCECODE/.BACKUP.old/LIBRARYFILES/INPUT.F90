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

	!*****************************************************************************************!
	!>>>>>>>> LIBRARIES <<<<<<<!
        !*****************************************************************************************!
	!INCLUDE 'VARIABLES.h'
	!INCLUDE 'UTIL.F90'
	!*****************************************************************************************!
	!>>>>>>>> LIBRARIES <<<<<<<!
        !*****************************************************************************************!

        !*****************************************************************************************!
	!>>>>>> SUBROUTINES <<<<<<<!
        !*****************************************************************************************!

	!INPUTS METADYNAMICS AND GRID PARAMETERS.
	SUBROUTINE INPUT
	USE GLOBAL_VAR
	IMPLICIT NONE
	INCLUDE 'mpif.h'
	CHARACTER (LEN=200) TOKEN, COMMENT, GRID
	INTEGER ITER1, ITER2
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!	
	!#IF DEFINED (_PARALLEL) 
	!READS FROM FILE ONLY IN PARENT PROCESSOR. 
	IF( MPIRANK .EQ. 0) THEN
	!#ENDIF
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
	!OPENS FILE.
		OPEN (UNIT = 1, FILE = INPUTFILE, STATUS ='OLD')
	!READS THE PACKAGE USED FOR METADYNAMICS.
		READ(1,*) TOKEN, PACKAGE
	!READS NUMBER OF DIMENSIONS AND CHECKS FOR ERROR.
		READ (1,*) TOKEN, NDIM
		IF (NDIM .LT. 1) THEN
			CALL PRINT_ERROR ('ERROR IN NUMBER OF DIMENSIONS.') 
			ERRORSTAT = .TRUE.
			RETURN
		ENDIF
		READ (1,*) TOKEN, NMTD
		IF (NMTD .LT. 1) THEN
			CALL PRINT_ERROR ('ERROR IN NUMBER OF METADYNAMICS STEPS.') 
			ERRORSTAT = .TRUE.
		RETURN
		ENDIF
	!READS TEMPERATURE.
		READ (1,*) TOKEN, TEMP
		IF (TEMP .LE. 0.0D0) THEN
			CALL PRINT_ERROR ('ERROR IN TEMPERATURE.')
			ERRORSTAT = .TRUE.
		ENDIF
	!CHANGES UNIT OF TEMPERATURE FROM K TO ATOMIC UNITS.
		TEMP = TEMP/3.1577464D5
	!ALLOCATES GRID VARIABLES..
		CALL ALLOCATE_ALL_GRID
		GRIDSIZE= 1
	!READS IF WELL TEMPERED.
		READ (1,*) TOKEN, WELLTEMP
	!READS THE BIASING TEMPERATURE IF WELL TEMPERED METADYNAMICS IS USED.
		IF (WELLTEMP .EQV. .TRUE.) THEN
			READ (1,*) WTDT
			IF (WTDT .LE. 0.0D0) THEN
				CALL PRINT_ERROR ('ERROR IN BIASING TEMPERATURE')
				ERRORSTAT = .TRUE.
			ENDIF
	!CHANGES UNIT OF BIASING TEMPERATURE FROM K TO ATOMIC UNITS.
			WTDT = WTDT/3.1577464D5
	!COMPUTES BIASING FACTOR.
			WTALPHA = (TEMP + WTDT) / WTDT
		ELSE 
			WTALPHA = 1.0D0
		ENDIF
	!INPUTS GRID PARAMETERS ABD CHECKS FOR ERROR. 
		READ (1,*) TOKEN
		IF (INDEX (TOKEN,'#') .EQV. 0) READ(1,*) TOKEN
		GRIDSIZE = 1
		DO ITER1 = 1,NDIM
			READ (1,*) GRIDMIN(ITER1), GRIDMAX(ITER1), GRIDDIF(ITER1) 
			NBIN(ITER1) =(GRIDMAX(ITER1) - GRIDMIN(ITER1)) / GRIDDIF(ITER1) + 1
			IF (NBIN(ITER1) .LE. 0) THEN
				CALL PRINT_ERROR ('ERROR IN GRID DETAILS.') 
				ERRORSTAT = .TRUE.
				RETURN
			ENDIF
	!CALCULATES SIZEOF BIAS.
			GRIDSIZE= GRIDSIZE*NBIN(ITER1) 
		ENDDO
	!READS IF C(T) IS TO BE PRINTED.
		READ (1,*) TOKEN, CDUMP
		IF (CDUMP .LE. 0) THEN
			CALL PRINT_ERROR ('ERROR IN PRINTING FREQUENCY OF REWEIGHTING FACTOR')
		ENDIF
	!ALLOCATES METADYNAMICS VARIABLES.
		CALL ALLOCATE_ALL_MTD
		CLOSE (1)
	ENDIF
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!  
	!#IF DEFINED(_PARALLEL) 
	!WAITS FOR OTHER PROCESSORS.. 
	CALL MPI_BARRIER (MPI_COMM_WORLD,MPIERROR)
	!SHARES ALL NON-ARRAY GRID VARIABLES WITH ALL PROCESSORS.
	!INTEGER TYPE VARIABLES.
	CALL MPI_BCAST (NDIM	, 1			, MPI_INTEGER		, 0, MPI_COMM_WORLD, MPIERROR)
	CALL MPI_BCAST (GRIDSIZE, 1			, MPI_INTEGER		, 0, MPI_COMM_WORLD, MPIERROR)
	IF (MPIERROR .NE. MPI_SUCCESS) THEN
                CALL PRINT_ERROR ('ERROR IN MPI BCAST OF ONE OF NON ARRAY GRID VARIABLES.')
                ERRORSTAT = .TRUE.
                CALL MPI_ABORT (MPI_COMM_WORLD, MPIBIND, MPIERROR)
		RETURN
        ENDIF
	!ALLOCATES ALL ARRAY GRID VARIABLES IN PROCESSORS OTHER THAN PARENT.
	IF (MPIRANK .NE. 0) THEN
	        CALL ALLOCATE_ALL_GRID
	ENDIF
	!SHARES ALL ARRAY GRID TYPE VARIABLES WITH ALL PROCESSORS.
	!INTEGER TYPE.
	CALL MPI_BCAST (NBIN	, NDIM			, MPI_INTEGER		, 0, MPI_COMM_WORLD, MPIERROR)
	!REAL TYPE.
	CALL MPI_BCAST (GRIDMIN	, NDIM			, MPI_DOUBLE_PRECISION	, 0, MPI_COMM_WORLD, MPIERROR)
	CALL MPI_BCAST (GRIDMAX	, NDIM			, MPI_DOUBLE_PRECISION	, 0, MPI_COMM_WORLD, MPIERROR)
	CALL MPI_BCAST (GRIDDIF	, NDIM			, MPI_DOUBLE_PRECISION	, 0, MPI_COMM_WORLD, MPIERROR)
	IF (MPIERROR .NE. MPI_SUCCESS) THEN
                CALL PRINT_ERROR ('ERROR IN MPI BCAST OF ONE OF ARRAY GRID VARIABLES.')
                ERRORSTAT = .TRUE.
                CALL MPI_ABORT (MPI_COMM_WORLD, MPIBIND, MPIERROR)
		RETURN
        ENDIF
	!#ENDIF
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!  
	!UPDATES LOG FILE IN PARENT PROCESSOR.
	IF (MPIRANK .EQ. 0) THEN
		OPEN (UNIT = 1, FILE = LOGFILE, ACCESS = 'APPEND', STATUS = 'OLD')
		WRITE(1,*) 
		WRITE(1,'(A,T40,A,I2)') 'NUMER OF DIMENSIONS', ':  ', NDIM
		GRID = ''
		DO ITER1 = 1,NDIM - 1
			WRITE (COMMENT,'(I5,A)') NBIN(ITER1), 'x'
			GRID =  TRIM (GRID) // COMMENT
		ENDDO
		WRITE (COMMENT,'(I5,A)') NBIN(ITER1), ' '
		WRITE(1,'(A,T40,A,A)') 'GRID', ':  ', ADJUSTL (TRIM (TRIM (GRID) // COMMENT))
		GRID = '['
		DO ITER1 = 1,NDIM - 1
			WRITE (COMMENT,'(F8.3,A)') GRIDMIN(ITER1), ', '
			GRID =  TRIM (GRID) // COMMENT
		ENDDO
		WRITE (COMMENT,'(F8.5,A)') GRIDMIN(ITER1), ' : '
		GRID =  TRIM (GRID) //COMMENT
		DO ITER1 = 1,NDIM - 1
        	        WRITE (COMMENT,'(F8.5,A)') GRIDMAX(ITER1), ', '
        	        GRID =  TRIM (GRID) // COMMENT
        	ENDDO
        	WRITE (COMMENT,'(F8.5,A)') GRIDMAX(ITER1), ']'
		WRITE(1,'(A,T40,A,A)') 'RANGE', ':  ', ADJUSTL (TRIM (TRIM (GRID) // COMMENT))
		CLOSE(1)
	ENDIF
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!
	
	!ALLOCATES VARIABLES.
	SUBROUTINE ALLOCATE_ALL_GRID
	USE GLOBAL_VAR
	IMPLICIT NONE
	!ALLOCATES GRID VARIABLES.
	ALLOCATE (GRIDMIN(NDIM), GRIDMAX(NDIM), GRIDDIF(NDIM), I(NDIM), NBIN(NDIM), S(NDIM))
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!	
        !#IF DEFINED (_PARALLEL)
	!ALLOCATES  
	ALLOCATE(MPIBIN(NDIM))	
        !#ENDIF
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
	RETURN
	END SUBROUTINE
 	!*****************************************************************************************!

	SUBROUTINE ALLOCATE_ALL_MTD
	USE GLOBAL_VAR
	!ALLOCATES METADYNAMICS VARIABLES.
	IMPLICIT NONE
	INTEGER ITER1, NGRID
	ALLOCATE (CENTER(NDIM,NMTD), WIDTH(NMTD), HEIGHT(NMTD), IMIN(NDIM,NMTD), IMAX(NDIM,NMTD))  
	ALLOCATE (BIAS(GRIDSIZE), C(NMTD / CDUMP), CNUMER(NMTD / CDUMP), CDENOM(NMTD / CDUMP))
	!INITIALIZES THE NUMERATOR AND DENOMENATOR OF THE REWEIGHTING FACTOR.
        !#IF DEFINED (_PARALLEL)
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!	
	ALLOCATE (TBIAS(GRIDSIZE), TCNUMER(NMTD / CDUMP), TCDENOM(NMTD / CDUMP + 1))
	NGRID = 1
	DO ITER1 = NDIM,2,-1
		NGRID = NGRID * NBIN(ITER1)
	ENDDO
	IF (MOD (NBIN(1),MPITASKS) .EQ. 0) THEN                                   	
		MPIBIN(1) = NBIN(1) / MPITASKS
		CNUM = DBLE (MPIBIN(1) * NGRID) * 1.0D0         	
        	CDEN = DBLE (MPIBIN(1) * NGRID) * 1.0D0
        	CNUMER(1) = CNUM
        	CDENOM(1) = CDEN 
        	C(1) = TEMP * DLOG(CNUMER(1)/CDENOM(1))
        ELSE
        	MPIBIN(1) = NBIN(1) / MPITASKS
        	IF (MPIRANK .NE. MPITASKS - 1) THEN
			CNUM = DBLE (MPIBIN(1) * NGRID) * 1.0D0        
			CDEN = DBLE (MPIBIN(1) * NGRID)* 1.0D0
                	CNUMER(1) = CNUM
                	CDENOM(1) = CDEN 
                	C(1) = TEMP * DLOG(CNUMER(1)/CDENOM(1))
        	ELSE
		CNUM = DBLE (NBIN(1) - MPIRANK * MPIBIN(1) * NGRID) * 1.0D0        
                CDEN = DBLE (NBIN(1) - MPIRANK * MPIBIN(1) * NGRID) * 1.0D0
                CNUMER(1) = CNUM
                CDENOM(1) = CDEN 
	        C(1) = TEMP * DLOG(CNUMER(1)/CDENOM(1))
        	ENDIF
        ENDIF
	!INITIALIZES THE REST OF THE ARRAY TYPE AVRIABLES.
	BIAS    = 0.0D0
        TBIAS	= 0.0D0
        TCNUMER	= 0.0D0
        TCDENOM	= 0.0D0
	!#ELSE
	!CNUM = DBLE (GRIDSIZE) * 1.0D0
	!CDEN = DBLE (GRIDSIZE) * 1.0D0
	!CNUMER(1) = CNUM
	!CDENOM(1) = CDEN 
	!C(1) = TEMP * DLOG(CNUMER(1)/CDENOM(1))
	!INITIALIZES THE REST OF THE ARRAY TYPE AVRIABLES.
	!BIAS = 0.0D0	
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
	!#ENDIF
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!
	
	!DEALLOCATES VARIABLES.
	SUBROUTINE DEALLOCATE_ALL
	USE GLOBAL_VAR
	IMPLICIT NONE
	DEALLOCATE(GRIDMIN, GRIDMAX, GRIDDIF, NBIN, S, I, CENTER, WIDTH, HEIGHT, BIAS, C) 
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!
	!>>>>>> SUBROUTINES <<<<<<<!
        !*****************************************************************************************!

