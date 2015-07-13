        !*************************************************************!
        !*      Meta-FES-Extract-1.1                                 *!
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
	!>>>>>> SUBROUTINES <<<<<<<!
	!*****************************************************************************************!
	
	!CHOOSES SUBROUTINE FOR CALCULATINGMETADYNAMICS BIAS.
	SUBROUTINE CALC_HILLS
	USE GLOBAL_VAR
	IMPLICIT NONE
	!FOR 1D METADYNAMICS.
	IF (NDIM .EQ. 1) CALL CALC_HILLS_1D
	RETURN
	END SUBROUTINE
	!*****************************************************************************************!

	!CALCULATES 1-D METADYNAMICS BIAS.
	SUBROUTINE CALC_HILLS_1D
	USE GLOBAL_VAR
	IMPLICIT NONE
	INCLUDE 'mpif.h'
	INTEGER ITER1, ITER2, IBIAS, ITERCV, IOSTATUS
	REAL*8 DIFSQ
	INTEGER, DIMENSION(:), ALLOCATABLE :: IBEGIN, IEND
	ALLOCATE (IBEGIN(NDIM), IEND(NDIM))
	!OPENS THE CV TIME SERIES FILE.
	OPEN (UNIT = 1, FILE = CVFILE, STATUS = 'OLD')
	READ(1,*) ITERCV, (GRIDPOINT(ITER2), ITER2 = 1,NDIM)
	!INTEGRATES OVER HILLS.
	DO ITER1 = 1,NHILL + 1
	!REWEIGHS THE CV TIME SERIES.
		DO
			EXIT
			IF (ITERCV .GT. TIME(ITER1)) EXIT
			S(1:NDIM) = GRIDPOINT(1:NDIM)
			CALL GETI()
			IF ((ITER1 .GE. NHILLMIN) .OR. ITER1 .LE. NHILLMAX ) PROB(IBIAS ()) = PROB(IBIAS ()) + DEXP((TBIAS(IBIAS ()) - C(ITER1 / CDUMP)) / TEMP)
			READ(1,*) ITERCV, (GRIDPOINT(ITER2), ITER2=1,NDIM)
		ENDDO
		CALL MPI_BARRIER (MPI_COMM_WORLD, MPIERROR)
	!BAILS OUT.
		IF (ITER1 .GT. NHILL) EXIT
	!FINDS RANGE OF INTEGERATION.
		!#IF DEFINED (_SERIAL) 
			IBEGIN(1) = IMIN(1,ITER1)
			IEND(1) = IMAX(1,ITER1)
		!#ELSE
			IF (MOD (NBIN(1),MPITASKS) .EQ. 0) THEN
				MPIBIN(1) = NBIN(1) / MPITASKS
				IBEGIN(1) = MAX (MPIRANK * MPIBIN(1) + 1, IMIN(1,ITER1)) 
				IEND(1)   = MIN ((MPIRANK + 1) * MPIBIN(1), IMAX(1,ITER1))
			ELSE
				MPIBIN(1) = NBIN(1) / MPITASKS
				IF (MPIRANK .NE. MPITASKS - 1) THEN
					IBEGIN(1) = MAX (MPIRANK * MPIBIN(1) + 1, IMIN(1,ITER1))
					IEND(1)   = MIN ((MPIRANK + 1) * MPIBIN(1), IMAX(1,ITER1))
				ELSE
					
					IBEGIN(1) = MAX (MPIRANK * MPIBIN(1) + 1, IMIN(1,ITER1))
					IEND(1)   = MIN (NBIN(1), IMAX(1,ITER1))
				ENDIF
			ENDIF
		!#ENDIF	
	!INTEGRATES OVER GRIDS.
		I(1) = IBEGIN(1) - 1
		DO 
	!UPDATES GRID POINT.
			I(1) = I(1) + 1
	!BAIL OUT.
			IF (I(1) .GT. IEND(1)) EXIT
	!UPDATES THE BIAS AND REWEIGHTING FACTOR.
			ITER2= IBIAS ()
			CNUM = CNUM - DEXP (BIAS(ITER2) / TEMP)
			CDEN = CDEN - DEXP (BIAS(ITER2) / TEMP * (WTALPHA-1.0D0))
			BIAS(ITER2) = BIAS(ITER2) + HEIGHT(ITER1) * DEXP (-0.50D0*DIFSQ (ITER1)/WIDTH(ITER1)/WIDTH(ITER1))
			CNUM = CNUM + DEXP (BIAS(ITER2) / TEMP)
			CDEN = CDEN + DEXP (BIAS(ITER2) / TEMP * (WTALPHA-1.0D0))
		ENDDO
	!STORES NUMERATOR AND DENOMENATOR OF REWEIGHTING FACTOR EVERY CDUMP STEPS.            	
        	IF (MOD(ITER1,CDUMP) .EQ. 0) THEN
        		CNUMER(ITER1 / CDUMP + 1) = CNUM
        		CDENOM(ITER1 / CDUMP + 1) = CDEN
        	ENDIF
	!#IF DEFINED(_PARALLEL)
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
	!PATCHES VARIOUS SUB GRIDS.
	!CALCULATES TOTAL BIAS.
		CALL MPI_BARRIER (MPI_COMM_WORLD, MPIERROR)
		CALL MPI_ALLREDUCE (BIAS	, TBIAS	 , NGRID		, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIERROR)
		IF (MPIERROR .NE. MPI_SUCCESS) THEN
			CALL PRINT_ERROR ('ERROR IN CALCULATION OF TOTAL BIAS.')
        	        ERRORSTAT = .TRUE.
        	        CALL MPI_ABORT (MPI_COMM_WORLD, MPIBIND, MPIERROR)
        	        RETURN
		ENDIF
	!CALCULATES TOTAL REWEIGHTING FACTOR.
		CALL MPI_ALLREDUCE (CNUMER	, TCNUMER, NHILL / CDUMP + 1 	, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIERROR)
		CALL MPI_ALLREDUCE (CDENOM	, TCDENOM, NHILL / CDUMP + 1	, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIERROR)
		IF(MPIERROR .NE. MPI_SUCCESS) THEN
			CALL PRINT_ERROR ('ERROR IN CALCULATION OF TOTAL REWEIGHTING FACTOR.')
        		ERRORSTAT = .TRUE.
        		CALL MPI_ABORT (MPI_COMM_WORLD, MPIBIND, MPIERROR)
        		RETURN
		ENDIF
		C(1:NHILL / CDUMP + 1) = TEMP * DLOG (TCNUMER(1:NHILL / CDUMP + 1) / TCDENOM(1:NHILL / CDUMP + 1))
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
	!#ENDIF
	ENDDO
	CLOSE(1)
	RETURN
	END SUBROUTINE
	!*****************************************************************************************!

	!CHOOSES SUBORUTINE FOR PRINTING METADYNAMICS FREE ENEGRY SURFACE AS A FUNCTION OF COLLECTIVE VARIABES.
	SUBROUTINE PRINT_HILLS
	USE GLOBAL_VAR
	IMPLICIT NONE
	IF (NDIM .EQ. 1) CALL PRINT_HILLS_1D 
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!

	!PRINTS THE 1-D METADYNAMICS FREE ENERGY AS A FUNCTON OF COLLECTIVE VARIABLES.
	SUBROUTINE PRINT_HILLS_1D
	USE GLOBAL_VAR
	IMPLICIT NONE
	INTEGER ITER1, IBIAS
	REAL*8 SMALL1, SMALL2
	!OPENS OUTPUT FILE.
	IF( MPIRANK .EQ. 0) THEN
		OPEN (UNIT = 1, FILE = OUTPUTFILE1, STATUS = 'OLD')
		OPEN (UNIT = 2, FILE = OUTPUTFILE2, STATUS = 'OLD')
		SMALL1 = BIAS(1)
		SMALL2 = BIAS(1)
		DO ITER1 = 1,NBIN(1)
	!SELECTS INDEX OF POINT ON GRID.
			I(1) = ITER1
	!FINDS THE GRID POINT.
			CALL GETS
	!FINDS THE SMALLEST FREE ENERGY VALUE.
			SMALL1 = MIN (SMALL1, -WTALPHA * TBIAS(IBIAS ()))
			IF ( -TEMP * DLOG (PROB( IBIAS())) + 1 .NE. -TEMP * DLOG (PROB( IBIAS()))) THEN
			SMALL2 = MIN (SMALL2, -TEMP * DLOG (PROB( IBIAS())))
			ENDIF
		ENDDO
		DO ITER1 = 1,NBIN(1)
	!SELECTS INDEX OF POINT ON GRID.
			I(1) = ITER1
	!FINDS THE GRID POINT.
			CALL GETS
	!PRINTS THE VALUE OF BIAS.
			WRITE(1,*) S(1), -WTALPHA * TBIAS(IBIAS ()) -SMALL1
			IF ( -TEMP * DLOG (PROB( IBIAS())) + 1 .EQ. -TEMP * DLOG (PROB( IBIAS()))) THEN
				WRITE(2,*) S(1), -SMALL2 
			ELSE
				WRITE(2,*) S(1), -TEMP * DLOG (PROB( IBIAS())) -SMALL2 
			ENDIF
		ENDDO
		CLOSE(1)
		CLOSE(2)
	ENDIF
	RETURN
	END SUBROUTINE
	
	!PRINTS REWEIGHTING FACTOR AS A FUNCTION OF METADYNAMICS STEPS.
	SUBROUTINE PRINT_C
	USE GLOBAL_VAR
	IMPLICIT NONE
	INTEGER ITER1
	OPEN (UNIT = 1, FILE = CFILE, STATUS = 'UNKNOWN')
	
	IF (MPIRANK .EQ. 0) THEN
		DO ITER1 = 1,NHILL / CDUMP
			WRITE(1,*) ITER1, C(ITER1)
		ENDDO
	ENDIF
	CLOSE(1)
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!
	!>>>>>> SUBROUTINES <<<<<<<!
	!*****************************************************************************************!

