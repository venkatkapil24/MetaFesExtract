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
	!*****************************************************************************************!
	!>>>>>>>> LIBRARIES <<<<<<<!
        !*****************************************************************************************!

        !*****************************************************************************************!
	!>>>>>> SUBROUTINES <<<<<<<!
        !*****************************************************************************************!

	!CREATES EMPTY FILES.
	SUBROUTINE CREATE_FILES
	USE GLOBAL_VAR
	IMPLICIT NONE
	OPEN (UNIT = 1, FILE = OUTPUTFILE1, STATUS ='REPLACE') 
	CLOSE (1) 
	OPEN (UNIT = 1, FILE = OUTPUTFILE2, STATUS ='REPLACE') 
	CLOSE(1)
	OPEN (UNIT = 1, FILE = ERRORFILE, STATUS ='REPLACE') 
	CLOSE (1) 
	OPEN (UNIT = 1, FILE = LOGFILE, STATUS ='REPLACE') 
	!PRINTS THE LABEL.
	WRITE (1,'(A)') "**********************************************"
	WRITE (1,'(A)') "*                 SUM_HILLS-1.1              *"
	WRITE (1,'(A)') "*                                            *"
	WRITE (1,'(A)') "*      A SUPERFAST FREE ENEGRY               *"
	WRITE (1,'(A)') "*      RECONSTRUCTION CODE.                  *"
	WRITE (1,'(A)') "*                                            *"
	WRITE (1,'(A)') "*      DEVELOPED BY VENKAT KAPIL             *"
	WRITE (1,'(A)') "*      CONTACT : venkat@iitk.ac.in           *"
	WRITE (1,'(A)') "*                                            *"
	WRITE (1,'(A)') "*      COPYRIGHT (C) 2015 VENKAT KAPIL       *"
	WRITE (1,'(A)') "*      SOME RIGHTS RESERVED.                 *"
	WRITE (1,'(A)') "*                                            *"
	WRITE (1,'(A)') "**********************************************"
	WRITE (1,*) 
	CLOSE (1) 
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!
	
	!PRINTS ERROR MESSAGE.
	SUBROUTINE PRINT_ERROR (MESSAGE) 
	USE GLOBAL_VAR
	IMPLICIT NONE
	CHARACTER (LEN = *) MESSAGE
	OPEN (1, FILE = ERRORFILE, STATUS = "OLD", POSITION = "APPEND", ACTION = "WRITE") 
	WRITE (1,'(A)') MESSAGE
	CLOSE (1) 
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!
	
	!PRINTS MESSAGE IN LOG.
	SUBROUTINE PRINT_LOG (MESSAGE) 
	USE GLOBAL_VAR
	IMPLICIT NONE
	CHARACTER (LEN = *) MESSAGE
	OPEN (1, FILE = LOGFILE, STATUS = "OLD", POSITION = "APPEND", ACTION = "WRITE") 
	WRITE (1,'(A)') MESSAGE
	CLOSE (1) 
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!
	
	!FINDS N-TUPLEGRID POINT FROM INDEX.
	SUBROUTINE GETS
	USE GLOBAL_VAR
	IMPLICIT NONE
	S(1:NDIM) = GRIDMIN(1:NDIM) + DBLE(I(1:NDIM) - 1) * GRIDDIF(1:NDIM)
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!
	
	!FINDS N-TUPLEGRID INDEX FROM POINT.
	SUBROUTINE GETI
	USE GLOBAL_VAR
	IMPLICIT NONE
	I(1:NDIM) = (S(1:NDIM) - GRIDMIN(1:NDIM)) / GRIDDIF(1:NDIM) + 1
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!

	!ALLOCATES GRID VARIABLES.
        SUBROUTINE ALLOCATE_ALL_GRID
        USE GLOBAL_VAR
        IMPLICIT NONE
        ALLOCATE (GRIDMIN(NDIM), GRIDMAX(NDIM), GRIDDIF(NDIM), I(NDIM), NBIN(NDIM), S(NDIM))
        RETURN
        END SUBROUTINE
        !*****************************************************************************************!
                                                                                                                
        !ALLOCATES METADYNAMICS VARIABLES.
        SUBROUTINE ALLOCATE_ALL_MTD01
        USE GLOBAL_VAR
        IMPLICIT NONE
        ALLOCATE (IMETACOLVAR(NMETACOLVAR), CENTER(NMETACOLVAR,NHILL), WIDTH(NHILL), HEIGHT(NHILL), TIME(NHILL), IMIN(NMETACOLVAR,NHILL), IMAX(NMETACOLVAR,NHILL))  
	CENTER = 0.0D0
	BIAS = 0.0D0
        !#IF DEFINED (_PARALLEL)
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!	
        ALLOCATE (MPIBIN(NDIM))	
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
	!#ENDIF 
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!

	!ALLOCATES METADYNAMICS VARIABLES.	
	SUBROUTINE ALLOCATE_ALL_MTD02
	USE GLOBAL_VAR
	IMPLICIT NONE
	INTEGER ITER1
	!CALCULATES THE TOTAL NUMBER OF GRID IN METADYNAMICS COLLECTIVE VARAIBLES SPACE.
	NMETACOLVARGRID  = 1
	DO ITER1 = 1,NMETACOLVAR
		NMETACOLVARGRID = NMETACOLVARGRID * NBIN(IMETACOLVAR(ITER1))
	ENDDO
	ALLOCATE (BIAS(NMETACOLVARGRID))
	BIAS = 0.0D0
        !#IF DEFINED (_PARALLEL)
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!	
        ALLOCATE (TBIAS(NMETACOLVARGRID))
	TBIAS = 0.0D0
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
	!#ENDIF 
	RETURN
	END SUBROUTINE
        !*****************************************************************************************!

        !*****************************************************************************************!
	SUBROUTINE ALLOCATE_ALL_REWT
	USE GLOBAL_VAR
	IMPLICIT NONE
	INTEGER ITER1, NOTHERGRID
        ALLOCATE (C(NHILL / CDUMP + 1), CNUMER(NHILL / CDUMP + 1), CDENOM(NHILL / CDUMP + 1))
        ALLOCATE (PROB(NGRID), GRIDPOINT(NDIM))
	C    = 0.0D0
	!INITIALIZES THE NUMERATOR AND DENOMENATOR OF THE REWEIGHTING FACTOR.
        !#IF DEFINED (_PARALLEL)
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!	
        ALLOCATE (TPROB(NGRID), TCNUMER(NHILL / CDUMP + 1), TCDENOM(NHILL / CDUMP + 1))
        NOTHERGRID = 1
        DO ITER1 = NDIM,2,-1
        	NOTHERGRID = NOTHERGRID * NBIN(ITER1)
        ENDDO
        IF (MOD (NBIN(1),MPITASKS) .EQ. 0) THEN                                   	
        	MPIBIN(1) = NBIN(1) / MPITASKS
        	CNUM = DBLE (MPIBIN(1) * NOTHERGRID) * 1.0D0         	
        	CDEN = DBLE (MPIBIN(1) * NOTHERGRID) * 1.0D0
        	CNUMER(1) = CNUM
        	CDENOM(1) = CDEN 
        	C(1) = TEMP * DLOG(CNUMER(1)/CDENOM(1))
        ELSE
        	MPIBIN(1) = NBIN(1) / MPITASKS
        	IF (MPIRANK .NE. MPITASKS - 1) THEN
        		CNUM = DBLE (MPIBIN(1) * NOTHERGRID) * 1.0D0        
        		CDEN = DBLE (MPIBIN(1) * NOTHERGRID)* 1.0D0
                	CNUMER(1) = CNUM
                	CDENOM(1) = CDEN 
                	C(1) = TEMP * DLOG(CNUMER(1)/CDENOM(1))
        	ELSE
        	CNUM = DBLE (NBIN(1) - MPIRANK * MPIBIN(1) * NOTHERGRID) * 1.0D0        
                CDEN = DBLE (NBIN(1) - MPIRANK * MPIBIN(1) * NOTHERGRID) * 1.0D0
                CNUMER(1) = CNUM
                CDENOM(1) = CDEN 
                C(1) = TEMP * DLOG(CNUMER(1)/CDENOM(1))
        	ENDIF
        ENDIF
        !INITIALIZES THE REST OF THE ARRAY TYPE AVRIABLES.
        TPROB   = 0.0D0
        TCNUMER	= 0.0D0
        TCDENOM	= 0.0D0
        !#ELSE
        !CNUM = DBLE (NGRID) * 1.0D0
        !CDEN = DBLE (NGRID) * 1.0D0
        !CNUMER(1) = CNUM
        !CDENOM(1) = CDEN 
        !C(1) = TEMP * DLOG(CNUMER(1)/CDENOM(1))
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

	!*****************************************************************************************!
        !>>>>>>> FUNCTIONS <<<<<<<<!
        !*****************************************************************************************!

	!CALCULATES DISTANCE BETWEEN THE CENTRE OF GAUSSIAN AND GRID POINT.
	REAL*8 FUNCTION DIFSQ(ITER1)
	USE GLOBAL_VAR
	IMPLICIT NONE
	INTEGER ITER1, ITER2
	DIFSQ = 0.0D0
	!UPDATES GRID COORDINATES.
	CALL GETS
	DO ITER2 = 1,NDIM
		DIFSQ = DIFSQ + (CENTER(ITER2,ITER1) - S(ITER2)) * (CENTER(ITER2,ITER1) - S(ITER2))
	ENDDO
	RETURN
	END FUNCTION
        !*****************************************************************************************!
	
	!FINDS THE LINEAR INDEX OF BIAS ARRAY.
	INTEGER*8 FUNCTION IBIAS
	USE GLOBAL_VAR
	IMPLICIT NONE
	INTEGER ITER1
	INTEGER*8 DUMMY1
	DUMMY1 = NGRID
	IBIAS = I(1)
	DO ITER1 = NDIM,2, - 1
		DUMMY1 = DUMMY1/NBIN(ITER1)
		IBIAS = IBIAS + DUMMY1 * (I(ITER1 - 1))
	ENDDO
	RETURN
	END FUNCTION
        !*****************************************************************************************!
        !>>>>>>> FUNCTIONS <<<<<<<<!
        !*****************************************************************************************!
