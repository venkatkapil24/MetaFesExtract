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
	CLOSe(1)
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
	DUMMY1 = GRIDSIZE
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
