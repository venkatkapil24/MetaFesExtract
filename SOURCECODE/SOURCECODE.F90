        !*************************************************************!
        !*      SUM_HILLS-1.1                                        *!
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
	INCLUDE 'LIBRARYFILES/VARIABLES.F90'
	INCLUDE 'LIBRARYFILES/UTIL.F90'
	INCLUDE 'LIBRARYFILES/INPUT.F90'
	INCLUDE 'LIBRARYFILES/HILLS.F90'
	!*****************************************************************************************!
	!>>>>>>>> LIBRARIES <<<<<<<!
	!*****************************************************************************************!

	!*****************************************************************************************!
	!>>>>>>>> PROGRAM <<<<<<<<!
	!*****************************************************************************************!
	PROGRAM MAIN
	USE GLOBAL_VAR
	!USES MPIONLY IF PARALLEL TAG IS READ DURING EXECUTION.
	!#IF DEFINED (_PARALLEL) 
	INCLUDE 'mpif.h'
	!ENDIF
	!DEFINES VARIABLES FOR RECORDING TIME.

	!DEFINES VARIABLES FOR PROGRAM.
	REAL*8 STARTTIME, PARAMINPUTTIME, MTDINPUTTIME, INTEGTIME, PRINTTIME
	INTEGER ITER1
	CHARACTER (LEN = 200) COMMENT, COMMENTFORMAT, GRID
	!#IF DEFINED (_PARALLEL) 	
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	!INITIALIZES MPITHREAD ONLY IF PARALLEL TAG IS READ DURING EXECUTION.
	CALL MPI_INIT (MPIERROR) 
	!CHECKS IF MPIINITIALIZATION IS SUCCESSFUL.
	IF (MPIERROR .NE. MPI_SUCCESS) THEN
		CALL PRINT_ERROR ('ERROR IN MPI INITIALIZATION') 
		CALL MPI_ABORT (MPI_COMM_WORLD, MPIBIND, MPIERROR) 
		STOP
	ENDIF
	!GIVES INDICES TO CHILD PROCESSORS.
	CALL MPI_COMM_RANK( MPI_COMM_WORLD, MPIRANK, MPIERROR) 
	CALL MPI_COMM_SIZE (MPI_COMM_WORLD, MPITASKS, MPIERROR)
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	!#ENDIF
	INPUTFILE 	= 'INPUTFILES/INPUTFILE'
	OUTPUTFILE1	= 'OUTPUTFILES/FES-SUMHILLS'
	OUTPUTFILE2	= 'OUTPUTFILES/FES-REWEIGHT'
	ERRORFILE 	= 'OUTPUTFILES/ERRORFILE'
	LOGFILE 	= 'OUTPUTFILES/LOGFILE'
	CFILE		= 'OUTPUTFILES/CFILE'
	CVFILE          = 'HILLFILES/CVFILE'
	CALL CREATE_FILES
	CALL CPU_TIME (STARTTIME)
	!INPUTS GRID PARAMETERS AND CHECKS FOR ERROR.
	CALL INPUT
	IF (ERRORSTAT .EQV. .TRUE.) CALL MPI_ABORT (MPI_COMM_WORLD, MPIBIND, MPIERROR)
	CALL CPU_TIME (PARAMINPUTTIME)
	!INPUTS METADYNAMICS PARAMETERS AND CHECKS FOR ERRORS.
	CALL INPUT_HILLS
	IF (ERRORSTAT .EQV. .TRUE.) CALL MPI_ABORT (MPI_COMM_WORLD, MPIBIND, MPIERROR)
	CALL CPU_TIME (MTDINPUTTIME) 
	!CALCULATES TOTAL BIAS AND CHECKS FOR ERRORS.
	CALL CALC_HILLS
	IF (ERRORSTAT .EQV. .TRUE.) CALL MPI_ABORT (MPI_COMM_WORLD, MPIBIND, MPIERROR)
	CALL CPU_TIME (INTEGTIME)
	!PRINTS METADYNAMICS FREE ENERGY.
	CALL PRINT_HILLS
	CALL PRINT_PROB
	CALL PRINT_C
	CALL CPU_TIME (PRINTTIME)	
	!WRITES THE TIME TAKEN BY THE PROGRAM.
	!#IF DEFINED (_PARALLEL) 	
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	CALL MPI_BARRIER (MPI_COMM_WORLD, MPIERROR)
	IF( MPIRANK .EQ. 0) THEN
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	!#ENDIF	
		OPEN (UNIT = 1, FILE = LOGFILE, POSITION = "APPEND",ACTION = "WRITE") 
		WRITE (1,*) 
		WRITE (1,'(A)')'TIME TAKEN BY CODE DURING:'
		WRITE (1,'(A,F16.8)')'INPUT                     (s) =', (PARAMINPUTTIME - STARTTIME) 
		WRITE (1,'(A,F16.8)')'READING HILLS FILE        (s) =', (MTDINPUTTIME - PARAMINPUTTIME) 
		WRITE (1,'(A,F16.8)')'CALCULATING TOTAL BIAS    (s) =',(INTEGTIME - MTDINPUTTIME) 
		WRITE (1,'(A,F16.8)')'PRINTING THE FREE ENERGY  (s) =', (PRINTTIME - INTEGTIME) 
		WRITE (1,*) 
		WRITE (1,'(A,F16.8)')'TOTAL CPU TIME            (s) =', (PRINTTIME - STARTTIME) 
		WRITE (1,'(A)') "**********************************************"
		CLOSE (1)
	!#IF DEFINED (_PARALLEL) 	
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	ENDIF
	!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$! 
	!#ENDIF	
	!#IF DEFINED (_PARALLEL) 
	!DEALLOCATE (BIASTOTAL)
	CALL DEALLOCATE_ALL
	CALL MPI_FINALIZE (MPIERROR) 
	!#ENDIF
	END PROGRAM
	!*****************************************************************************************!
	!>>>>>>>> PROGRAM <<<<<<<<!
	!*****************************************************************************************!
