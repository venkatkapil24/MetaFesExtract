#WRITES THE WALLTIME AS A FUNCTION OF NUMBER OF PROCESSORS IN DIFFERENT DIRECTORIES.

for I in `seq 1 1 32` ; do
#CREATES DISTICE FOLDER AND COPIES FILES MADATORY FOR EXECUTION.
 	if test $I -lt 10
	then
#ADDS A ZERO BEFORE INDEX IF SINGLE DIGIT.
		tail -2 NP.EQ.0$I/LOGFILE | head -1 > FILE01
		awk '{print $5 }' FILE01 >> FILE
		
	else
		tail -2 NP.EQ.$I/LOGFILE  | head -1 > FILE01
		awk '{print $5 }' FILE01 >> FILE
	fi
	rm FILE01
done
