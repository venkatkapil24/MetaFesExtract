#CREATES INPUTFILE FOR VARYING NUMBER OF PROCESSORS IN DIFFERENT DIRECTORIES.

for I in `seq 1 1 32` ; do
#CREATES DISTICE FOLDER AND COPIES FILES MADATORY FOR EXECUTION.
 	if test $I -lt 10
	then
#ADDS A ZERO BEFORE INDEX IF SINGLE DIGIT.
		cp -r TEMPLATE NP.EQ.0$I
#ENTERS THE DIRECTORY.
		cd NP.EQ.0$I
	else
		cp -r TEMPLATE NP.EQ.$I
		cd NP.EQ.$I
	fi
#CHANGES THE NUMBER OF PROCESSORS IN SUBMIT FILE.
	sed -e s/numberofprocessors/$I/g SUBMIT > SUBMIT.sh
#REMOVES BUIFFER FILES.
	rm SUBMIT
	cd ..
done
