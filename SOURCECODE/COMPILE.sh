#CREATES BACKUP OF SOURCECODE.
cp -r .BACKUP .BACKUP.old
cp -r LIBRARYFILES .BACKUP/
cp -r SOURCECODE.F90 .BACKUP/
#COMPILES
mpif90 -fp-stack-check -check all -traceback -gen-interfaces -warn interfaces SOURCECODE.F90 -o META-FES-EXTRACT.x
rm *.f90 *.mod 
