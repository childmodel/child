!Small utility to convert CHILD output to VTK format (vizualisation with softwares such as Paraview)
!usage: just type child2vtk in a terminal and you will be prompted for the base name of the output files
!Vincent Godard, CEREGE, Aix-Marseille University, France (godard@cerege.fr)
PROGRAM child2vtk
IMPLICIT NONE
! definitions
INTEGER, PARAMETER :: NTIME=1000
INTEGER :: I,J,K
REAL :: TIME
CHARACTER(40) :: pref,fileout
CHARACTER(3) :: meter
INTEGER :: centaine,dizaine,unite
CHARACTER(40) :: filenodes,filetri,filez,filetau, fileq, fileslp,filearea
REAL,dimension(:,:), ALLOCATABLE  :: TABNODES, TABZ, TABTAU, TABQ, TABSLP,TABAREA
INTEGER,dimension(:,:), ALLOCATABLE  :: TABTRI
INTEGER :: NN, NE, NN2

!lecture du prefixe
write(*,*)'Enter name of the run'
      read(*,'(a40)') pref

! ouverture fichier .nodes
filenodes=TRIM(TRIM(pref) // '.nodes')
write(*,*) 'Ouverture ', filenodes
OPEN(30,FILE=filenodes,STATUS='OLD')

! ouverture fichier .tri
filetri=TRIM(TRIM(pref) // '.tri')
write(*,*) 'Ouverture ', filetri
OPEN(40,FILE=filetri,STATUS='OLD')

! ouverture fichier .z
filez=TRIM(TRIM(pref) // '.z')
write(*,*) 'Ouverture ', filez
OPEN(50,FILE=filez,STATUS='OLD')

! ouverture fichier .tau
filetau=TRIM(TRIM(pref) // '.tau')
write(*,*) 'Ouverture ', filetau
OPEN(60,FILE=filetau,STATUS='OLD')

! ouverture fichier .q
fileq=TRIM(TRIM(pref) // '.q')
write(*,*) 'Ouverture ', fileq
OPEN(70,FILE=fileq,STATUS='OLD')

! ouverture fichier .slp
fileslp=TRIM(TRIM(pref) // '.slp')
write(*,*) 'Ouverture ', fileslp
OPEN(80,FILE=fileslp,STATUS='OLD')

! ouverture fichier .area
filearea=TRIM(TRIM(pref) // '.area')
write(*,*) 'Ouverture ', filearea
OPEN(90,FILE=filearea,STATUS='OLD')

DO K=1,NTIME


READ(30,*,end=1000) TIME
READ(30,*) NN
IF(K .EQ. 1) THEN 
   ALLOCATE(TABNODES(NN,2)) 
ENDIF 
DO I=1,NN
   READ(30,*) TABNODES(I,1),TABNODES(I,2)
ENDDO

READ(40,*) TIME
READ(40,*) NE
IF(K .EQ. 1) THEN 
   ALLOCATE(TABTRI(NE,3)) 
ENDIF 
DO I=1,NE
   READ(40,*) TABTRI(I,1),TABTRI(I,2),TABTRI(I,3)
ENDDO
TABTRI=TABTRI+1

READ(50,*) TIME
READ(50,*) NN
IF(K .EQ. 1) THEN 
   ALLOCATE(TABZ(NN,1)) 
ENDIF 
DO I=1,NN
   READ(50,*) TABZ(I,1)
ENDDO

READ(60,*) TIME
READ(60,*) NN
IF(K .EQ. 1) THEN 
   ALLOCATE(TABTAU(NN,1)) 
ENDIF 
DO I=1,NN
   READ(60,*) TABTAU(I,1)
ENDDO

READ(70,*) TIME
READ(70,*) NN
IF(K .EQ. 1) THEN 
   ALLOCATE(TABQ(NN,1)) 
ENDIF 
DO I=1,NN
   READ(70,*) TABQ(I,1)
ENDDO

READ(80,*) TIME
READ(80,*) NN
IF(K .EQ. 1) THEN 
   ALLOCATE(TABSLP(NN,1)) 
ENDIF 
DO I=1,NN
   READ(80,*) TABSLP(I,1)
ENDDO

READ(90,*) TIME
READ(90,*) NN2
IF(K .EQ. 1) THEN 
   ALLOCATE(TABAREA(NN,1)) 
ENDIF 
TABAREA=0.
DO I=1,NN2
   READ(90,*) TABAREA(I,1)
ENDDO


 centaine=K/100
 dizaine=(K-centaine*100)/10
 unite=(K-centaine*100-dizaine*10)
 meter(1:1)=CHAR(centaine+ICHAR('0'))
 meter(2:2)=CHAR(dizaine+ICHAR('0')) 
 meter(3:3)=CHAR(unite+ICHAR('0'))
 fileout=TRIM(TRIM(pref) // meter // '.vtk')

WRITE(*,*) K, fileout

OPEN (20,file=fileout,status='unknown')
write(20,'(a)')'# vtk DataFile Version 3.0'
write(20,'(a)')'CHILD'
write(20,'(a)')'ASCII'
write(20,'(a)')'DATASET UNSTRUCTURED_GRID'
write(20,'(a7,i10,a6)')'POINTS ',NN,' float'
DO I=1,NN
   write(20,'(3f15.5)')  TABNODES(I,1), TABNODES(I,2), TABZ(I,1)
ENDDO
write(20,'(A6, 2I10)') 'CELLS ',NE,(3+1)*NE
DO I=1,NE
   write(20,'(9I10)')3,TABTRI(I,1)-1,TABTRI(I,2)-1,TABTRI(I,3)-1
ENDDO
write(20,'(A11, I10)') 'CELL_TYPES ',NE
DO I=1,NE
   write(20,'(I2)')5
ENDDO
write(20,'(a11,i10)')'POINT_DATA ',NN
      write(20,'(a)')'SCALARS Altitude float 1'
      write(20,'(a)')'LOOKUP_TABLE default'
      DO I=1,NN
      write(20,*)TABZ(I,1)
      ENDDO
      write(20,'(a)')'SCALARS Shear_stress float 1'
      write(20,'(a)')'LOOKUP_TABLE default'
      DO I=1,NN
      write(20,*)TABTAU(I,1)
      ENDDO
      write(20,'(a)')'SCALARS Discharge float 1'
      write(20,'(a)')'LOOKUP_TABLE default'
      DO I=1,NN
      write(20,*)TABQ(I,1)
      ENDDO
      write(20,'(a)')'SCALARS Slope float 1'
      write(20,'(a)')'LOOKUP_TABLE default'
      DO I=1,NN
      write(20,*)TABSLP(I,1)
      ENDDO
      write(20,'(a)')'SCALARS Area float 1'
      write(20,'(a)')'LOOKUP_TABLE default'
      DO I=1,NN
      write(20,*)TABAREA(I,1)
      ENDDO

CLOSE(20)


ENDDO
1000 CONTINUE

write(*,*) 'Fermeture fichiers'
CLOSE(30)
CLOSE(40)
CLOSE(50)
CLOSE(60)
CLOSE(70)
CLOSE(80)
CLOSE(90)

RETURN
END