! $ Id$

      PROGRAM validationBTF
!
      USE definition
      IMPLICIT NONE
!     
!     --- Matrice au format CR (Compressed Rows) 
!          IPTR(I) : Pointeur sur le debut de la ligne I
!          JCN(IPTR(I): IPTR(I)-1): indices de col dans la ligne I
      INTEGER  :: N,NE
      INTEGER, DIMENSION (:), ALLOCATABLE :: IPTR, JCN
!
!     --- Declaration du graphe associe a la matrice
      type(cell_graphe),dimension(:),allocatable ::g
!
!     -- Donnees produites par l'algorithme de BTF
      INTEGER, DIMENSION (:), ALLOCATABLE :: PERM, PTDebCF
      INTEGER :: NbCF
!
!
      INTEGER  :: IERR
      CHARACTER :: DATA*80   ! nom du fichier
      INTEGER  :: INFILE,LP,TLEVEL
      LOGICAL  :: TRACE  
      DOUBLE PRECISION temps, trop_petit
      PARAMETER (trop_petit = 1.0D-6)

      INFILE=5   ! unité logique pour lire la matrice
      LP    =6   ! unité logique pour imprimer les messages
      TLEVEL=0   ! niveau de trace (0=Pas de trace)
!                  0   pas de trace
!                  >0  affichage de traces pendant le code BTF
!                  > 3 affichage complet
!                

      IERR=0
!
!     -------------------------------
!     -- LECTURE HEADER et ALLOCATION
!     -------------------------------
!     --- lecture de N (ordre) et de NE (nb de nonzeros)
      WRITE(LP,'(A)') ' ** Read Header of Input Matrix '
      CALL READ_N_NE(N,NE,INFILE,ierr)
      REWIND(INFILE)
      
      if(ierr.gt.0) goto 500
      WRITE(6,*) ' ** Read header: N, NE = ', N,NE
!
!     -- Allocation des données pour la matrice au format Compressed Row (CR)
      WRITE(6,*) ' ** Allocate matrix in Compressed Row (CR) format '
      IF (ALLOCATED(IPTR)) THEN 
         DEALLOCATE(IPTR)
      ENDIF
      IF (ALLOCATED(JCN)) THEN 
         DEALLOCATE(JCN)
      ENDIF
      ALLOCATE( IPTR ( N+1 ), stat=ierr )
      IF (ierr.gt.0) GOTO 500
      ALLOCATE( JCN ( NE ), stat=ierr )
      IF (ierr.gt.0) GOTO 500
!
!     -------------------------------------------
!     -- LECTURE de la matrice au format SIMPLE 
!     -- et PASSAGE au format Compressed Row (CR)
!     -------------------------------------------
      CALL READ_MATRIX_SF(N,NE,INFILE,IPTR,JCN,IERR)
      REWIND(INFILE)
      
      IF(IERR.GT.0)THEN
         WRITE(*,*)'Problem when reading the matrix'
         GOTO 500
      ENDIF
!
!     --------------------------------------
!     -- AlLOCATION ET CONTRUCTION DU GRAPHE 
!     --------------------------------------
      ALLOCATE(g(N), stat=ierr)
      IF (ierr.gt.0) GOTO 500
!     Construction du graphe
      CALL CRTO_graphe (g, N, IPTR, NE, JCN)
      
      if (TLEVEL.GT.3) call EDITION_GRAPHE(g,n)

!     ---------------------------------
!     -- AlLOCATION et calcul de la BTF
!     ---------------------------------
       ALLOCATE( PERM ( N ), stat=ierr )
       IF (ierr.gt.0) GOTO 500
       ALLOCATE( PTDebCF ( N ), stat=ierr )
       IF (ierr.gt.0) GOTO 500
!     
       NbCF       = -99999
       PTDebCF    = -88888
       PERM       = -77777
       IERR       = -12345
       TRACE = (TLEVEL.GT.0)  
          WRITE (LP,'(//A)') ' Appel de mon code BTF .... '

           call secdeb(temps)  ! Prise de debut du temps elapsed
           CALL TRIANG (N, g, TRACE, PERM, NbCF,        &
                              PTDebCF, IERR)
           call secfin(temps)  ! Prise de fin du temps elapsed

       write(6,*) '    '
       write(6,*) ' -------------------------------'
       write (6,*) ' IERREUR en retour de BTF     =', IERR
       write(6,*) ' N                             =', N
       write(6,*) ' Nb de composantes fortes      =', NbCF
       IF (temps.GT.trop_petit)   &
       write(6,*) ' Temps de calcul               =', temps
       write(6,*) ' -------------------------------'
       write(6,*) '    '
!
!     ----------------------------
!     -- Nettoyage du graphe
!     -- deallocation des données
!     ----------------------------
      CALL CLEAN_graphe (g, N, TRACE)
      DEALLOCATE(g)
      if (ALLOCATED (PERM)) DEALLOCATE (PERM)
      if (ALLOCATED (PTDebCF)) DEALLOCATE (PTDebCF)

      if (ALLOCATED (g)) DEALLOCATE (g)
      if (ALLOCATED (JCN)) DEALLOCATE (JCN)
      if (ALLOCATED (IPTR)) DEALLOCATE (IPTR)
      STOP
 500  CONTINUE
      WRITE(LP,'(A)') 'End of validation BTF'
      STOP
      END

      
      SUBROUTINE CRTO_graphe(g,dim,PTROW,NZ,INDCOL)
      USE definition
      IMPLICIT NONE
!     Creation du graphe g associe a une  matrice d'ordre dim
!     a partir d'une matrice intiale au format CR (Compressed Rows)
!     Le graphe g a ete alloue dans le programme
!     principal. 
      integer, intent(in)  :: dim,NZ
      integer, intent(in)  :: PTROW(DIM+1), INDCOL(NZ)
      type (cell_graphe),dimension(dim), intent(out) ::g
!     variables de travail
      type (cell_graphe),pointer::courant,deb
      integer :: retour
      integer  :: i,j
      print *, ' ** DEBUT creation_graphe **'
      do i=1,dim
!     Creation de la liste associee au noeud j
!     initialisation du degre
         g(i)%indice=0
         nullify(deb)
         do j=PTROW(i+1)-1,PTROW(i),-1
!     Insertion d'un nouvel element non nul sauf le terme diagonal
            if (INDCOL(j)/=i) then
               allocate(courant,stat=retour)
               if ( retour > 0 ) return
               courant%indice=INDCOL(j)
               courant%suivant=>deb
               deb=>courant
               g(i)%indice=g(i)%indice+1
            end if
         end do
!     Mise a jour du pointeur d'entree pour la liste associee au noeud i
         g(i)%suivant=>deb
      end do
      print *, ' ** FIN   creation_graphe **'
      return
      END SUBROUTINE CRTO_graphe
!
!

      SUBROUTINE READ_N_NE(N,NE,INFILE,ierr)
!
!   Objet:
!   -----
!     Read N and NZ from the File Descriptor INFILE
!
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  INFILE
      INTEGER, INTENT(OUT)::  N,NE,ierr
      ierr=0
      READ(INFILE,*,ERR=600) N
      READ(INFILE,*,ERR=600) NE
      RETURN
 600  CONTINUE
      WRITE(*,*)'File Not in SIMPLE format'
      ierr=1
      RETURN
      END SUBROUTINE READ_N_NE
      
      SUBROUTINE READ_MATRIX_SF(N,NE,INFILE,IPTR,JCN,IERR)
!     ---------------------------------------------
!
!   Objet:
!   -----
!     Lecture sur l'unité logique INFILE d'un fichier contenant
!     une matrice au format SIMPLE 
!     (Matrices du repertoire MAT.simple) a partir de son
!     descripteur (INFILE).
!     La matrice est alors transformée au format CR (Compressed Row).
!     
!   Parametres d'entree
!   --------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N,NE,INFILE
      INTEGER, DIMENSION(NE) , INTENT(OUT) :: JCN 
      INTEGER, DIMENSION(N+1), INTENT(OUT) :: IPTR
      INTEGER, INTENT(OUT) :: IERR
!
!      IN: Parametres fournis en entree 
!      --------------------------------
!      INFILE: INPUT de type entier, descripteur de fichier contenant la matrice
!           (open effectué dans le programme principal)
!        Rappel du format SIMPLE de matrice:
!        ligne 1     : N          ( l'ordre de la matrice N )
!        ligne 2     : NE         ( le nombre d'élements nonnuls )
!        Ligne 3     : I1   J1    ( élément non nul a(I1,J1))
!        Ligne 4     : I2   J2    ( élément non nul a(I2,J2))
!        .
!        .
!        .
!        Ligne NZ+2  : Inz  Jnz   ( élément non nul a(Inz,Jnz))
!
!   
!       N, NE : respectivement l'ordre et le nombres d'entrees non nulles de la
!               matrice
!
!      INOUT: Parametres fournis en entree et modifies en sortie
!      ---------------------------------------------------------
!      La matrice au format Compressed row:
!        IPTR: tableau de taille N+1 d'entiers tel que
!          IPTR(I) : position du début de la ligne I dans le tableau JCN. 
!        JCN : tableau de taille NE d'entiers tel que
!          JCN(IPTR(I): IPTR(I)-1): indices de colonne de la ligne I (non necessairement
!                                   ordonnés).
!     
!      INOUT: aucun
!      ------------
!
!  Variables locales
!   -----------------
      INTEGER, DIMENSION(:),ALLOCATABLE :: JCN_LOC,IRN_LOC
      INTEGER I,J,K, Nloc, NEloc
      LOGICAL TRACE
!
!     Verification N, NE associe au fichier
!     et allocation des tableaux pour lire 
!     la matrice locale au format SIMPLE.
!
      IERR=0
      READ(INFILE,*,ERR=600) Nloc
      READ(INFILE,*,ERR=600) NEloc
      IF ((N.NE.Nloc).OR.(NE.NE.NEloc)) GOTO 500

      IF(ALLOCATED(IRN_LOC))DEALLOCATE(IRN_LOC)
      ALLOCATE( IRN_LOC(NE), stat=ierr )
      IF (ierr.gt.0) GOTO 500
      IF(ALLOCATED(JCN_LOC))DEALLOCATE(JCN_LOC)
      ALLOCATE( JCN_LOC (NE), stat=ierr )
      IF (ierr.gt.0) GOTO 500

      DO I=1,NE
         READ(INFILE,*,ERR=600) IRN_LOC(I),JCN_LOC(I)
      ENDDO
      

      TRACE=.FALSE.

!     
!     Convertir la matrice du format SIMPLE en 
!     une matrice au format CR (Compressed Row)
!
      CALL SIMPLEtoCR ( N, NE, JCN_LOC, IRN_LOC,  &
                        IPTR, JCN, TRACE)

! 
!
      DEALLOCATE(JCN_LOC)
      DEALLOCATE(IRN_LOC)
      RETURN
 500  CONTINUE
      WRITE(*,*)'Allocation problem in READ_MATRIX_SF'
      IERR=1
      RETURN
 600  CONTINUE
      IERR=2
      RETURN
      END SUBROUTINE READ_MATRIX_SF

!     -----------------------------------------
      SUBROUTINE SIMPLEtoCR( N, NZ, JCN, IRN,    &
                             IPTR, JCNcr, TRACE)
!     -----------------------------------------
!
!   Objet:
!   -----
!     Convertir une matrice du format SIMPLE (voir repertoire MAT.simple)
!     au format Compressed Row (CR).
!
      IMPLICIT NONE
! -- input 
!     -- matrice stockée au format SIMPLE (IRN, JCN)
      INTEGER, INTENT(IN) ::  N, NZ
      INTEGER, INTENT(IN) :: JCN(NZ), IRN(NZ)
      LOGICAL, INTENT(IN) :: TRACE
! -- output 
!     -- matrice stockee au format Compressed Row (CR)
      INTEGER, INTENT(OUT):: IPTR(N+1), JCNcr(NZ)
      
!     -- working variables
      INTEGER I,J, NERR, MP, J1, J2, K, JPOS
      
!     ----------------------------------------------------------
! Check if there are zero columns.
      NERR =0
      MP = 6                    !  output unit for messages
!     
!     --- Set to zero the row pointers of A.
      IPTR = 0
!     
!     Scan the input, count the elements in each row and
!     ignore elements with indices out of range.
      NERR =0                   ! number of out of range indices
      DO K=1,NZ
         I =JCN(K)
         J =IRN(K)
         IF ( (I.LE.0 .OR. I.GT.N) .OR.    &
              (J.LE.0 .OR. J.GT.N) ) THEN
!     --- out of range row or column index
            NERR = NERR+1
            CYCLE               ! skip entry
         ENDIF
         IPTR(J) = IPTR(J) +1
      END DO

      IF (NERR.GT.0) THEN
         WRITE(MP,*) ' ********** ERROR detected in  SIMPLEtoCR'
         WRITE(MP,*) NERR, ' out of range indices in matrice '
      ENDIF
      
!     Prepare the row pointers in such a way that pointer IPTR(J)
!     points to the first position of column J+1.
      IPTR(1)=IPTR(1)+1
      DO K=2,N+1
         IPTR(K)= IPTR(K)+IPTR(K-1)
      ENDDO

!     Set the column pointers and store the matrix by column
      DO K=1,NZ
         I=JCN(K)
         J=IRN(K)
         JPOS=IPTR(J)
         JCNcr(JPOS-1)=I
         IPTR(J)=JPOS-1
      ENDDO

!     
      if (TRACE) then
         write(6,*) ' Input : Matrix au format SIMPLE :'
         write(6,*) '... IRN=', IRN
         write(6,*) '... JCN=', JCN
         
         write(6,*) ' Output : Matrix au format Compressed Row (CR) :'
         write(6,*) '... IPTR =', IPTR
         write(6,*) '... JCNcr=', JCNcr
      endif
!
      RETURN
      END SUBROUTINE SIMPLEtoCR
      

!*******************************************************
      SUBROUTINE secdeb(t)
!***************************************************
        DOUBLE PRECISION :: t,t1
        DOUBLE PRECISION :: ELAPSE8
        EXTERNAL ELAPSE8
!
        t=ELAPSE8(t1)
        return
        END SUBROUTINE SECDEB
!*******************************************************
       SUBROUTINE secfin(t)
!****************************************************
        DOUBLE PRECISION :: t
        DOUBLE PRECISION :: ELAPSE8, ZERO
        EXTERNAL ELAPSE8
!
        ZERO = 0.0
        t=ELAPSE8(ZERO)-t
        return
        END SUBROUTINE SECFIN


