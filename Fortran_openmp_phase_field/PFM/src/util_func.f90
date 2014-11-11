Module UTIL_FUNC
! THIS MODULE CONTINS UTILITY FUNCTIONS AND SUBROUTINES THAT USED BY OTHER MODULES
! T_PROFILE()		THE THERMAL PATH, CONSISTING OF A SERIES OF LINEAR SEGMENTS, BEING CALLED IN EVERY TIME STEP TO UPDATE TEMPERATURE
! MYNGB(A,B,N)		RETURN THE INDEX OF INDEX A'S N-TH NEIGHBOUR IN ONE DIMENSION WITH B POINTS
! PHASE_INDEX(N)	RETURN THE PHASE INDEX OF PHASE FIELD N
! RANDOM_INT(MINV,MAXV) RETURN A RANDOM INTEGER BETWEEN MINV AND MAXV
! DIS(I,J,K,X,Y,Z)	RETURN THE DISTANCE BETWEEN POINT(I,J,K) AND POINT(X,Y,Z)
! DIS1(I,J,K,X,Y,Z)	RETURN THE DISTANCE BETWEEN POINT(I,J,K) AND POINT(X,Y,Z) CONSIDERING THE THICKNESS REDUCTION
! MATCH(NN,M)		RETURN THE INDEX 'I' OF VECTOR 'M'THAT M[I]==NN, OTHERWISE RETURN 0
! GET_NGB1(I,J,K,NGB)	GET THE INDICES OF NEIGBOURS OF POINT(I,J,K) IN THE MATRIX NGB)
! SSORT(X,Y,N)		QUICK SORT X AND Y BASED ON X
! SSORT(X,Y,Z,N)	QUICK SORT X,Y AND Z BASED ON X

      Use COMMON_VARIABLE
      Implicit None
Contains
      Subroutine T_PROFILE ()
!
         Real Lowertemp, Uppertemp, Lowertime, Uppertime
         Integer :: I = 0
!
         If (Time >= T_path(1, I+1)) I = I + 1
!
         Lowertemp = T_path (2, I)
         Uppertemp = T_path (2, I+1)
         Lowertime = T_path (1, I)
         Uppertime = T_path (1, I+1)
!
         Temperature = (Uppertemp-Lowertemp) / (Uppertime-Lowertime) * &
        & (Time-Lowertime) + Lowertemp
!
      End Subroutine
      Integer Function MYNGB (Coord, Lcoord, N)
!
         Integer Coord, Lcoord, N
         Integer Nn
!        Myngb=Modulo(Coord-1+Mod(N,Lcoord),Lcoord)+1
!
         Nn = Modulo (N, Lcoord)
         Myngb = Coord + Nn
         If (Myngb > Lcoord) Then
            Myngb = Myngb - Lcoord
         Else If (Myngb < 1) Then
            Myngb = Myngb + Lcoord
         End If
!
      End Function
!
!
      Elemental Integer Function PHASE_INDEX (I)
         Integer, Intent (In) :: I
         If (I <= 0) Then
            Phase_index = - 10000
         Else If (I > Npp1) Then !!! Austenite
            Phase_index = 2
         Else If (I > Stoich_n) Then !!! Ferrite
            Phase_index = 1
         Else !If(I<Npp2)Then  !!! Pearlite Or Cementite
            Phase_index = 3
         End If
      End Function
!
      Function RANDOM_INT (Minv, Maxv)
         Real Xx
         Integer Random_int, Minv, Maxv
         Call Random_number (Xx)
         Random_int = Nint (Minv+Xx*(Maxv-Minv+1)-0.5)
!
      End Function
!
      Pure Real Function DIS (I, J, K, X, Y, Z)! Ix,Iy,Iz,Width Are Used.
!
         Integer, Intent (In) :: I, J, K
         Real, Intent (In) :: X, Y, Z
         Real Xx, Yy, Zz
!
         Xx = Min (((I-1)*Width-X)**2, &
        & (Ix*Width-Abs((I-1)*Width-X))**2)
         Yy = Min (((J-1)*Width-Y)**2, &
        & (Iy*Width-Abs((J-1)*Width-Y))**2)
         Zz = Min (((K-1)*Width-Z)**2, &
        & (Iz*Width-Abs((K-1)*Width-Z))**2)
         Dis = Sqrt (Xx+Yy+Zz)
!
      End Function
!
      Real Function DIS1 (I, J, K, X, Y, Z, Rt)! Ix,Iy,Iz,Width Are Used.
!
         Integer I, J, K
         Real X, Y, Z
         Real Xx, Yy, Zz
         Real Rt
!
         Xx = Min (((I-1)*Width-X)**2, &
        & (Ix*Width-Abs((I-1)*Width-X))**2)
         Yy = Min (((J-1)*Width-Y)**2, &
        & (Iy*Width-Abs((J-1)*Width-Y))**2)
         Zz = Min (((K-1)*Width-Z)**2, &
        & (Iz*Width-Abs((K-1)*Width-Z))**2)
         Dis1 = Sqrt (Rt**2*Xx+Yy/Rt/Rt+Zz)
!
      End Function
!
      Function MATCH (Nn, Matrix)
!
         Integer Matrix (:), Nn, Match, I
         Match = 0
         Do I = 1, Size (Matrix)
            If (Matrix(I) == Nn) Then
               Match = I
               Return
            End If
         End Do
!
      End Function
!
      Subroutine GET_NGB1 (I, J, K, Ngb)
         Integer I, J, K, Ngb (3, 0:6), Ii, Jj, Kk
!Dir$ Vector
         Do Jj = 0, Nngb
!
            Ngb (1:3, Jj) = (/ I, J, K /)
         End Do
!Dir$ Ivdep
         Do Ii = 1, Dms
            Ngb (Ii, (2*Ii-1) :2*Ii) = Ngb (Ii, (2*Ii-1) :2*Ii) + (/ - &
           & 1, 1 /)
         End Do
         If (I == 1) Ngb (1, 1) = Ix
         If (J == 1) Ngb (2, 3) = Iy
         If (K == 1) Ngb (3, 5) = Iz
         If (I == Ix) Ngb (1, 2) = 1
         If (J == Iy) Ngb (2, 4) = 1
         If (K == Iz) Ngb (3, 6) = 1
!
      End Subroutine
! SUBROUTINE GET_NGB2(I,J,K,NGB)
!  integer i,j,k,ngb(3,0:nngb),ii
!
!  ngb(:,0)=(/i,j,k/)
!  ngb(:,1)=(/mod(i+ix-2,ix)+1,j,k/)
!  ngb(:,2)=(/mod(i+ix  ,ix)+1,j,k/)
!  ngb(:,3)=(/i,mod(j+iy-2,iy)+1,k/)
!  ngb(:,4)=(/i,mod(j+iy  ,iy)+1,k/)
!  ngb(:,5)=(/i,j,mod(k+iz-2,iz)+1/)
!  ngb(:,6)=(/i,j,mod(k+iz  ,iz)+1/)
!
! END SUBROUTINE
!!!!***********************************************
      Subroutine SSORT (x, y, N)
!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!!     .. Scalar Arguments ..
         Integer N
!     .. Array Arguments ..
         Real x (:)
         Integer y (:)!,Z(:)
!     .. Local Scalars ..
         Real R, T, TT, TTY, TY !,TZ,TTZ
         Integer i, IJ, j, k, L, M
!     .. Local Arrays ..
         Integer IL (21), IU (21)
!
         x (1:N) = - x (1:N)
!
!     Sort X and carry Y along
!
         M = 1
         i = 1
         j = N
         R = 0.375E0
!
110      If (i .Eq. j) Go To 150
         If (R .Le. 0.5898437E0) Then
            R = R + 3.90625E-2
         Else
            R = R - 0.21875E0
         End If
!
120      k = i
!
!     Select a central element of the array and save it in location T
!
         IJ = i + Int ((j-i)*R)
         T = x (IJ)
         TY = y (IJ)
     ! TZ=Z(IJ)
!
!     If first element of array is greater than T, interchange with T
!
         If (x(i) .Gt. T) Then
            x (IJ) = x (i)
            x (i) = T
            T = x (IJ)
            y (IJ) = y (i)
            y (i) = TY
            TY = y (IJ)
        ! Z(IJ)=Z(I)
        ! Z(I)=TZ
        ! TZ=Z(IJ)
         End If
         L = j
!
!     If last element of array is less than T, interchange with T
!
         If (x(j) .Lt. T) Then
            x (IJ) = x (j)
            x (j) = T
            T = x (IJ)
            y (IJ) = y (j)
            y (j) = TY
            TY = y (IJ)
        ! Z(IJ)=Z(J)
        ! Z(J)=TZ
        ! TZ=Z(IJ)
!
!        If first element of array is greater than T, interchange with T
!
            If (x(i) .Gt. T) Then
               x (IJ) = x (i)
               x (i) = T
               T = x (IJ)
               y (IJ) = y (i)
               y (i) = TY
               TY = y (IJ)
           ! Z(IJ)=Z(I)
           ! Z(I)=TZ
            !TZ=Z(IJ)
            End If
         End If
!
!     Find an element in the second half of the array which is smaller
!     than T
!
130      L = L - 1
         If (x(L) .Gt. T) Go To 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
140      k = k + 1
         If (x(k) .Lt. T) Go To 140
!
!     Interchange these elements
!
         If (k .Le. L) Then
            TT = x (L)
            x (L) = x (k)
            x (k) = TT
            TTY = y (L)
            y (L) = y (k)
            y (k) = TTY
        ! TTZ=Z(L)
        ! Z(L)=Z(K)
        ! Z(K)=TTZ
            Go To 130
         End If
!
!     Save upper and lower subscripts of the array yet to be sorted
!
         If (L-i .Gt. j-k) Then
            IL (M) = i
            IU (M) = L
            i = k
            M = M + 1
         Else
            IL (M) = k
            IU (M) = j
            j = L
            M = M + 1
         End If
         Go To 160
!
!     Begin again on another portion of the unsorted array
!
150      M = M - 1
         If (M .Eq. 0) Go To 190
         i = IL (M)
         j = IU (M)
!
160      If (j-i .Ge. 1) Go To 120
         If (i .Eq. 1) Go To 110
         i = i - 1
!
170      i = i + 1
         If (i .Eq. j) Go To 150
         T = x (i+1)
         TY = y (i+1)
     ! TZ=Z(I+1)
         If (x(i) .Le. T) Go To 170
         k = i
!
180      x (k+1) = x (k)
         y (k+1) = y (k)
     ! Z(K+1)=Z(K)
         k = k - 1
         If (T .Lt. x(k)) Go To 180
         x (k+1) = T
         y (k+1) = TY
     ! Z(K+1)=TZ
         Go To 170
!
!     !lean up
!
190      Do 200 i = 1, N
            x (i) = - x (i)
200      Continue
!
         Return
!
      End Subroutine
!!!!***********************************************
      Subroutine SSORT2 (x, y, z, N)
!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and optionally make the same interchanges in
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!!     .. Scalar Arguments ..
         Integer N
!     .. Array Arguments ..
         Real x (:), z (:)
         Integer y (:)!,Z(:)
!     .. Local Scalars ..
         Real R, T, TT, TTY, TY, TZ, TTZ
         Integer i, IJ, j, k, L, M
!     .. Local Arrays ..
         Integer IL (21), IU (21)
!
         x (1:N) = - x (1:N)
!
!     Sort X and carry Y along
!
         M = 1
         i = 1
         j = N
         R = 0.375E0
!
110      If (i .Eq. j) Go To 150
         If (R .Le. 0.5898437E0) Then
            R = R + 3.90625E-2
         Else
            R = R - 0.21875E0
         End If
!
120      k = i
!
!     Select a central element of the array and save it in location T
!
         IJ = i + Int ((j-i)*R)
         T = x (IJ)
         TY = y (IJ)
         TZ = z (IJ)
!
!     If first element of array is greater than T, interchange with T
!
         If (x(i) .Gt. T) Then
            x (IJ) = x (i)
            x (i) = T
            T = x (IJ)
            y (IJ) = y (i)
            y (i) = TY
            TY = y (IJ)
            z (IJ) = z (i)
            z (i) = TZ
            TZ = z (IJ)
         End If
         L = j
!
!     If last element of array is less than T, interchange with T
!
         If (x(j) .Lt. T) Then
            x (IJ) = x (j)
            x (j) = T
            T = x (IJ)
            y (IJ) = y (j)
            y (j) = TY
            TY = y (IJ)
            z (IJ) = z (j)
            z (j) = TZ
            TZ = z (IJ)
!
!        If first element of array is greater than T, interchange with T
!
            If (x(i) .Gt. T) Then
               x (IJ) = x (i)
               x (i) = T
               T = x (IJ)
               y (IJ) = y (i)
               y (i) = TY
               TY = y (IJ)
               z (IJ) = z (i)
               z (i) = TZ
               TZ = z (IJ)
            End If
         End If
!
!     Find an element in the second half of the array which is smaller
!     than T
!
130      L = L - 1
         If (x(L) .Gt. T) Go To 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
140      k = k + 1
         If (x(k) .Lt. T) Go To 140
!
!     Interchange these elements
!
         If (k .Le. L) Then
            TT = x (L)
            x (L) = x (k)
            x (k) = TT
            TTY = y (L)
            y (L) = y (k)
            y (k) = TTY
            TTZ = z (L)
            z (L) = z (k)
            z (k) = TTZ
            Go To 130
         End If
!
!     Save upper and lower subscripts of the array yet to be sorted
!
         If (L-i .Gt. j-k) Then
            IL (M) = i
            IU (M) = L
            i = k
            M = M + 1
         Else
            IL (M) = k
            IU (M) = j
            j = L
            M = M + 1
         End If
         Go To 160
!
!     Begin again on another portion of the unsorted array
!
150      M = M - 1
         If (M .Eq. 0) Go To 190
         i = IL (M)
         j = IU (M)
!
160      If (j-i .Ge. 1) Go To 120
         If (i .Eq. 1) Go To 110
         i = i - 1
!
170      i = i + 1
         If (i .Eq. j) Go To 150
         T = x (i+1)
         TY = y (i+1)
         TZ = z (i+1)
         If (x(i) .Le. T) Go To 170
         k = i
!
180      x (k+1) = x (k)
         y (k+1) = y (k)
         z (k+1) = z (k)
         k = k - 1
         If (T .Lt. x(k)) Go To 180
         x (k+1) = T
         y (k+1) = TY
         z (k+1) = TZ
         Go To 170
!
!     !lean up
!
190      Do 200 i = 1, N
            x (i) = - x (i)
200      Continue
!
         Return
!
      End Subroutine
      Function poly4 (x, a, b, c, d)
!
         Real a, b, c, d, x, poly4
!
         poly4 = a * x + b * x * x + c * x * x * x + d * x * x * x * x
      End Function
      Function diff_poly4 (x, a, b, c, d)
!
         Real x, a, b, c, d, diff_poly4
!
         diff_poly4 = a + 2 * b * x + 3 * c * x * x + 4 * d * x * x * x
!
      End Function
!
      Subroutine phase_list
!
         Integer austenite_vol, ferrite_vol, carbide_vol
!
         austenite_vol = count (ind_ph1(1, :, :, :) > npp1)
         carbide_vol = count (ind_ph1(1, :, :, :) <= stoich_n)
         ferrite_vol = ix * iy * iz - austenite_vol - carbide_vol
         Write (12,*) 'Phases:'
!
         If (ferrite_vol > 0) write (12,*) '    Ferrite'
         If (austenite_vol > 0) write (12,*) '    Austenite'
         If (carbide_vol > 0) write (12,*) '    Carbide'
!
         Write (12,*) '====================='
      End Subroutine
!
      Function avg_size (phase)
         Integer phase
         Integer i, j, min_num, max_num, total_num
         Real total_size, single_size, avg_size
!
!
         If (phase == 2) Then
            min_num = npp1 + 1
            max_num = tpp
         Else If (phase == 1) Then
            min_num = stoich_n + 1
            max_num = npp0
         Else If (phase == 3) Then
            min_num = 1
            max_num = stoich_n
         Else If (phase == 4) Then
            min_num = npp0 + 1
            max_num = npp1
         Else
            avg_size = - 1
            Return
         End If
!
         total_num = 0
!
         Do i = min_num, max_num
            single_size = count (ind_ph1(1, :, :, :) == i)
!
            If (dms == 1) Then
               If (single_size < 2*i_intf) Cycle
               total_size = total_size + single_size * width
            Else If (dms == 2) Then
               If (single_size < pi*i_intf**2) Cycle
               total_size = total_size + Sqrt (4.*single_size/pi) * &
              & width
            Else If (dms == 3) Then
               If (single_size < 4*i_intf**3) Cycle
               total_size = total_size + (6.*single_size/pi) ** (1./3.) &
              & * width
            End If
            total_num = total_num + 1
         End Do
         If (total_num == 0) Then
            avg_size = 0
         Else
            avg_size = total_size / total_num
         End If
!
      End Function
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!DATA OUTPUT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Subroutine WRITE_GRAIN_DATA (T, STAT)
!
         Real T, tmp, vf_stoich, vf_deff, vf_rexf, vf_aust
         Integer i, j, k, grain_vol
         Character (Len=10) :: String_time
         Logical STAT
!
!
         Write (String_time, "(F9.1)") T
!//VOLUME FRACTION
         vf_stoich = count (ind_ph1(1, :, :, :) <= stoich_n) * 1. / &
        & (ix*iy*iz)!!!STOICH FRACTION
         vf_deff = count (ind_ph1(1, :, :, :) <= npp0) * 1. / &
        & (ix*iy*iz) - vf_stoich
         vf_rexf = count (ind_ph1(1, :, :, :) <= npp1) * 1. / &
        & (ix*iy*iz) - vf_deff - vf_stoich !!!REX ALPHA FRACTION
         vf_aust = 1 - vf_stoich - vf_deff - vf_rexf !!!AUSTENITE FRACTION
!
         If (all(phase_index(ind_ph1) /= 3)) Then
            mobl0 (:, 3) = 0.
            mobl0 (3, :) = 0.
            Write (*,*) "pearlite all dissolved"
!
         End If
!   if(all(phase_index(ind_ph1)/=1))then
!        diff0(1)=0.
!        mobl0(1,:)=0.
!        mobl0(:,1)=0.
!  endif
!
         Write (10, '(E10.3,F6.1,F6.3,F6.3,F6.3,I5,F5.1,F5.1,F6.2)') T, &
        & temperature - 273.15, vf_rexf, vf_aust, vf_stoich, tpp - &
        & npp1, avg_size (4), avg_size (2), sum (val_rate, val_rate > &
        & 0.) / count (val_rate > 0.)
!!!!screen output
         Write (*, "('Time Step:',E10.2,'s',E10.2,'s')") dt_diff, del_t
         Write (12, "('Time Step:',E10.2,'s',E10.2,'s')") dt_diff, &
        & del_t
!   write(*,*)"        PARTITION_CONEFFICIENT  "
!   DO I=1,N_PHASE
!    ! DO J=1,N_PHASE
!     WRITE(*,'(3E9.2)')(PARTITION_K(I,J),J=1,N_PHASE)
!   ENDDO
         Write (*,*) ' Time ', 'T(C) ', '  AVE CONC'
         Write (12,*) ' Time ', 'T(C) ', 'AUS_SIZE  AVE CONC'
!
         Write (*, "(E10.3,F6.1,F6.3)") T, temperature - 273., sum &
        & (val_con1*val_phase) / (ix*iy*iz)
         Write (12, "(E10.3,F6.1,F6.3)") T, temperature - 273., sum &
        & (val_con1*val_phase) / (ix*iy*iz)
         Write (*,*)
!
         If (STAT) Then
!
!//STRUCTURE
            command = foldername
            command ((len_trim(foldername)+1) :) = '/micr_time' // trim &
           & (adjustl(String_time)) // '.vti'
!
            Call structure_out (command)
!
         End If
!
!
!
         If (vf_aust > 0.999) Then
!
            command ((len_trim(foldername)+1) :) = '/grain_size' // &
           & trim (adjustl(String_time)) // 's.txt'
            Open (11, File=command)
            If (iz == 1) Then
               Write (11,*) '  Diameter (um)  '
!
               Do i = npp1 + 1, tpp
                  grain_vol = count (ind_ph1(1, :, :, :) == i)
                  Write (11, '((f15.6))') (grain_vol*dms*2.0/pi) ** &
                 & (1./dms) * width
               End Do
            Else If (iz > 1) Then
               Write (11,*) '  Diameter (um)  '
               k = 0
               Do while (k <  2*alpha_nuc_N)
                  Call random_number (tmp)
                  Do i = npp0 + 1, npp1
                     grain_vol = count (ind_ph1(1, :, :, Nint(tmp*iz)) &
                    & == i)
                     If (grain_vol > 0) Then
                        Write (11, '(2(f15.6))') (grain_vol*dms*2.0/pi) &
                       & ** (1./dms) * width
                        k = k + 1
                     End If
                     grain_vol = count (ind_ph1(1, :, Nint(tmp*iy), :) &
                    & == i)
                     If (grain_vol > 0) Then
                        Write (11, '(2(f15.6))') (grain_vol*dms*2.0/pi) &
                       & ** (1./dms) * width
                        k = k + 1
                     End If
                     grain_vol = count (ind_ph1(1, Nint(tmp*ix), :, :) &
                    & == i)
                     If (grain_vol > 0) Then
                        Write (11, '(2(f15.6))') (grain_vol*dms*2.0/pi) &
                       & ** (1./dms) * width
                        k = k + 1
                     End If
                  End Do
               End Do
            End If
         End If
!
!
!
      End Subroutine
!
      Subroutine structure_out (filename)
         Character (Len=80) :: filename
         Integer i, j, k
         Open (9, File=filename)
         If (dms > 1) Then
            Write (9, "(A)") '<?xml version="1.0"?>'
            Write (9, "(A)") '<VTKFile type="ImageData" version="0.1" b&
           &yte_order="LittleEndian">'
            Write (9,*) '<ImageData WholeExtent="1 ', ix, ' 1 ', iy, ' &
           &1 ', iz, '" Origin="0 0 0" Spacing="', width, ' ', width, '&
           & ', width, '">'
            Write (9,*) '<Piece Extent="1 ', ix, ' 1 ', iy, ' 1 ', iz, &
           & '">'
            Write (9,*) '<PointData Scalars="energy">'
            Write (9,*) '<DataArray type="Int16" Name="korn" format="as&
           &cii">'
!
            Do k = 1, iz
               Do j = 1, iy
                  Do i = 1, ix
!                                  write(9,"(I3.1)",advance="no")count(phase_index(ind_ph1(:,i,j,k))==1)
!
!
!                                if(ind_ph1(2,i,j,k)>0)then
!                                   if(val_phase(2,i,j,k)>phasemin.and. val_phase(3,i,j,k)>phasemin.and.val_phase(1,i,j,k)>phasemin)
!                                            write(9,"(I2.1)",advance="no")6
!                                   elseif(val_phase(2,i,j,k)>phasemin.and. val_phase(3,i,j,k)>phasemin) then
!                                            write(9,"(I2.1)",advance="no")7
!                                   elseif (val_phase(2,i,j,k)>phasemin.and. val_phase(1,i,j,k)>phasemin)then
!                                            write(9,"(I2.1)",advance="no")8
!                                    else
!                                         write(9,"(I2.1)",advance="no")9
!                                    endif
                     If (ind_ph1(2, i, j, k) > 0 .And. val_ph1(1, i, j, &
                    & k) < 0.9) Then
                        Write (9, "(I2.1)", Advance="no") 0
                     Else If (ind_ph1(1, i, j, k) <= stoich_n) Then
                        Write (9, "(I2.1)", Advance="no") 1
                     Else If (ind_ph1(1, i, j, k) <= npp0) Then
                        Write (9, "(I2.1)", Advance="no") 2
                     Else If (ind_ph1(1, i, j, k) <= npp1) Then
                        Write (9, "(I2.1)", Advance="no") 3
                     Else If (ind_ph1(1, i, j, k) <= tpp) Then
                        Write (9, "(I2.1)", Advance="no") 4 !0.4+0.5*rex_energy(ind_ph1(1,i,j,k))/maxval(rex_energy)
                     Else
                        Write (9, "(I2.1)", Advance="no") 5
                     End If
!
!
                  End Do
               End Do
            End Do
!
            Write (9,*)
            Write (9, "(A)") '</DataArray>'
            Write (9,*) '<DataArray type="Float32" Name="conc" format="&
           &ascii">'
         End If
         Do k = 1, iz
            Do j = 1, iy
               Do i = 1, ix
                              !  if(ind_ph1(2,i,j,k)>0)then
                              !          write(9,*)1.
                              !  else
                  Write (9, "(F7.3)", Advance="no") dot_product &
                 & (val_con1(:, i, j, k), val_phase(:, i, j, k))
                              !  endif
!
               End Do
            End Do
         End Do
         If (dms > 1) Then
            Write (9,*)
            Write (9, "(A)") '</DataArray>'
            Write (9, "(A)") '</PointData>'
            Write (9, "(A)") '</Piece>'
            Write (9, "(A)") '</ImageData>'
            Write (9, "(A)") '</VTKFile>'
         End If
         Write (*,*) 'Microstructure is output.'
         Write (12,*) 'Microstructure is output.'
         Close (9)
      End Subroutine
!
      Subroutine GSIZE2D
         Real T, area_fraction, volume_fraction, intf_aa, intf_ab
         Integer i, j, k, grain_area, tmp
         Character (Len=10) :: String_time
         Logical STAT
         command ((len_trim(foldername)+1) :) = '/grain_size' // trim &
        & (adjustl(String_time)) // 's.txt'
!
         Close (11)
         Open (11, File=command)
         Write (11,*) '  Diameter (um)  '
         Do i = npp0 + 1, tpp
            grain_area = count (ind_ph1(1, :, :, :) == i)
            Write (11, '((f15.6))') (grain_area*dms*2.0/pi) ** &
           & (1./(1.*dms)) * width
         End Do
         Close (11)
!
         If (iz > 1) Then
            command ((len_trim(foldername)+1) :) = '/2dcut_grain_size' &
           & // trim (adjustl(String_time)) // 's.txt'
            Open (11, File=command)
            Write (11,*) '    Diameter-via-area (um)'
            k = 0
!
            Do while (k <  500)
               tmp = random_int (1, iz)
               Do i = npp0 + 1, tpp
                  grain_area = count (ind_ph1(1, :, :, tmp) == i)
                  If (grain_area > i_intf**2) Then
                     Write (11, '(f15.6)') Sqrt (grain_area*4.0/pi) * &
                    & width
                     k = k + 1
                  End If
               End Do
            End Do
!
         End If
!
         Close (11)
      End Subroutine
End Module
