Module PRE_PROCESS
!!THIS MODULE IS TO READ IN THE VALUES OF ALL THE PARAMETERS FROM INPUT FILE;
!!AND TO CONSTRUCT THE INITIAL STRUCTURE BY TESSELLATION;
!!AND TO ASSIGN A STORED ENERGY VALUE,CONCENTRATION,TO EACH GRAIN (PHASE);
!!AND TO GET DIFFUSE INTERFACE;
!
      Use UTIL_FUNC
  
      Implicit None
Contains

! THE MAIN SUBROUTINE WHICH CALLS THE OTHER SUBROUTINES IN THE FILE
! READ THE INPUT FILE, INITIALIZE VARAIBLES, POINTERS, CONSTRUCT INITIAL MICROSTRUCTURES.
      Subroutine INPUT_AND_INITIALIZE
! SYMBOLS NOTATION
! REDUCTION	USED TO READ THE DEFORMATION REDUCTION OF THE MICROSTRUCTURE
! NUM_FERRITE	NUMBER OF FERRITE GRAINS
! P_FRAC 	FRACTION OF PEARLITE IN THE INITIAL MICROSTRUCTURE
!OTHERS		SOME VARIABLES FOR TEMPERARY USE
         Implicit None
         Integer Ini_random1, Ini_random2, Seedput (8), I, J
         Integer Num_ferrite
         Real Reduction, Temp, P_frac
         
         integer number_input_arg
         character(len=32)::input_file
         integer name_length
         integer error_code
!
!************************
         number_input_arg=COMMAND_ARGUMENT_COUNT() 
         if(number_input_arg<1)then
             write(*,*)"Please input a file"
             stop
         endif
         call GET_COMMAND_ARGUMENT(1,input_file,name_length,error_code)
         if(error_code>0)then
             write(*,*) "cannot open the input file"
             stop
         elseif(error_code<0)then
             write(*,*)"file name is too long"
             stop
         endif
         Open (8, File=input_file)
         Read (8,*)
         Read (8,*) Foldername
         Read (8,*)
         Read (8,*)
         Read (8,*)
         Read (8,*) Ix
         Read (8,*)
         Read (8,*) Iy
         Read (8,*)
         Read (8,*) Iz
         Read (8,*)
         Read (8,*) Width
         Read (8,*)
         Read (8,*) I_intf
         Read (8,*)
         Read (8,*) Phasemin
         Read (8,*)
         Read (8,*) Stoich_n, Num_ferrite
         Read (8,*)
         Read (8,*) P_frac
         Read (8,*)
         Read (8,*) Matrix_conc
         Read (8,*)
         Read (8,*) Stoich_conc
         Read (8,*)
         Read (8,*) Input_flag
         Read (8,*)
         Read (8,*) Final_t
         Read (8,*)
         Read (8,*) Temperature0
         T_path (1, :) = 1.E100
         T_path (1, 0) = 0.0
         T_path (2, 0) = Temperature0
         Read (8,*)
         Read (8,*) I
         Read (8,*) T_path (1:2, 1:I)
         T_path (2, :) = T_path (2, :) + 273.15
         Read (8,*)
         Read (8,*) Reduction
         Read (8,*)
         Read (8,*) N_iteration
         If (N_iteration < 10) Write (*,*) "Warning: Number Of Iteratio&
        &Ns Is Insufficient. At Leat Ten Is Recommended."
         Read (8,*)
         Read (8,*) Write_t
         Read (8,*)
         Read (8,*) Micro_out
         Read (8,*)
         Read (8,*) Ini_random1
         Read (8,*)
         Read (8,*) Ini_random2
         Read (8,*)
         Read (8,*) Radius
         Read (8,*)
         Read (8,*) Energy_file
         Read (8,*)
         Read (8,*) Nucleation_type
         Read (8,*)
         Read (8,*) Alpha_nuc_n
         Read (8,*)
         Read (8,*) Theta_star
         Read (8,*)
         Read (8,*) Gamma_nuc_t, Gamma_nuc_n, Aus_nuc_ferrite_t, &
        & Aus_nuc_ferrite_rate
         Read (8,*)!!!//Ferrite-Austenite
         Read (8,*) Mobl0 (1, 2), Mobl_q (1, 2)
         Read (8,*)!!!//Ferrite-Carbide
         Read (8,*) Mobl0 (1, 3), Mobl_q (1, 3)
         Read (8,*)!!!//Austenite-Carbide
         Read (8,*) Mobl0 (2, 3), Mobl_q (2, 3)
         Read (8,*)!!!///Ferrite-Ferrite
         Read (8,*) Mobl0 (1, 1), Mobl_q (1, 1)
         Read (8,*)!!!///Austenite-Austenite
         Read (8,*) Mobl0 (2, 2), Mobl_q (2, 2)
         Read (8,*)!!!///Carbide-Carbide
         Read (8,*) Mobl0 (3, 3), Mobl_q (3, 3)
         Read (8,*)!!!//Ferrite-Austenite
         Read (8,*) Intf_e (1, 2)
         Read (8,*)!!!//Ferrite-Carbide
         Read (8,*) Intf_e (1, 3)
         Read (8,*)!!!//Austenite-Carbide
         Read (8,*) Intf_e (2, 3)
         Read (8,*)!!!///Ferrite-Ferrite
         Read (8,*) Intf_e (1, 1)
         Read (8,*)!!!///Austenite-Austenite
         Read (8,*) Intf_e (2, 2)
         Read (8,*)!!!///Carbide-Carbide
         Read (8,*) Intf_e (3, 3)
!	 DIFFUSIVITY PRE-FACTOR AND ACTIVATION ENERGY
         Do I = 1, 2
            Read (8,*)
            Read (8,*) Diff0 (I), Diff_q (I)
         End Do

         Ref_c_para = 0.0
         Slope_para = 1.0
         Ref_c_nple = 0.0

         Ref_c_nple_r = 0.

!

         Partition_k_para = 1.

!
         Read (8,*)!!!//Solute Drag Model
         Read (8,*) Binding_energy, V_max0, Q_sd
         Close (8)
!
!************************
!***********************
!!!##########initialize Parameters###########
         If (Min(Ix, Iy, Iz) < 1) Then
            Write (*,*) "Input Error Of Domain Size"
            Stop
         End If
!
!
         Do I = 1, 3
            Do J = I, 3
               Mobl0 (J, I) = Mobl0 (I, J)
               Mobl_q (J, I) = Mobl_q (I, J)
               Intf_e (J, I) = Intf_e (I, J)

            End Do
         End Do

         Temperature = Temperature0

  !//Make Directory For Output
         Command = 'mkdir '
         Command (7:) = Foldername
         Call System (Command)
!
         Eta = Width * I_intf
         Factor1 = 1. / Width ** 2
         Factor2 = Pi ** 2 / 2. / Eta ** 2
         Factor3 = Pi / Eta
         Domain_size = (/ Ix, Iy, Iz /)
         Time = 0.
!
         Flag_diffusion = .False.
         Flag_aus_nuc_pearlite = .True. !!!!Flag Of Aus Nucleation In Pearlite
         Flag_aus_nuc_ferrite = .True.
!
         Dms = Count (Domain_size > 1)!!!Dimension
         Nngb = 2 * Dms !!!!Number Of Nearest Neighbors
!!!!Allocate Matrix
         Allocate (Ind_ph1(Npp, Ix, Iy, Iz), Val_ph1(Npp, Ix, Iy, Iz))
         Allocate (Ind_ph2(Npp, Ix, Iy, Iz), Val_ph2(Npp, Ix, Iy, Iz))
         Allocate (Val_con1(N_phase, Ix, Iy, Iz), Val_con2(N_phase, Ix, &
        & Iy, Iz), Val_phase(N_phase, Ix, Iy, Iz))
         Allocate (Flag(Ix, Iy, Iz), If_potential_site(Ix, Iy, Iz), &
        & If_interface(Ix, Iy, Iz), If_interface2(Ix, Iy, Iz))!!!Nucleation Matrix--Boolian
         If_potential_site = .True. !!!! Flag Of Potential Nucleation Sites
         Val_ph1 = 0.
         Val_ph2 = 0.
         Ind_ph1 = 0
         Ind_ph2 = 0
         Val_con1 = 0.
         Val_con2 = 0.
         Val_phase = 0.
         If_interface = .False.
         If_interface2 = .False.
 !!!Solute Drag Model Appended
         Allocate (Val_rate(Ix, Iy, Iz))
         Val_rate = Naf

 !!!
!
!
!
!!!############################################
!
!!!#########INITIALIZE COMPUTATION DOMAIN######
!// INITIALIZE RANDOM SEED
         Seedput = Ini_random1
         Call RANDOM_SEED (Put=Seedput)
         Call RANDOM_NUMBER (Temp)
!//Voronoi Tesselation
!
         If (Input_flag == 0) Then
            Npp0 = Stoich_n
            Call INI_STRUCT (Reduction)
            Call P_CONSTRUCT (P_frac)!!!Select Pearlite
!
   !  Call Write_grain_data(-2.,.True.)
            Npp0 = Num_ferrite + Stoich_n
            Tpp = Npp0
            Npp1 = Npp0 + Alpha_nuc_n
            Call FERR_CONSTRUCT (Reduction)
!
            Allocate (Rex_energy(Npp0))!!!Stored Energy Of Deformed Ferrite
            Allocate (Coherent_grain(Npp1+1:10000, 2))
     !//Rex_energy Randomized
            Call REX_EXP ()
         Else
            Npp0 = Num_ferrite + Stoich_n
            Call READ_STRUCT ()
!
         End If
!
!//Diffuse Interface At Boundaries
         Call INITIALIZE (N_iteration)
         Write (*, "('Fraction Of Stoichiometric Phase Is',F7.3)") Sum &
        & (Val_ph1, Ind_ph1 <= Stoich_n) / Product (Domain_size)
         Call CON_INI
!
!
!
!
!
!!!#################################
!
!//Re-Initialize For Nucleation
         Seedput = 0
         Seedput = Ini_random2
         Call RANDOM_SEED (Put=Seedput)
         Call RANDOM_NUMBER (Temp)!!!!Get Rid Of The Beginning Data, Usually Not Good Randomness
!
!
!
!!!Results File Handle
         Command = Foldername
         Command ((Len_trim(Foldername)+1) :) = '/frac_time.txt'
         Open (10, File=Command)
         Command = Foldername
         Command ((Len_trim(Foldername)+1) :) = '/nuclei.txt'
         Open (11, File=Command)
         Command ((Len_trim(Foldername)+1) :) = '/simulation.log'
         Open (12, File=Command)
!
         Write (10,*) '       Time(S)    Fraction Of Rex_alpha    Fract&
        &Ion Of Gamma    Fraction Of Stoich'
         Write (11,*) '       Time(S)  X(Grid Point)  Y(Grid Point)  Z(&
        &Grid Point)'
        
        Call T_PROFILE()
        Call WRITE_GRAIN_DATA (0., .True.)
      End Subroutine
!
! ASSIGN EACH GRAIN A STORED ENERGY BASED ON THE DISTRIBUTION IN STORE.TXT
! USING A COMBINATORIAL OPTIMIZATION SIMILAR TO SIMULATED ANNEALING.
      Subroutine Rex_exp ()
         Real Grain_vol (Npp0), X_energy (E_level), F_targt (E_level), &
        & Cumu (E_level), F_sampl (E_level), Generator, Error
         Integer I, J, K, Ind_energy (Npp0), Ii
!
 !!Cc// Now Rex_energy Is The Volume Of Each Grain
         Forall (I=1:Npp0) Rex_energy (I) = Count (Ind_ph1(1, :, :, :) &
        & == I) * 1.0
!
!!Cc// Read The Energy Distribution Density And Get The Cumulative Function
         Open (9, File=Energy_file)
 !Read(9,*)!E_level
!
         Do I = 1, E_level
            Read (9,*) X_energy (I), F_targt (I)
         End Do
         Close (9)
         F_targt = F_targt / Sum (F_targt)
         Do I = 1, E_level
            Cumu (I) = F_targt (I)
            If (I > 1) Cumu (I) = Cumu (I-1) + Cumu (I)
         End Do
!!//Assign A Energy No. To Each Grain
         Do I = 1, Npp0
            Call RANDOM_NUMBER (Generator)
            Ind_energy (I) = Minloc (Cumu, 1, Cumu >= Generator)
         End Do
!!Cc// Calculate The Volume Fraction Of Each Energy Level
         Do I = 1, E_level
            F_sampl (I) = Sum (Rex_energy((Stoich_n+1) :Npp0), &
           & Ind_energy == I)
         End Do
!
         F_targt = F_targt * Count (Ind_ph1(1, :, :, :) > Stoich_n)
!
         Do Ii = 1, 50000 !!!!! Simulated Annealing Algorithm
            K = Random_int (Stoich_n+1, Npp0)
            J = Random_int (1, E_level)
            I = Ind_energy (K)
            If (I == J) Cycle
            Error = F_sampl (J) - F_sampl (I) - F_targt (J) + F_targt &
           & (I) + Rex_energy (K)
!
            If (Error < 0.) Then
               Ind_energy (K) = J
               F_sampl (I) = F_sampl (I) - Rex_energy (K)
               F_sampl (J) = F_sampl (J) + Rex_energy (K)
            End If
         End Do
!
         Rex_energy = X_energy (Ind_energy) * 1.E-12 !!
         Write (*,*) "Optimization Of Stored Energy Assignment Is Done!&
        &"
         Write (12,*) "Optimization Of Stored Energy Assignment Is Done&
        &!"
!
         Open (9, File='rex_energy.txt')
         Write (9,*) '  Strain  Energy    ', '   Volume    ', '   Grain&
        & No.  '
         Do I = 1, Npp0
            Write (9,*) Rex_energy (I), Count (Ind_ph1(1, :, :, :) == &
           & I) * (Width**Dms), I
         End Do
         Close (9)
!
      End Subroutine
!
! CONSTRUCT A SINGLE-PHASE COLD-ROLLED MICROSTRUCTURE, SIMILAR TO A UNIFORM STREACH OF A PICTURE
      Subroutine Ini_struct (Reduction)
!
         Integer I, J, K, L, N
         Real Reduction
         Real Temp (Npp0), Temp2 (Npp0, 3), Ini_grain_size
         If (Npp0 == 0) Then
            Write (*,*) 'At Least 1 Grain Should Be Set'
            Write (12,*) 'At Least 1 Grain Should Be Set'
            Stop
         End If
         If (Npp0 == 1) Then
            Ind_ph1 (1, :, :, :) = 1
            Val_ph1 (1, :, :, :) = 1.
            Return
         End If
!
         Call RANDOM_NUMBER (Temp2)
         Temp2 (:, 1) = Ix * Width * Temp2 (:, 1)
         Temp2 (:, 2) = Iy * Width * Temp2 (:, 2)
         Temp2 (:, 3) = Iz * Width * Temp2 (:, 3)
!
         If (Iz == 1) Temp2 (:, 3) = 0.
         If (Iy == 1) Temp2 (:, 2) = 0.
         If (Ix == 1) Temp2 (:, 1) = 0.
         
!$omp Parallel Do &
!$omp Shared(Ind_ph1,Val_ph1,Temp2) &
!$omp Private(I,J,K,L,Temp)
!
         Do J = 1, Iy
            Do K = 1, Iz
               Do I = 1, Ix
!Dir$ Ivdep
                  Do L = 1, Npp0
                     Temp (L) = Dis1 (I, J, K, Temp2(L, 1), Temp2(L, &
                    & 2), Temp2(L, 3), 1.-Reduction)
                  End Do
                  Ind_ph1 (1, I, J, K) = Minloc (Temp, 1)
                  Val_ph1 (1, I, J, K) = 1.
               End Do
            End Do
         End Do
!$omp End Parallel Do
!
         Write (*,*) "Initial Microstructure Done By Tesselation!"
         Write (12,*) "Initial Microstructure Done By Tesselation!"
!
!
      End Subroutine
!
! CONSTRUCT A DUAL-PHASE COLD-ROLLED MICROSTRUCTURE
      Subroutine P_construct (Percent)
      !

!
         Real Percent
         Logical P_flag (Npp0), B_flag (Npp0)
         Integer P_index, P_pct, P_vol, Grain_vol, P_num, I, J, K, N, &
        & Jj
         Integer Band_spacing, Band_center, Band_num
         Integer, Pointer :: Temp (:)
         Band_spacing = nint(14.0 / Width)
         Band_num = Iy / Band_spacing
         P_pct = Ix * Iy * Iz * Percent
         P_vol = Ix * Iy * Iz
         P_flag = .True. !!!All Grains Are Pearlite Now
         Write (*,*) 'Optimization Of Volume Fraction Of Stoichiometric&
        & Phase Start.......'
         Write (12,*) 'Optimization Of Volume Fraction Of Stoichiometri&
        &C Phase Start.......'
!
  !!!!!!Banded Structure!!!!!
  !!!!!!Grains In Bands Are Selected!!!!!
         B_flag = .False. !!!No Grains Are In Bands
!
         Do J = nint(0.5 * Band_spacing), nint((Band_num-0.5) * Band_spacing), &
        & (Band_spacing)
!
            Do K = 1, Iz
!
!
               Do I = 1, Ix
                  Do N = 0, 1
                     B_flag (Ind_ph1(1, I, J+N, K)) = .True.
                  End Do
               End Do
            End Do
         End Do
!
!
!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!Randome Structure!!!!!!!!!!
 !B_flag=.True.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
         Do While (Abs(P_vol-P_pct) > P_pct*0.02)
            P_index = Random_int (1, Npp0)!!!!Randomly Selcet A Candidate
            If (B_flag(P_index) .And. (Count(B_flag) < Count(P_flag))) &
           & Cycle!!!.And.(Count(B_flag)<Count(P_flag))
            Grain_vol = Count (Ind_ph1(1, :, :, :) == P_index)!!!!Calculate The Volume
!
            If (P_flag(P_index)) Then
!
               If ((2*(P_vol-P_pct)-Grain_vol) > 0) Then !!If The Grain Is Pearlite And If Opting It Out Can Diminish The Error
                  P_flag (P_index) = .False. !!!!! This Grain Is Ferrite Now
                  P_vol = P_vol - Grain_vol
               End If
            Else If ((2*(P_vol-P_pct)+Grain_vol) < 0) Then
               P_flag (P_index) = .True.
               P_vol = P_vol + Grain_vol
            End If
!
         End Do
!
         Stoich_n = Count (P_flag)
!
         Do I = Stoich_n + 1, Npp0
            If (P_flag(I)) Then
               Do J = 1, Stoich_n
                  If ( .Not. P_flag(J)) Then
                     Where (Ind_ph1(1, :, :, :) == J) Ind_ph1 (1, :, :, &
                    & :) = Npp0 + 1
                     Where (Ind_ph1(1, :, :, :) == I) Ind_ph1 (1, :, :, &
                    & :) = J
                     Where (Ind_ph1(1, :, :, :) == Npp0+1) Ind_ph1 (1, &
                    & :, :, :) = I
                     P_flag (J) = .True.
                     P_flag (I) = .False.
                     Exit
                  End If
!
               End Do
            End If
         End Do
!
         Write (*, "('Fraction Of Stoichiometric Phase Is',F7.3)") &
        & P_vol * 1. / Product (Domain_size)
         Write (12, "('Fraction Of Stoichiometric Phase Is',F7.3)") &
        & P_vol * 1. / Product (Domain_size)
!
         Write (*,*) 'Number Of Stoichiometric Phase Particle Is', &
        & Stoich_n
         Write (12,*) 'Number Of Stoichiometric Phase Particle Is', &
        & Stoich_n
!
      End Subroutine
!
! CONSTRUCT A DUAL-PHASE COLD-ROLLED MICROSTRUCTURE-STAGE II 
      Subroutine Ferr_construct (Reduction)
!
         Integer N
!
         Integer I, J, K, L
!
         Real Temp ((Stoich_n+1) :Npp0), Temp2 ((Stoich_n+1) :Npp0, 3), &
        & Ini_grain_size, Reduction
         If (Npp0 == 0) Then
            Write (*,*) 'At Least 1 Grain Should Be Set'
            Write (12,*) 'At Least 1 Grain Should Be Set'
            Stop
         End If
         If (Npp0 == 1) Then
            Ind_ph1 (1, :, :, :) = 1
            Val_ph1 (1, :, :, :) = 1.
            Return
         End If
         If (Any(Ind_ph1(1, :, :, :) == 0)) Stop
!
         Call RANDOM_NUMBER (Temp2)
         Temp2 (:, 1) = Ix * Width * Temp2 (:, 1)
         Temp2 (:, 2) = Iy * Width * Temp2 (:, 2)
         Temp2 (:, 3) = Iz * Width * Temp2 (:, 3)
!
         If (Iz == 1) Temp2 (:, 3) = 0.
         If (Iy == 1) Temp2 (:, 2) = 0.
         If (Ix == 1) Temp2 (:, 1) = 0.
!$omp Parallel Do &
!$omp Shared(Ind_ph1,Val_ph1,Temp2) &
!$omp Private(I,J,K,L,Temp)
!
         Do J = 1, Iy
            Do K = 1, Iz
               Do I = 1, Ix
                  If (Ind_ph1(1, I, J, K) <= Stoich_n) Cycle
                  Do L = Stoich_n + 1, Npp0
                     Temp (L) = Dis1 (I, J, K, Temp2(L, 1), Temp2(L, &
                    & 2), Temp2(L, 3), 1.-Reduction)
                  End Do
!
                  Ind_ph1 (1, I, J, K) = Minloc (Temp, 1) + Stoich_n !!!!Minloc From 1 To Size(Temp).
                  Val_ph1 (1, I, J, K) = 1.
               End Do
            End Do
         End Do
!$omp End Parallel Do
         Write (*,*) 'Rex_construct Done'
         Write (12,*) 'Rex_construct Done'
!
      End Subroutine

! READ INITAL MICROSTRUCTURE FROM A FILE
      Subroutine Read_struct ()
         Integer I, J, K, N
         Open (8, File='micro.vti')
         Do I = 1, 7
            Read (8,*)
         End Do
         Do K = 1, Iz
            Do J = 1, Iy
               Do I = 1, Ix
                  Read (8,*) Ind_ph1 (1, I, J, K)
               End Do
            End Do
         End Do
         Close (8)
         Write (*,*) "Initial Microstructure Read In From File 'Micro.V&
        &Ti'"
         Write (12,*) "Initial Microstructure Read In From File 'Micro.&
        &Vti'"
!
         If (Maxval(Ind_ph1) > Tpp) Then
            Write (*,*) "Phase Index Flow Out"
            Stop
         End If
         Val_ph1 (1, :, :, :) = 1.0
!
         Open (8, File='rex_energy.txt')
         Do I = 1, Npp0
            Read (8,*) Rex_energy (I)
         End Do
!
         Close (8)
!
      End Subroutine
! 
! GENERATE DIFFUSE INTERFACES BY RUNNING SOME ITERATIONS
      Subroutine Initialize (N)
!
         Integer I, J, K, M, N
!
         Del_t = 0.5 / Dms * Width ** 2
         Do M = 1, N
!$omp Parallel Do &
!$omp Shared(Ind_ph1,Ind_ph2,Val_ph1,Val_ph2) &
!$omp Private(I,J,K)
            Do K = 1, Iz
!
               Do J = 1, Iy
                  Do I = 1, Ix
                     Call PURE_GROWTH (I, J, K)
                  End Do
               End Do
!
            End Do
!$omp End Parallel Do
            Ind_tmp => Ind_ph1
            Ind_ph1 => Ind_ph2
            Ind_ph2 => Ind_tmp
!
            Val_tmp => Val_ph1
            Val_ph1 => Val_ph2
            Val_ph2 => Val_tmp
!
         End Do
!
!
         Do K = 1, Iz
            Do J = 1, Iy
               Do I = 1, Ix
!
                  If (Val_ph1(1, I, J, K) < 1.-Phasemin) If_interface &
                 & (I, J, K) = .True.
!
!
               End Do
            End Do
         End Do
!
         If (Any(Val_ph1(2, :, :, :) > Phasemin)) Then
            Write (*,*) "Diffuse Interface Initialized"
            Write (12,*) "Diffuse Interface Initialized"
         Else
            Write (*,*) "No Diffuse Interface"
            Write (12,*) "No Diffuse Interface"
         End If
!
         Write (*,*) "The Average Ferrite Size Is", Avg_size (1)
         Write (12,*) "The Average Ferrite Size Is", Avg_size (1)
!
      End Subroutine
!
! SIMULATING PURE GRAIN GROWTH, CALLED BY SUBROUTINE INITIALIZE()
      Subroutine Pure_growth (I, J, K)
!
         Integer I, J, K, Ngb (3, 0:6), Ii, Jj, Kk, N, Error, Iix, Iiy, &
        & Iiz, Lst_name (Nmax), Phase1, Phase2
         Real Lst_phi (0:Nngb, Nmax), Lst_tmp (Nmax)
         Real Temp1
         Lst_name = 0
         Lst_phi = 0.
!
         Lst_tmp = 0.
         Error = 0
         N = 0
!
!!!!!!Get Neighbors
         Call Get_ngb1 (I, J, K, Ngb)
!
         Do Jj = 0, Nngb
!
            Iix = Ngb (1, Jj)
            Iiy = Ngb (2, Jj)
            Iiz = Ngb (3, Jj)
!
            Do Ii = 1, Count (Val_ph1(:, Iix, Iiy, Iiz) > Phasemin)
               Error = Match (Ind_ph1(Ii, Iix, Iiy, Iiz), Lst_name)
               If (Error == 0) Then
                  N = N + 1
                  Lst_name (N) = Ind_ph1 (Ii, Iix, Iiy, Iiz)
                  Error = N
               End If
               Lst_phi (Jj, Error) = Val_ph1 (Ii, Iix, Iiy, Iiz)
            End Do
         End Do
 !!!!Solve Pfm Equations
         Do Ii = 1, N
!
            Do Jj = 1, N
               If (Ii == Jj) Cycle
               If (Min(Lst_name(Jj), Lst_name(Ii)) <= Stoich_n) Cycle!!!!No Interation Between Ferrite And Pearlite
               Temp1 = Factor1 * (Lst_phi(0, Jj)*Sum(Lst_phi(0:Nngb, &
              & Ii))-Lst_phi(0, Ii)*Sum(Lst_phi(0:Nngb, Jj))) + Factor2 &
              & * (Lst_phi(0, Ii)-Lst_phi(0, Jj))
!
               Lst_tmp (Ii) = Lst_tmp (Ii) + Temp1
!
            End Do
!
            Lst_tmp (Ii) = Lst_tmp (Ii) * Del_t + Lst_phi (0, Ii)
!
            If (Lst_tmp(Ii) < Phasemin) Then
               Lst_tmp (Ii) = 0.
            Else If (Lst_tmp(Ii) > 1.) Then
               Lst_tmp (Ii) = 1.
            End If
!
         End Do
!
         Call Ssort (Lst_tmp, Lst_name, Nmax)
         Val_ph2 (:, I, J, K) = Lst_tmp (1:Npp) / Sum (Lst_tmp(1:Npp))
         Ind_ph2 (:, I, J, K) = Lst_name (1:Npp)
!
      End Subroutine
!
! ASSIGN CARBON CONCENTRATION FOR EACH POINTS
      Subroutine Con_ini ()
!
         Integer I, J, K, N
!
         Val_con1 (1, :, :, :) = Matrix_conc
         Val_con1 (2, :, :, :) = 0. !!!Meaningless
         Val_con1 (3, :, :, :) = Stoich_conc
         Val_phase = 0.
         Do K = 1, Iz
            Do J = 1, Iy
               Do I = 1, Ix
!
                  Do N = 1, N_phase
!
                     Val_phase (N, I, J, K) = Sum (Val_ph1(:, I, J, K), &
                    & Phase_index(Ind_ph1(:, I, J, K)) == N)
                  End Do
               End Do
            End Do
         End Do
         Write (*,*) 'Carbon Concentration Is Initilized.'
         Write (12,*) 'Carbon Concentration Is Initilized.'
         Write (*, "('The Average Solute Content Is',F7.3)") Sum &
        & (Val_phase*Val_con1) / Product (Domain_size)
         Write (12, "('The Average Solute Content Is',F7.3)") Sum &
        & (Val_phase*Val_con1) / Product (Domain_size)
         Write (*,*) "Phase Fraction", Sum (Val_phase(1, :, :, :)) / &
        & (Ix*Iy*Iz), Maxval (Val_phase(2, :, :, :)), Sum (Val_phase(3, &
        & :, :, :)) / (Ix*Iy*Iz)
         Domain_carbon = Sum (Val_phase*Val_con1)
      End Subroutine
!
End Module
