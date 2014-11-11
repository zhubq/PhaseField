!
!THIS IS A MODULE FOR MULTIPHASE FIELD MODEL.
!SINCE THERE ARE MANY FIELDS AT EACH POINT BUT ONLY A FEW (USUALLY FEWER THAN 6) WITH NONZERO VALUES,SO IT IS POSSIBLE TO STORE ONLY
!ACTUALLY IT IS BETTER TO DYNAMICALLY STORE THE NONZERO VALUES BUT THAT NEEDS DYNAMICALLY CHANGING AN ARRAY.
Module SPARSE
!! SOLVER() 			LOOP THE GRID POINTS IN THE DOMAIN TO SOLVE PHASE FIELD AND DIFFUSION EQUATIONS
!! UPDATE_PFM(I,J,K)		SOLVE PHASE FIELD EQUATIONS ON POINT(I,J,K)
!! PHASE_PARTITIONING(I,J,K)	PARTITION CARBON CONCENTRATION AMONG PHASES SUCH THAT QUASI-EQUILIBRIUM IS REACHED
!! UPDATE_DIFF(I,J,K)		SOLVE DIFFUSION EQUATIONS ON POINT(I,J,K)
!! DRIVING_FORCE(I,J,II,JJ,CI,CJ) RETURN DRIVING FORCE BETWEEN GRAIN I AND J, WITH PHASE INDEX II AND JJ, CARBON CONCNETRATION CI AND CJ
!! DISSIPATION_ITERATION	GET THE INTERFACE VELOCITY BY ITERATIVELY SOLVE V=DRIVING_FORCE-SOLUTE_DRAG(V) FOR AUSTENITE FORMATION
!! DISSIPATION_ITERATION_R	GET THE INTERFACE VELOCITY BY ITERATIVELY SOLVE V=DRIVING_FORCE-SOLUTE_DRAG(V) FOR FERRITE FORMATION
!! DISSIPATION_ITERATION2	AN DIFFERENT VERSION OF DISSIPATION_ITERATION, FOR BACKUP
!! DISSIPATION_ITERATION_R2	AN DIFFERENT VERSION OF DISSIPATION_ITERATION_R, FOR BACKUP
!!THERMO-DATA			UPDATE THE THERMODYNAMIC DATA FOR EACH TIME STEP SINCE THEY ARE TEMPERATURE-DEPENDENT
!! TIMESTEP			UPDATE 'TIME',AND TIMESTEP-'DEL_T' AND TIME/TEMPERATURE DEPEDENT PARAMETERS
!! PARAMETER_DRAG_PB		UPDATE TEMPERATURE-DEPENDENT PARAMTERS IN HILLERT(PURDY-BRECHET) SOLUTE DRAG MODEL
!! SOLUTE_DRAG_PB(V0, E0, De, U0)		RETURN THE SOLUTE DRAG FORCE GIVEN THE INPUT PARAMTERS
!! SOLUTE_DRAG_PB0(V0, E0, U_nple, U0)		RETURN THE SOLUTE DRAG FORCE GIVEN THE INPUT PARAMETERS
!!
      Use OMP_LIB
      Use UTIL_FUNC
      Implicit None
      Logical coh_flag, dif_fl
!
Contains
!
      Subroutine SOLVER ()
         Integer I, J, K, N, C1, C2, Cr, Cm
         Dif_fl = .False.
         If_interface2 = .False.
 !$omp Parallel Do &
 !$omp Shared(Ind_ph1,Ind_ph2,Val_ph1,Val_ph2,Val_con1,Val_con2,Val_phase,If_interface) &
 !$omp Private(I,J,K)
!
         Do J = 1, Iy
            Do K = 1, Iz
               Do I = 1, Ix
                  Call UPDATE_PFM (I, J, K)
!
               End Do
            End Do
         End Do
  !$omp End Parallel Do
!
         Ind_tmp => Ind_ph1
         Ind_ph1 => Ind_ph2
         Ind_ph2 => Ind_tmp
!
         Val_tmp => Val_ph1
         Val_ph1 => Val_ph2
         Val_ph2 => Val_tmp
!
         If_intf_tmp => If_interface
         If_interface => If_interface2
         If_interface2 => If_intf_tmp
!
!   Val_con1=Val_con1*Domain_carbon/Sum(Val_phase*Val_con1)
!
         If ( .Not. Flag_diffusion) Return
! If(All(Phase_index(Ind_ph1)/=1))Return
!
         Do N = 1, Time_ratio
!
    !$omp Parallel Do &
    !$omp Shared(Ind_ph1,Ind_ph2,Val_ph1,Val_ph2,Val_con1,Val_con2,Val_phase) &
    !$omp Private(I,J,K)
!
            Do J = 1, Iy
               Do K = 1, Iz
                  Do I = 1, Ix
!
                     Call Update_diff (I, J, K)
!
                  End Do
               End Do
            End Do
     !$omp End Parallel Do
!
            Val_tmp => Val_con1
            Val_con1 => Val_con2
            Val_con2 => Val_tmp
!
         End Do
!
! ENDIF
      End Subroutine
!
      Subroutine UPDATE_PFM (I, J, K)
         Integer I, J, K, Ngb (3, 0:6), Ii, Jj, Kk, N, N_nonzero, &
        & Error, Iix, Iiy, Iiz, Lst_name (Nmax), Phaseii, Phasejj, &
        & Lst_phase (Nmax), Co_index, P
         Real Lst_phi (0:Nngb, Nmax), Lst_tmp (Nmax), Driv_matrix &
        & (Nmax), Drag_matrix (Nmax)
         Real Totalc, Temp1, Temp2
         Logical Flag_coexist (2), Rex_intf, Coh_intf, Has_ferrite
         Real Drag_ij, Normal_vec, Driv_ij, Phase2_prev !!!Solute Drag Appended
!   Real Grad(3),Distance,Sum_v,Avg_v,Lambda1,Ref_tmp(N_phase,N_phase),Slope_tmp(N_phase,N_phase),Pc_tmp(N_phase,N_phase),Pk_tmp(N_p
!   Integer Avging_range,Num
!   Avging_range=I_intf
         Integer Avging_range
!   Avging_range=I_intf
         Real Avg_con (N_phase), Grad (3), Num, Distance
!
!
         If ( .Not. If_interface(I, J, K)) Then
            Val_ph2 (:, I, J, K) = Val_ph1 (:, I, J, K)
            Ind_ph2 (:, I, J, K) = Ind_ph1 (:, I, J, K)
            Val_rate (I, J, K) = 0.
!    Avg_velocity(I,J,K)=0.
            Return
         End If
         Driv_matrix = 0.
         Drag_matrix = 0.
         Lst_name = 0
         Lst_phi = 0.
         Lst_tmp = 0.
         Error = 0
         N = 0
         Lst_phase = 0
!
         Ind_ph2 (:, I, J, K) = 0
         Val_ph2 (:, I, J, K) = 0.
!
!
!
!!!!!!Get Neighbors
         Call Get_ngb1 (I, J, K, Ngb)
         N_nonzero = Count (Val_ph1(:, I, J, K) > Phasemin)
         N = N_nonzero
         Lst_name (1:N) = Ind_ph1 (1:N, I, J, K)
         Lst_phase (1:N) = Phase_index (Lst_name(1:N))
         Lst_phi (0, 1:N) = Val_ph1 (1:N, I, J, K)
!
         Flag_coexist = .False.
!
         Rex_intf = .False.
         Coh_intf = .False.
         Has_ferrite = .False.
         Co_index = 0
!
         If (Val_phase(1, I, J, K) > Phasemin) Then !Flag_coexist=.True.
            If (Val_phase(3, I, J, K) > Phasemin) Flag_coexist (1) = &
           & .True. !!Ferrite+Pearlite
            If (Val_phase(2, I, J, K) > Phasemin) Then
               Flag_coexist (2) = .True. !!Ferrite+Austenite
!
!
            End If
         End If
!
         Outer: Do Jj = 1, Nngb !!!!!!Retrieve All Non-Zero Variables
!
            Do Ii = 1, Count (Val_ph1(:, Ngb(1, Jj), Ngb(2, Jj), Ngb(3, &
           & Jj)) > Phasemin)
               Error = Match (Ind_ph1(Ii, Ngb(1, Jj), Ngb(2, Jj), &
              & Ngb(3, Jj)), Lst_name)
               If (Error == 0) Then
                  N = N + 1
                  If (N > Nmax) Exit Outer
                  Lst_name (N) = Ind_ph1 (Ii, Ngb(1, Jj), Ngb(2, Jj), &
                 & Ngb(3, Jj))
                  Lst_phase (N) = Phase_index (Lst_name(N))
                  Error = N
               End If
               Lst_phi (Jj, Error) = Val_ph1 (Ii, Ngb(1, Jj), Ngb(2, &
              & Jj), Ngb(3, Jj))
            End Do
         End Do Outer
!
!!!********************* Change 23/03/2012
         If (N == 1) Then !!!!Within One Grain
            Ind_ph2 (:, I, J, K) = Ind_ph1 (:, I, J, K)
            Val_ph2 (:, I, J, K) = Val_ph1 (:, I, J, K)
            Val_rate (I, J, K) = 0.
!
            Return
         End If
!!!!*************************
         Totalc = Dot_product (Val_phase(:, I, J, K), Val_con1(:, I, J, &
        & K))
! !!!!!! Average Carbon Concentration
!
         Avging_range = I_intf + 1
         Num = 0.
         Avg_con = 0.
!
!
         If (Flag_coexist(2)) Then
            Grad (1) = Val_phase (1, Myngb(I, Ix, 1), J, K) - Val_phase &
           & (1, Myngb(I, Ix,-1), J, K)
            Grad (2) = Val_phase (1, I, Myngb(J, Iy, 1), K) - Val_phase &
           & (1, I, Myngb(J, Iy,-1), K)
            Grad (3) = Val_phase (1, I, J, Myngb(K, Iz, 1)) - Val_phase &
           & (1, I, J, Myngb(K, Iz,-1))
            Grad = Grad / 2.
            Do Kk = - Min ((Iz-1)/2, Avging_range), Min ((Iz-1)/2, &
           & Avging_range)
               Iiz = Myngb (K, Iz, Kk)
               Do Jj = - Min ((Iy-1)/2, Avging_range), Min ((Iy-1)/2, &
              & Avging_range)
                  Iiy = Myngb (J, Iy, Jj)
                  Do Ii = - Min ((Ix-1)/2, Avging_range), Min &
                 & ((Ix-1)/2, Avging_range)
                     Iix = Myngb (I, Ix, Ii)
                   !!!!! Two Interfaces May Be Close To Each Other:
!
!
                     Distance = ((Jj*Grad(3)-Kk*Grad(2))**2+(Kk*Grad(1)-&
                    & Ii*Grad(3))**2+(Ii*Grad(2)-Jj*Grad(1))**2)
!
!
                     If (Distance > Sum(Grad**2)*1.0) Cycle
                     If (Val_phase(1, Iix, Iiy, Iiz) > 1.-Phasemin) &
                    & Has_ferrite = .True.
                     If (Any(Val_phase(1:2, Iix, Iiy, Iiz) < Phasemin)) &
                    & Cycle
                     Rex_intf = Rex_intf .Or. &
                    & (Any(Phase_index(Ind_ph1(:, Iix, Iiy, Iiz)) == 1 &
                    & .And. Ind_ph1(:, Iix, Iiy, Iiz) <= Npp0) .And. &
                    & Any(Phase_index(Ind_ph1(:, Iix, Iiy, Iiz)) == 1 &
                    & .And. Ind_ph1(:, Iix, Iiy, Iiz) > Npp0))
                     Do P = 1, Npp
                        If (Phase_index(Ind_ph1(P, Iix, Iiy, Iiz)) == &
                       & 2) Then
                           If (Any(Ind_ph1(:, Iix, Iiy, Iiz) == &
                          & Coherent_grain(Ind_ph1(P, Iix, Iiy, Iiz), &
                          & 1))) Coh_intf = .True.
                        End If
                     End Do
                     If (Dot_product(Grad, (/ Val_phase(1, Myngb(Iix, &
                    & Ix, 1), Iiy, Iiz)-Val_phase(1, Myngb(Iix, Ix,-1), &
                    & Iiy, Iiz), Val_phase(1, Iix, Myngb(Iiy, Iy, 1), &
                    & Iiz)-Val_phase(1, Iix, Myngb(Iiy, Iy,-1), Iiz), &
                    & Val_phase(1, Iix, Iiy, Myngb(Iiz, Iz, &
                    & 1))-Val_phase(1, Iix, Iiy, Myngb(Iiz, Iz,-1)) /)) &
                    & < 0.0) Cycle
                     If (Val_phase(1, Iix, Iiy, Iiz) < &
                    & (0.4*Val_phase(2, Iix, Iiy, Iiz)) .Or. &
                    & Val_phase(1, Iix, Iiy, Iiz) > (2.5*Val_phase(2, &
                    & Iix, Iiy, Iiz))) Cycle!!! Only Ferrit And Austenite
                     Avg_con = Avg_con + Val_con1 (:, Iix, Iiy, Iiz)
                     Num = Num + 1.
!
                  End Do
               End Do
            End Do
!
            Avg_con = Avg_con / Dble (Num)
            If (Nint(Num) == 0) Avg_con = Val_con1 (:, I, J, K)
!
!
         Else
            Avg_con = Val_con1 (:, I, J, K)
         End If
! Avg_con=Val_con1(:,I,J,K)
!
!
!
!
!!!!!!
!!!Solve Pfm Equations
         Do Ii = 1, N
!
            Phaseii = Lst_phase (Ii)
!
!
            Normal_vec = Factor3 * Sqrt (Lst_phi(0, Ii)*(1.-Lst_phi(0, &
           & Ii)))!Sqrt(Temp1)/2./Width
            Do Jj = 1, N
!
               Phasejj = Lst_phase (Jj)
!
               If (Ii == Jj .Or. Max(Lst_name(Ii), Lst_name(Jj)) <= &
              & Npp0 .Or. Abs(Phaseii-Phasejj) == 2) Cycle!!!No Interaction: 1) Same Grian, 2) Between Stoich And Ferrite, 3) Between Deformed Ferrites
           !  If(Flag_dissolve_pearlite==.False..And.Phaseii/=Phasejj.And.Any((/Phaseii,Phasejj/)==1))Cycle
               If (Flag_coexist(1) .And. (Phaseii+Phasejj == 3)) Cycle!!!Only Pearlite-Austenite Reaction
!              If(Flag_coexist(2).And.Phaseii==Phasejj)Cycle
               If (Flag_coexist(2) .And. Phaseii == Phasejj .And. &
              & Min(Lst_name(Ii), Lst_name(Jj)) > Npp0) Cycle!!! No Grain Growth If There Is Phase Transformation, This Causes Sudden Growth Immediately After Full Austenitization.
!             If(Phaseii+Phasejj==3.And.(Any(Lst_name(1:N_nonzero)<=Npp0).And.Any(Lst_name(1:N_nonzero)>Npp0.And.Lst_name(1:N_nonzer
!
               Driv_ij = (Factor1*(Lst_phi(0, Jj)*Sum(Lst_phi(0:Nngb, &
              & Ii))-Lst_phi(0, Ii)*Sum(Lst_phi(0:Nngb, &
              & Jj)))+Factor2*(Lst_phi(0, Ii)-Lst_phi(0, Jj))) * Intf_e &
              & (Phaseii, Phasejj) + Factor3 * Sqrt (Lst_phi(0, &
              & Ii)*Lst_phi(0, Jj)) * Driving_force (Lst_name(Ii), &
              & Lst_name(Jj), Phaseii, Phasejj, Avg_con(Phaseii), &
              & Avg_con(Phasejj))
!                          (Driving_force(Lst_name(Ii),Lst_name(Jj),Phaseii,Phasejj,Val_con1(Phaseii,I,J,K),Val_con1(Phasejj,I,J,K),
                        !!!!!Bug:N_phases(3) Is Not Defined!!!!
               If (Phaseii+Phasejj == 3) Then
!                        If(Coh_intf.And.Rex_intf)Cycle
                  Drag_ij = Factor3 * Sqrt (Lst_phi(0, Ii)*Lst_phi(0, &
                 & Jj))!!!!Solute Drag Model For Alpha And Gamma
               Else
                  Drag_ij = 0.
               End If
!
               If (((((Phaseii == 2 .And. &
              & Any(Coherent_grain(Lst_name(Ii), :) == Lst_name(Jj))) &
              & .Or. (Phasejj == 2 .And. &
              & Any(Coherent_grain(Lst_name(Jj), :) == &
              & Lst_name(Ii)))))) .And. &
              & Count(Phase_index(Ind_ph1(1:N_nonzero, I, J, K)) == 1) &
              & == 1 .And. ( .Not. Rex_intf) .And. Has_ferrite) Then !.And.Count(Phase_index(Ind_ph1(1:N_nonzero,I,J,K))==1)==1   Count(Phas
!
                  Driv_ij = 0.02 * Mobl (Phaseii, Phasejj) * Driv_ij
                  Drag_ij = 0.02 * Mobl (Phaseii, Phasejj) * Drag_ij
!
               Else
                  Driv_ij = Mobl (Phaseii, Phasejj) * Driv_ij
                  Drag_ij = Mobl (Phaseii, Phasejj) * Drag_ij
!
               End If
!
               Driv_matrix (Ii) = Driv_matrix (Ii) + Driv_ij
!              Driv_matrix(Jj)=Driv_matrix(Jj)-Driv_ij
               Drag_matrix (Ii) = Drag_matrix (Ii) + Drag_ij
!              Drag_matrix(Jj)=Drag_matrix(Jj)+Drag_ij
            End Do
        !!!! Dphi/Dt=Drag_i-Drag_j*Dragforce
        !!!!Interation To Solve Drag Force Equation.
!
            Lst_tmp (Ii) = Driv_matrix (Ii)
!
            If ( .Not. Flag_coexist(1) .And. Flag_coexist(2) .And. Ii &
           & <= N_nonzero) Then
!
               If (Driv_matrix(Ii) > 0.) Then
                  If (Phaseii == 2) Then
!                    If(Temperature<900.)Write(*,*)"Note:",Driv_matrix(Ii),Val_con1(2,I,J,K)
                     Call Dissipation_iteration (Lst_tmp(Ii), &
                    & Driv_matrix(Ii), Drag_matrix(Ii), Normal_vec, &
                    & Val_con1(2, I, J, K))
                  Else If (Phaseii == 1) Then
                     Call Dissipation_iteration_r (Lst_tmp(Ii), &
                    & Driv_matrix(Ii), Drag_matrix(Ii), Normal_vec, &
                    & Val_con1(2, I, J, K))
                  End If
!
!
               Else If (Driv_matrix(Ii) < 0.) Then
!
                  If (Phaseii == 1) Then
                     Call Dissipation_iteration &
                    & (Lst_tmp(Ii),-Driv_matrix(Ii), Drag_matrix(Ii), &
                    & Normal_vec, Val_con1(2, I, J, K))!!!!!   Newton_iteration(X0,X,P,K,Norm)
                  Else If (Phaseii == 2) Then
                     Call Dissipation_iteration_r &
                    & (Lst_tmp(Ii),-Driv_matrix(Ii), Drag_matrix(Ii), &
                    & Normal_vec, Val_con1(2, I, J, K))!!!!!   Newton_iteration(X0,X,P,K,Norm)
                  End If
                  Lst_tmp (Ii) = - Lst_tmp (Ii)
               Else
                  Lst_tmp (Ii) = 0.
!
               End If
!
            End If
!
!
            If (N == 2) Then
               Lst_tmp (2) = - Lst_tmp (1)
               Exit
            End If
!
         End Do
         Lst_tmp (1:N) = Lst_tmp (1:N) * Del_t + Lst_phi (0, 1:N)
!
         Where (Lst_tmp(1:N) < Phasemin)
            Lst_tmp (1:N) = 0.
            Elsewhere (Lst_tmp(1:N) > 1.)
            Lst_tmp (1:N) = 1.
         End Where
!
         Call Ssort2 (Lst_tmp, Lst_name, Lst_phi(0, :), Nmax)
         Val_ph2 (1:Npp, I, J, K) = Lst_tmp (1:Npp) / Sum &
        & (Lst_tmp(1:Npp))
         N = Count (Val_ph2(:, I, J, K) > Phasemin)!!!Number Of Nonzero Order Parameter
         Val_ph2 (1:N, I, J, K) = Val_ph2 (1:N, I, J, K) / Sum &
        & (Val_ph2(1:N, I, J, K))
!
         Ind_ph2 (1:N, I, J, K) = Lst_name (1:N)
!
         If (N > 1) Then
            Do Ii = 0, Nngb
               If_interface2 (Ngb(1, Ii), Ngb(2, Ii), Ngb(3, Ii)) = &
              & .True.
            End Do
         End If
!
!
!!!!!*********************************************************************************************************************
  !!Partition Totalc To Composition C
  !Partition Among Non-Zero Phases.
!
         Phase2_prev = Val_phase (2, I, J, K)
!
         Do Ii = 1, N_phase
!
            Val_phase (Ii, I, J, K) = Sum (Val_ph2(1:N, I, J, K), &
           & Phase_index(Ind_ph2(1:N, I, J, K)) == Ii)
         End Do
!
         If (Val_phase(2, I, J, K) > Phasemin .And. Val_phase(1, I, J, &
        & K) > Phasemin) Then
!
            Val_rate (I, J, K) = (Val_phase(2, I, J, K)-Phase2_prev) / &
           & (Del_t*Factor3*Sqrt(Val_phase(2, I, J, K)*(1.-Val_phase(2, &
           & I, J, K))))
!
         Else
            Val_rate (I, J, K) = 0.
!
         End If
!
         Do Ii = 1, N_phase
            Temp1 = Dot_product (Val_phase(:, I, J, K), &
           & Partition_c_para(:, Ii))
            Temp2 = Dot_product (Val_phase(:, I, J, K), &
           & Partition_k_para(:, Ii))
!
            Val_con1 (Ii, I, J, K) = (Totalc-Temp1) / Temp2
!
         End Do
!
!
!
!
!
      End Subroutine
!
      Subroutine PHASE_PARTITIONING (I, J, K)
         Integer I, J, K, Ii, Jj, Kk, Iix, Iiy, Iiz
         Real Grad (3), Distance, Sum_v, Avg_v, Lambda1, Ref_tmp &
        & (N_phase, N_phase), Slope_tmp (N_phase, N_phase), Pc_tmp &
        & (N_phase, N_phase), Pk_tmp (N_phase, N_phase)
         Real Totalc, Temp1, Temp2
         Integer Avging_range, Num
!
         Totalc = Dot_product (Val_phase(:, I, J, K), Val_con1(:, I, J, &
        & K))
!
         Do Ii = 1, N_phase
            Temp1 = Dot_product (Val_phase(:, I, J, K), &
           & Partition_c_para(:, Ii))
            Temp2 = Dot_product (Val_phase(:, I, J, K), &
           & Partition_k_para(:, Ii))
!
            Val_con1 (Ii, I, J, K) = (Totalc-Temp1) / Temp2
!
         End Do
      End Subroutine
      
      Subroutine UPDATE_DIFF (I, J, K)
!
         Integer I, J, K, Ngb (3, 0:6), Ii, Jj, Kk, Iix, Iiy, Iiz
         Real Totalc, Del_c, Temp1, Temp2
! Real Lambda1,Slope_tmp(N_phase,N_phase),Ref_tmp(N_phase,N_phase),Pc_tmp(N_phase,N_phase),Pk_tmp(N_phase,N_phase)
!
!
         Call Get_ngb1 (I, J, K, Ngb)
!
!
         Del_c = 0.
!
!
         Do Kk = 1, Nngb
            Iix = Ngb (1, Kk)
            Iiy = Ngb (2, Kk)
            Iiz = Ngb (3, Kk)
!
            Del_c = Del_c + Sum (Diff(1:2)*Sqrt(Val_phase(1:2, Iix, &
           & Iiy, Iiz)*Val_phase(1:2, I, J, K))*(Val_con1(1:2, Iix, &
           & Iiy, Iiz)-Val_con1(1:2, I, J, K)))
         End Do
!
         Totalc = Dot_product (Val_phase(1:N_phase, I, J, K), &
        & Val_con1(1:N_phase, I, J, K)) + Factor1 * Dt_diff * Del_c
!
         Do Ii = 1, N_phase
            Temp1 = Dot_product (Val_phase(:, I, J, K), &
           & Partition_c_para(:, Ii))
            Temp2 = Dot_product (Val_phase(:, I, J, K), &
           & Partition_k_para(:, Ii))
!
            Val_con2 (Ii, I, J, K) = (Totalc-Temp1) / Temp2
!
         End Do
!
!
!
!
!
      End Subroutine
!
      Real Function DRIVING_FORCE (I, J, II, JJ, CI, CJ)
!//I,J------PARAMETER INDEX;
!//II,JJ----PHASE INDEX;
!//CI,CJ----CONCENTRATION.
         Integer I, J, Ii, Jj
         Real Ci, Cj
!   Real,Optional::Lambda
!   If(.Not.Present(Lambda))Lambda=0.
         Driving_force = 0.
         If (I <= Npp0) Driving_force = - Rex_energy (I)
         If (J <= Npp0) Driving_force = Rex_energy (J)
!
!!!##########linear Phase Diagram###########
!
         If (Ii /= Jj .And. (Ii == 2 .Or. Jj == 2)) Then
!
            Driving_force = Driving_force - Ratio_k_para (Ii, Jj) * &
           & (Max(Ci, Cj)-Ref_c_para(2, 1))
! !      Driving_force=Driving_force-Entropy(Ii,Jj)*(0.5*(Slope_para(Ii,Jj)*(Ci-Ref_c_para(Ii,Jj))+Slope_para(Jj,Ii)*(Cj-Ref_c_para(
!
!  Elseif(Ii==3)Then
!      Driving_force=Entropy(Ii,Jj)*( Temperature-964.)!-Slope_para(Jj,Ii)*(Cj-Ref_c_para(Jj,Ii)))
!  Elseif(Jj==3)Then
!      Driving_force=Entropy(Ii,Jj)*( Temperature-964.)!-Slope_para(Ii,Jj)*(Ci-Ref_c_para(Ii,Jj)))
         End If
!
!!!########################################
!!!#####QUADRATIC APPROXIMATION############
!     if(max(ii,jj)<3)then
!
!         driving_force=driving_force-0.5*entropy(ii,jj)*(slope(ii,jj)*(Ci-ref_C(ii,jj))+slope(jj,ii)*(cj-ref_c(jj,ii)))
!     elseif(ii==3)then
!         driving_force=driving_force+entropy(ii,jj)*slope(jj,ii)*(cj-ref_c(jj,ii))
!     else
!         driving_force=driving_force+entropy(ii,jj)*slope(ii,jj)*(ci-ref_c(ii,jj))
!     endif
 !!!#######################################
!
      End Function
!
      Subroutine DISSIPATION_ITERATION (Solution, Driv, Drag, Norm, C)
!
         Real, Intent (In) :: Driv, Drag, Norm, C
         Real Solution, Dv, Tmp_v, Tmp_err, Err1, Err2, X1, X2, V
         Integer I, J
!
!
!
         Err2 = Driv - Drag * (Drag0(1))
!
         If (Err2 < 0.) Then
            Solution = 0.
            Return
         End If
!
! Solution=Err2
!     Return
!
!
!
!      Dv=0.02*Min(Driv/Norm,5.*V_max)
         Dv = 0.1 * V_max
         V = 0.
         Do J = 1, 2
            Do I = 1, 40
               Err1 = Err2
               V = V + Dv
!
               Err2 = Driv - V * Norm - Drag * (Solute_drag_pb(V, Ex0, &
              & De_heat(1), Ux0(1))*Correct_factor(1))!-(1.-Correct_factor(1))*Solute_drag_pb(V,0.,De_heat(1),Ux0(1) )&
!                      +Correct_factor(1)*Sum(Solute_drag_pb(I*Dv,0.,M_nple(1),(/1.887e-2,3.65e-3,3.12e-3/) )) )
!                                   +Correct_factor(1)*( Solute_drag_pb( V,0.,De_heat(2),Ux0(2) )+Solute_drag_pb( V,0.,De_heat(3),Ux
!
               If (Err1*Err2 <= 0.) Then
                  X1 = V - Dv
                  X2 = V
!
                  Do While ((X2-X1) > 0.01*X2)
                     Tmp_v = 0.5 * (X1+X2)
!
                     Tmp_err = Driv - Tmp_v * Norm - Drag * &
                    & (Solute_drag_pb(Tmp_v, Ex0, De_heat(1), &
                    & Ux0(1))*Correct_factor(1))!-(1.-Correct_factor(1))*Solute_drag_pb(Tmp_v,0.,De_heat(1),Ux0(1))&
!                                                 + Correct_factor(1)*Sum( Solute_drag_pb(Tmp_v,0.,M_nple,(/1.887e-2,3.65e-3,3.12e-3
!                                               +Correct_factor(1)*(Solute_drag_pb(Tmp_v,0.,De_heat(2),Ux0(2))+Solute_drag_pb(Tmp_v,
                     If (Err1*Tmp_err <= 0.) Then
                        X2 = Tmp_v
                        Err2 = Tmp_err
                     Else
                        X1 = Tmp_v
                        Err1 = Tmp_err
                     End If
                  End Do
                  Solution = 0.5 * (X1+X2) * Norm
                  Return
               End If
            End Do
            V = Driv / Norm
            Dv = - 0.2 * V
            Err2 = - 1.
         End Do
!
!
!  Write(*,*)"  T(C)         Carbon In Gamma   Chem Drivforce       Total Driv         V-Max"
! Write(*,*)Temperature-273.15,C ,Ratio_k_para(2,1)*(C-Ref_c_para(2,1)),Driv/Norm,Drag,Norm,V_max,Err1,Err2,X1,X2
! Stop "Iteration Doesn'T Converge"
!      If(Err2<0.)Then
!          Solution=0.
!          Else
!            Solution=Driv
!      Endif
      End Subroutine
!
      Subroutine DISSIPATION_ITERATION_R (Solution, Driv, Drag, Norm, &
     & C)
!
         Real, Intent (In) :: Driv, Drag, Norm, C
         Real Solution, Dv, Tmp_v, Tmp_err, Err1, Err2, X1, X2, V
         Integer I, J
!
!
!
!
!
         Err2 = Driv - Drag * (Drag0(2))
         If (Err2 < 0.) Then
            Solution = 0.
            Return
         End If
!  Solution=Err2
!     Return
!
!
!
         Dv = 0.1 * V_max
         V = 0.
         Do J = 1, 2
            Do I = 1, 40
               V = V + Dv
               Err1 = Err2
               Err2 = Driv - V * Norm - Drag * (Solute_drag_pb(V, Ex0, &
              & De_cool(1), Ux0(1))*Correct_factor(2))!-(1.-Correct_factor(2))*Solute_drag_pb(V,0.,De_cool(1),Ux0(1) )&
                            !      +Correct_factor(2)*( Solute_drag_pb(V,0.,De_cool(2),Ux0(2))+Solute_drag_pb(V,0.,De_cool(3),Ux0(3)
               If (Err1*Err2 <= 0.) Then
                  X1 = V - Dv
                  X2 = V
!
                  Do While ((X2-X1) > 0.01*X2)
                     Tmp_v = 0.5 * (X1+X2)
!
                     Tmp_err = Driv - Tmp_v * Norm - Drag * &
                    & (Solute_drag_pb(Tmp_v, Ex0, De_cool(1), &
                    & Ux0(1))*Correct_factor(2))!-(1.-Correct_factor(2))*Solute_drag_pb(Tmp_v,0.,De_cool(1),Ux0(1) )&
                                          !   +Correct_factor(2)*(Solute_drag_pb(Tmp_v,0.,De_cool(2),Ux0(2) )+Solute_drag_pb(Tmp_v,0
!
                     If (Err1*Tmp_err <= 0.) Then
                        X2 = Tmp_v
                        Err2 = Tmp_err
                     Else
                        X1 = Tmp_v
                        Err1 = Tmp_err
                     End If
                  End Do
                  Solution = 0.5 * (X1+X2) * Norm
                  Return
               End If
            End Do
            V = Driv / Norm
            Dv = - 0.2 * V
            Err2 = - 1.
         End Do
!
         Stop "Iteration Failed"
      End Subroutine
      Subroutine DISSIPATION_ITERATION2 (Solution, Driv, Drag, Norm, C)
!
         Real, Intent (In) :: Driv, Drag, Norm, C
         Real Solution, Dv, Tmp_v, Tmp_err, Err1, Err2, X1, X2, K_nple, &
        & Lambda
         Integer I, J, Minx
!
!
!
         Dv = 0.005 * Min (Driv/Norm, 50.*V_max)
!
         Err2 = Driv - Drag * (Drag0(1))
         If (Err2 <= 0.) Then
            Solution = 0.
            Return
         End If
!
         Do I = 1, 200
            Err1 = Err2
!
            Err2 = Driv - I * Dv * Norm - Drag * (Solute_drag_pb(I*Dv, &
           & Binding_energy, M_nple(1), &
           & Ux0(1))-(1.-Correct_factor(1))*Solute_drag_pb(I*Dv, 0., &
           & M_nple(1), Ux0(1))+Correct_factor(1)*(Solute_drag_pb(I*Dv, &
           & 0., M_nple(2), Ux0(2))+Solute_drag_pb(I*Dv, 0., M_nple(3), &
           & Ux0(3))))!-Lambda*(K_nple-Abs(Ratio_k_para(2,1)))*Abs(C-Ref_c_para(2,1)) )
            If (Err1 >= 0. .And. Err2 <= 0.) Then
               X1 = (I-1) * Dv
               X2 = I * Dv
!
               Do While ((X2-X1) > 0.01*X2)
                  Tmp_v = 0.5 * (X1+X2)
                  Tmp_err = Driv - Tmp_v * Norm - Drag * &
                 & (Solute_drag_pb(Tmp_v, Binding_energy, M_nple(1), &
                 & Ux0(1))-(1.-Correct_factor(1))*Solute_drag_pb(Tmp_v, &
                 & 0., M_nple(1), &
                 & Ux0(1))+Correct_factor(1)*(Solute_drag_pb(Tmp_v, 0., &
                 & M_nple(2), Ux0(1))+Solute_drag_pb(Tmp_v, 0., &
                 & M_nple(3), Ux0(1))))!-Lambda*(K_nple-Abs(Ratio_k_para(2,1)))*Abs(C-Ref_c_para(2,1)) )
!
                  If (Err1*Tmp_err <= 0.) Then
                     X2 = Tmp_v
                     Err2 = Tmp_err
                  Else
                     X1 = Tmp_v
                     Err1 = Tmp_err
                  End If
               End Do
               Solution = 0.5 * (X1+X2) * Norm
               Return
            End If
         End Do
!
         If (Err2 < 0.) Then
            Solution = 0.
         Else
            Solution = Driv
         End If
!
      End Subroutine
!
      Subroutine DISSIPATION_ITERATION_R2 (Solution, Driv, Drag, Norm, &
     & C)
!
         Real, Intent (In) :: Driv, Drag, Norm, C
         Real Solution, Dv, Tmp_v, Tmp_err, Err1, Err2, X1, X2, &
        & K_nple_r, Lambda
         Integer I, J, Minx
!
!
         Dv = 0.005 * Min (Driv/Norm, 50.*V_max)
!
!     Drag0=Solute_drag_pb(0.,Binding_energy,N_nple(1),1.887e-2)+Solute_drag_pb(0.,0.,N_nple(2),3.65e-3)+Solute_drag_pb(0.,0.,N_nple
!
         Err2 = Driv - Drag * (Drag0(2))
         If (Err2 <= 0.) Then
            Solution = 0.
            Return
         End If
!
         Do I = 1, 200
            Err1 = Err2
            Err2 = Driv - I * Dv * Norm - Drag * (Solute_drag_pb(I*Dv, &
           & Binding_energy, N_nple(1), &
           & Ux0(1))-(1.-Correct_factor(2))*Solute_drag_pb(I*Dv, 0., &
           & N_nple(1), Ux0(1))+Correct_factor(2)*(Solute_drag_pb(I*Dv, &
           & 0., N_nple(2), Ux0(1))+Solute_drag_pb(I*Dv, 0., N_nple(3), &
           & Ux0(3))))
            If (Err1 >= 0. .And. Err2 <= 0.) Then
               X1 = (I-1) * Dv
               X2 = I * Dv
!
               Do While ((X2-X1) > 0.01*X2)
                  Tmp_v = 0.5 * (X1+X2)
!
                  Tmp_err = Driv - Tmp_v * Norm - Drag * &
                 & (Solute_drag_pb(Tmp_v, Binding_energy, N_nple(1), &
                 & Ux0(1))-(1.-Correct_factor(2))*Solute_drag_pb(Tmp_v, &
                 & 0., N_nple(1), &
                 & Ux0(1))+Correct_factor(2)*(Solute_drag_pb(Tmp_v, 0., &
                 & N_nple(2), Ux0(2))+Solute_drag_pb(Tmp_v, 0., &
                 & N_nple(3), Ux0(3))))!-Lambda*(K_nple-Abs(Ratio_k_para(2,1)))*Abs(C-Ref_c_para(2,1)) )
!
                  If (Err1*Tmp_err <= 0.) Then
                     X2 = Tmp_v
                     Err2 = Tmp_err
                  Else
                     X1 = Tmp_v
                     Err1 = Tmp_err
                  End If
               End Do
               Solution = 0.5 * (X1+X2) * Norm
               Return
            End If
         End Do
!
         If (Err2 < 0.) Then
            Solution = 0.
         Else
            Solution = Driv
         End If
!
      End Subroutine
!
!
      Subroutine THERMO_DATA ()
!
!
         If (Temperature < 1086.52) Then
!     Ref_c_para(1,2)=Poly4(Temperature-1086.52,-1.1428e-4,-1.78883e-7,-3.66186e-10,2.53049e-12)
!     Ref_c_para(2,1)=Poly4(Temperature-1086.52,-0.00325,1.13712e-5,-2.88193e-8,1.05482e-10)
!
!     Slope_para(1,2)=1./Diff_poly4(Temperature-1086.52,-1.1428e-4,-1.78883e-7,-3.66186e-10,2.53049e-12)
!     Slope_para(2,1)=1./Diff_poly4(Temperature-1086.52,-0.00325,1.13712e-5,-2.88193e-8,1.05482e-10)
!
            Ref_c_para (1, 2) = Poly4 (Temperature-273.15,-8.309e-5, &
           & 0., 0., 0.) + 6.883e-2
            Ref_c_para (2, 1) = Poly4 (Temperature-273.15, &
           & 2.607e-2,-6.58e-5, 3.984e-8, 0.) + 0.9115
!
            Slope_para (1, 2) = 1. / Diff_poly4 &
           & (Temperature-273.15,-8.309e-5, 0., 0., 0.)
            Slope_para (2, 1) = 1. / Diff_poly4 (Temperature-273.15, &
           & 2.607e-2,-6.58e-5, 3.984e-8, 0.)
!
!
!
         Else
            Ref_c_para (1, 2) = Poly4 (Temperature-1086.52,-1.1428e-4, &
           & 0., 0., 0.)
            Ref_c_para (2, 1) = Poly4 (Temperature-1086.52,-0.00325, &
           & 0., 0., 0.)
!
            Slope_para (1, 2) = 1. / Diff_poly4 &
           & (0.,-1.1428e-4,-1.78883e-7,-3.66186e-10, 2.53049e-12)
            Slope_para (2, 1) = 1. / Diff_poly4 (0.,-0.00325, &
           & 1.13712e-5,-2.88193e-8, 1.05482e-10)
!
         End If
!
!
! !!!****$$$$%%%%%
! !!!!Nple Poly4=A*X+B*X*X+C*X*X*X+D*X*X*X*X
         If (Temperature < 1110.64) Then
 !!!!!!Nple
            Ref_c_nple (1, 2) = Poly4 &
           & (Temperature-1110.64,-1.293e-4,-2.5283e-7,-9.489e-10, 0.0)
            Ref_c_nple (2, 1) = Poly4 (Temperature-1110.64,-3.3441e-3, &
           & 1.0406e-5,-5.6345e-8, 0.0)
!             Slope_nple (1, 2) = - 9034. !66500.+50.*Temperature
!             Slope_nple (2, 1) = 1. / Diff_poly4 &
!            & (Temperature-1110.64,-3.3441e-3, 1.0406e-5,-5.6345e-8, &
!            & 0.0)
         Else
!!!!!Nple
!             Slope_nple (1, 2) = - 9034.
!             Slope_nple (2, 1) = 1. / Diff_poly4 (0.,-3.23154e-3, &
!            & 1.39693e-5, 5.73906e-9, 3.25344e-10)
!
            Ref_c_nple (1, 2) = Poly4 (Temperature-1110.64,-1.23698e-4, &
           & 0., 0., 0.)
            Ref_c_nple (2, 1) = Poly4 (Temperature-1110.64,-3.23154e-3, &
           & 0., 0., 0.)
         End If
!
         If (Temperature < 1041.3) Then
!    Ref_c_nple_r(2,1)=Poly4(Temperature-1041.3,-2.5194e-3,1.0661e-5,-2.3132e-8,0.)
!    Ref_c_nple_r(1,2)=Poly4(Temperature-1041.3,-6.E-5,-2.E-7,-9.43e-12,0.)
!    Slope_nple_r(1,2)=1./Diff_poly4(Temperature-1041.3,-6.E-5,-2.E-7,-9.43e-12,0.)
!    Slope_nple_r(2,1)=1./Diff_poly4(Temperature-1041.3,-2.5194e-3,1.0661e-5,-2.3132e-8,0.)
!
!             Ref_c_nple_r (2, 1) = Poly4 (Temperature-273.15,-6.605e-2, &
!            & 8.137e-5,-3.56e-8, 0.0) + 18.82
!             Ref_c_nple_r (1, 2) = Poly4 &
!            & (Temperature-1041.3,-2.E-7,-9.43e-12, 0., 0.) - 6.E-5
!             Slope_nple_r (1, 2) = 1. / Diff_poly4 &
!            & (Temperature-1041.3,-6.E-5,-2.E-7,-9.43e-12, 0.)
!             Slope_nple_r (2, 1) = 1. / Diff_poly4 &
!            & (Temperature-273.15,-6.605e-2, 8.137e-5,-3.56e-8, 0.)
!
!
         Else
            Ref_c_nple_r (1, 2) = - 6.15e-5 * (Temperature-1041.3)
            Ref_c_nple_r (2, 1) = - 2.39e-3 * (Temperature-1041.3)
!             Slope_nple_r (1, 2) = 1. / (-6.15e-5)
!             Slope_nple_r (2, 1) = 1. / (-2.518e-3)
         End If
!
!
! !!!!****####
         Ref_c_para (1, 3) = Ref_c_para (1, 2)
         Ref_c_para (3, 1) = Ref_c_para (2, 1)
         Slope_para (1, 3) = Slope_para (1, 2)
         Slope_para (3, 1) = Slope_para (2, 1)
!
         Ref_c_para (2, 3) = Ref_c_para (2, 1)
         Ref_c_para (3, 2) = Ref_c_para (2, 1)
         Slope_para (2, 3) = Slope_para (2, 1)
         Slope_para (3, 2) = Slope_para (2, 1)
         Partition_k_para = Transpose (Slope_para) / Slope_para
         Partition_c_para = Ref_c_para - Partition_k_para * Transpose &
        & (Ref_c_para)
!
         Ref_c_nple (1, 3) = Ref_c_nple (1, 2)
         Ref_c_nple (3, 1) = Ref_c_nple (2, 1)
!          Slope_nple (1, 3) = Slope_nple (1, 2)
!          Slope_nple (3, 1) = Slope_nple (2, 1)
!          Partition_k_nple = Transpose (Slope_nple) / Slope_nple
!          Partition_c_nple = Ref_c_nple - Partition_k_nple * Transpose &
!         & (Ref_c_nple)
!
         Ref_c_nple_r (1, 3) = Ref_c_nple_r (1, 2)
         Ref_c_nple_r (3, 1) = Ref_c_nple_r (2, 1)
!          Slope_nple_r (1, 3) = Slope_nple_r (1, 3)
!          slope_nple_r (3, 1) = slope_nple_r (2, 1)
!          partition_k_nple_r = transpose (slope_nple_r) / slope_nple_r
!          partition_c_nple_r = ref_c_nple_r - partition_k_nple_r * &
!         & transpose (ref_c_nple_r)
!
!   ref_c=avg_lambda*ref_c_nple+(1.-avg_lambda)*ref_c_para
!   slope=avg_lambda*slope_nple+(1.-avg_lambda)*slope_para
!   partition_K=transpose(slope)/slope
!   partition_C=ref_C-partition_k*transpose(ref_c)
!
!
         Ratio_k_para (1, 2) = (0.000264*(Temperature-273.15)**2-&
        & 0.4466*(Temperature-273.15)+243.4713) * 1.03e-12 !58.0e-12
         Ratio_k_para (2, 1) = - Ratio_k_para (1, 2)
         Ratio_k_para (3, 2) = Ratio_k_para (1, 2)
         Ratio_k_para (2, 3) = - Ratio_k_para (3, 2)
         Ratio_k_para (2, 2) = 0.
      End Subroutine
!
!
      Subroutine TIMESTEP ()
!
  !  real time_factor
         Real Dt1, Dt2
!
!
         Call Thermo_data ()
!
         Mobl = Mobl0 * Exp (-Mobl_q/(R*Temperature))!!!Mobility Matrix
!
!        Mobl(2,1)=Mobl(1,2)
         Diff0 (3) = 0.
         Diff = Diff0 * Exp (-Diff_q/(R*Temperature))!!!Diffusion Array
!
!
!
         Dt1 = Width ** 2 / Maxval (Mobl*Intf_e) / Nngb * 0.75
         If (Temperature < Gamma_nuc_t .And. Flag_diffusion) Dt1 = Dt1 &
        & * 0.02
!          Dt1 = Min (Dt1, 1./Heat_rate)
         Dt_diff = Width ** 2 / Maxval (Diff) / Nngb
         Dt_diff = Min (Dt_diff, Dt1*0.999)
!
         Time_ratio = Floor (Dt1/Dt_diff)
         Del_t = Time_ratio * Dt_diff
!
         Call Parameter_drag_pb()
!
!
      End Subroutine
!
!
      Subroutine PARAMETER_DRAG_PB ()
!
         Real Driv_f (2)
!
         V_max = V_max0 * Exp (-Q_sd/(R*Temperature))!!!!!Solute Drag Parameter
!
!
         M_nple (1) = Max (1.E-10, &
        & 3.9793e-5*(Temperature-273.15)-2.3163e-2)
         M_nple (2) = Max (1.E-10,-2.928e-8*(Temperature-273.15)**2+&
        & 5.629e-5*(Temperature-273.15)-2.35e-2)
         M_nple (3) = Max (1.E-10, 2.253e-8*(Temperature-273.15)**2-&
        & 4.255e-5*(Temperature-273.15)+2.351e-2)
!
         N_nple (1) = Max (1.E-10,-3.818e-9*(Temperature-273.15)**3+&
        & 9.055e-6*(Temperature-273.15)**2-7.408e-3*(Temperature-&
        & 273.15)+2.118)
         N_nple (2) = Max (1.E-10,-7.1227e-10*(Temperature-273.15)**3+&
        & 1.6307e-6*(Temperature-273.15)**2-1.2584e-3*(Temperature-&
        & 273.15)+3.3166e-1)
         N_nple (3) = 3.E-3
         Ux0 = (/ 1.887e-2, 3.65e-3, 3.12e-3 /)
!
!
         De_heat = - 0.5 * Log (M_nple*(1.-Ux0)/(Ux0*(1.-M_nple)))
         De_cool = - 0.5 * Log (N_nple*(1.-Ux0)/(Ux0*(1.-N_nple)))
!
!
!     De_heat=M_nple
!     De_cool=N_nple
         Ex0 = Binding_energy / (8.314*Temperature)
         Driv_f (1) = Abs (Ratio_k_para(2, 1)*(Ref_c_para(2, &
        & 1)-Ref_c_nple(2, 1)))
         Driv_f (2) = Abs (Ratio_k_para(2, 1)*(Ref_c_para(2, &
        & 1)-Ref_c_nple_r(2, 1)))
!
         Drag0 (1) = Solute_drag_pb (0., Ex0, De_heat(1), Ux0(1)) + &
        & Solute_drag_pb (0., 0., De_heat(2), Ux0(2)) + Solute_drag_pb &
        & (0., 0., De_heat(3), Ux0(3))
!
         Drag0 (2) = Solute_drag_pb (0., Ex0, De_cool(1), Ux0(1)) + &
        & Solute_drag_pb (0., 0., De_cool(2), Ux0(2)) + Solute_drag_pb &
        & (0., 0., De_cool(3), Ux0(3))
         Correct_factor (1) = Driv_f (1) / Drag0 (1)
         Correct_factor (2) = Driv_f (2) / Drag0 (2)
!
         Drag0 (1) = Solute_drag_pb (0., Ex0, De_heat(1), Ux0(1)) * &
        & Correct_factor (1)
         Drag0 (2) = Solute_drag_pb (0., Ex0, De_cool(1), Ux0(1)) * &
        & Correct_factor (2)
!     if(temperature<970.)drag0(1)=drag0(1)*10.!!!!*******&&&&&&%%%%^$#^$%^%$^&%^$&%^&%^&%^&%^&%^&%^
      End Subroutine
!
!
      Real Function SOLUTE_DRAG_PB (V0, E0, De, U0)
 !!!! V0-um/s: normalized
 !!!! T-K
 !!!! DE--J/MOL /RT
 !!!! E0--J/MOL /RT
 !!!! C0--MOL/um3
!
         Real, Intent (In) :: V0, E0, U0, De
         Real V, A, B, U1, D0, D1, D1_d0, D2_d1, Vu0
!
!
         V = V0 / V_max
         A = (De-E0) + V !
         B = (De+E0) + V !
         D0 = 1. - U0
         Vu0 = V * U0
         D1_d0 = (A-Vu0) / (A*D0+(A-V)*U0*Exp(-(A-Vu0)))
         D1 = D0 * D1_d0
         U1 = 1. - D1
         D2_d1 = (B-Vu0) / (B*D1+(B*U1-Vu0)*Exp(-(B-Vu0)))
!  D2=D1*D2_d1
!  U2=1.-D2
!
!
!
         Solute_drag_pb = (2.*De*U0-((A-V)/A*(Vu0+Log(D1_d0))+(B-&
        & V)/B*(Vu0+Log(D2_d1)))) * 8.314 * Temperature / 7.2e12
!
!
      End Function
!
!       Elemental REAL Function SOLUTE_DRAG_PB0 (V0, E0, U_nple, U0)
!  !!!! V0-um/s
!  !!!! T-K
!  !!!! DE--J/MOL
!  !!!! E0--J/MOL
!  !!!! C0--MOL/um3
! !
!          Real, Intent (In) :: V0, E0, U0, U_nple
!          Real De, V, A, B, Tmp1, Tmp2, Tmp3, U1, U2
! !
!          De = - 0.5 * Log (U_nple*(1.-U0)/(U0*(1.-U_nple))) * &
!         & (8.314*Temperature)
!          A = (De-E0) / (8.314*Temperature)
!          B = (De+E0) / (8.314*Temperature)
!          V = V0 / V_max
! !
!          A = A + V
!          B = B + V
!          U1 = 1. - (1.-U0) * (A-V*U0) / &
!         & (A*(1.-U0)+(A*U0-V*U0)*Exp(-(A-V*U0)))
!          U2 = 1. - (1.-U1) * (B-V*U0) / &
!         & (B*(1.-U1)+(B*U1-V*U0)*Exp(-(B-V*U0)))
! !
!          Solute_drag_pb0 = (2.*De*U0-((A-V)/A*(V*U0+Log((1.-U1)/(1.-&
!         & U0)))+(B-V)/B*(V*U0+Log((1.-U2)/(1.-U1))))*8.314*Temperature) &
!         & / 7.2e12
! !
!       End Function
     
!       Elemental REAL Function PARTITION_K_MN (V0, E0, Ratio_nple)
! !
!          Real, Intent (In) :: V0, E0, Ratio_nple !!!***Interface Velocity
!          Real Diff_intf, Diff_intf0, Diff_intf_q, Delta, M !!!Trans-Interface Diffusion, Interface Thickness,Bulk Conc Of Mn,Ferrite Side;
!          Real V, Md
!          Real A, B, De
!          V = V0 / V_max
!          De = - 0.5 * Log (Ratio_nple) * (8.314*Temperature)
!          A = (De-E0) / (8.314*Temperature)
!          B = (De+E0) / (8.314*Temperature)
!          M = 1. + A * Exp (-(V+B)) * (Exp(-(V+A))-1.) / (V+A) + B * &
!         & (Exp(-(V+B))-1.) / (V+B)
! !
!          Partition_k_mn = (M-1.) / (Ratio_nple-1.)
!          If (Isnan(Partition_k_mn)) Partition_k_mn = 1. !!!! At Low T, M_nple May Be Negative, So That Return Nan;
! !
!       End Function
 
!
End Module
!
!
