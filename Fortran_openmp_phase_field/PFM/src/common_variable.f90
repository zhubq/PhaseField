Module Common_variable
! GLOBAL OR COMMON VARIABLES SHARED BY ALL SUBROUTINES IN THE PROGRAM.
      Use Omp_lib

      Implicit None
      

! IX,IY,IZ               NUMBER OF GRID POINTS IN EACH DIMENSION, FOR 2D SIMULATION, SET IZ=1, FOR 1D SIMULATION, SET IY=1 AND IZ=1.
! DOMAIN_SIZE            A VECTOR, I.E. [IX,IY,IZ].
! I_INTF                 INTERFACE THICKNESS IN NUMBER OF GRID POINTS
! WIDTH                  GRID SPACING
! ETA                    INTERFACE THICKNESS IN MICRON, =WIDTH*I_INTF
! PHASEMIN               INTERFACE CUTOFF, VALUE BELOW IS REGARDED OUTSIDE INTERFACE
! NPP                    NUMBER OF PHASE FIELD VARIABLE STORED FOR EACH GRID POINTS
! NMAX                   MAXIMUM NONZERO PHASE FIELD VARAIBLES AT THE INTERFACE
! N_PHASE                NUMBER OF PHASES: FERRITE, PEARLITE,AUSTENITE
! DMS                    DIMENSIONALITY OF THE DOMAIN_SIZE 1,2 OR 3.
! NNGB                   NUMBER OF NEIGHBORS FOR A GRID POINTS,=2*DMS
! N_INTERATION           NUMBER OF INTERATIONS TO GENERATE DIFFUSE INTERFACE IN THE INITIALIZATION SUBROUTINES.
! TPP                    TOTAL NUMBER OF PHASE FIELD VARAIBLES
! NPP0                   NUMBER OF PHASE FIELD VARIABLES IN THE INITIAL MICROSTRUCTURE
! E_LEVEL                LEVELS OF STORED ENERGY IN THE STORE.TXT FILE, USED TO ASSIGN A LEVEL OF STORED ENERGY TO EACH GRAIN.
! PI                     THE PI CONSTANT,3.141592654
! NAF                    INFINITELY NEGATIVE NUMBER
! R                      GAS CONSTANT
! IND_PH1                4-DIMENSIONAL MATRIX, STORES THE INDICE OF PHASE FIELD VARIABLES FOR A 3D GRID,AND STORES NPP INDICE ON EACH GRID POINT.
! IND_PH2                AN ALTERNATIVE OF IND_PH1 TO SAVE UPDATED IND_PH1 IN EACH ITERATION.
! IND_TMP                A TEMPERARY POINTER TO SWITH IND_PH1 AND IND_PH2
! VAL_PH1                4-DIMENSIONAL MATRIX, STORES THE VALUES OF PHASE FIELD VARIABLES IN 3D SPACE.
! VAL_PH2                AN ALTERNATIVE OF VAL_PH1 TO SAVE UPDATED VAL_PH1 IN EACH ITERATION.
! VAL_TMP                A TEMPERARY POINTER TO SWITH VAL_PH1 AND VAL_PH2
! VAL_PHASE              4-DIMENSIONAL MATRIX, STORES VALUES OF PHASES
! VAL_CONC1              4-DIMENSIONAL MATRIX, STORES CARBON CONCENTRATION OF PHASES
! VAL_CON2               AN ALTERNATIVE OF VAL_CON1
! IF_POTENTIAL_SITE      POTENTIAL SITES FOR NUCLEATION
! IF_INTERFACE           LABEL OF THE INTERFACE
! IF_INTERFACE2          ALTERNATIVE OF IF_INTERFACE
! IF_INTF_TMP            TEMPERARY POINTER TO SWITCH IF_INTERFACE WITH IF_INTERFACE2
! FLAG                   RECORD NUCLEI POSITION
! REX_ENERGY             STORE THE STORED ENERGY OF EACH GRAIN
! MATRIX_CONC            INITIAL CARBON CONCENTRATION IN FERRITE
! STOICH_CONC            CARBON CONCENTRATION IN PEARLITE
! AUSTE_CONC             CARBON CONCENTRATION IN AUSTENITE
! DEL_T                  TIME STEP FOR PHASE FIELD SOLVER
! DT_DIFF                TIME STEP FOR DIFFUSION SOLVER
! TIME_RATIO             NUMBER OF CALLING DIFFUSION SOLVER FOR EACH CALL OF PHASE FIELD SOVLER.
! TIME                   TIME IN THE SIMULATION
! TEMPERATURE            TRANSIENT TEMPERATURE AT 'TIME' IN THE SIMULATION
! TEMPERATURE0           INITIAL TEMPERATURE
! FINAL_T                SIMULATION END TIME
! WRITE_T                TIME INTERVAL FOR DATA OUTPUT
! COHERENT_GRAIN         STORE THE INDICE OF FERRITE WITH WHICH THE AUSTENITE HAS A KS RELATIONSHIP
! STOICH_N               NUMBER OF PEARLITE PARTICLES
! NPP1                   NUMBER OF FERRITE+PEARLITE GRAINS
! DIFF                   VECTOR, CARBON DIFFUSIVITY IN EACH PHASE
! DIFF0                  VECTOR, PRE-FACTOR OF CARBON DIFFUSIVITY
! DIFF_Q                 VECTOR, ACTIVATION ENERGY OF CARBON DIFFUSIVITY
! MOBL                   MATRIX, INTERFACE MOBILITY
! MOBL0                  MATRIX, PRE-FACTOR OF INTERFACE MOBILITY
! MOBL_Q                 MATRIX, ACTIVATION ENERGY OF INTERFACE MOBILITY
! INTF_E                 MATRIX, INTERFACE ENERGY
! REF_C_PARA             MATRIX, LIQUIDUS AND SOLIDUS LINES IN THE PARA PHASE DIAGRAM
! SLOPE_PARA             MATRIX, SLOPES OF LIQUIDUS AND SOLIDUS LINES IN THE PARA PHASE DIAGRAM
! PARTITION_K_PARA       MATRIX, PARTITION COEFFICENTS BETWEEN EACH PAIR OF PHASES
! PARTITION_C_PARA       MATRIX, PARTITIONED CARBON BETWEEN EACH PAIR OF PHASES
! ENTROPY_PARA           MATRIX, ENTROPY BETWEEN EACH PAIR OF PHASES
! REF_C_NPLE             MATRIX, LIQUIDUS AND SOLIDUS LINES IN THE NPLE PHASE DIAGRAM FOR AUSTENITE FORMATION
! REF_C_NPLE_R           MATRIX, LIQUIDUS AND SOLIDUS LINES IN THE NPLE PHASE DIAGRAM FOR FERRITE FORMATION
! ALPHA_NUC_N            NUMBER OF NUCLEI FOR FERRITE RECRYSTALLIZATION
! GAMMA_NUC_N            NUMBER OF AUSTENITE NUCLEI IN PEARLITE
! GAMMA_NUC_T            START TEMPERATURE OF AUSTNEITE NUCLEATION IN PEARLITE
! GAMMA_NUC_FERRITE_T    START TEMPERATURE OF AUSTENITE NUCLEATION AT FERRITE GRAIN BOUNDARIES
! GAMMA_NUC_FERRITE_RATE AUSTENITE NUCLEATION RATE AT FERRITE GRAIN BOUNDARIES
! FLAG_AUS_NUC_PEARLITE  FLAG INDICATING IF NUCLEATION IN PEARITE IS FEASIBLE
! FLAG_AUS_NUC_FERRITE   FLAG INDICATING IF NUCLEATION AT FERRITE GRAIN BOUNDARIES IS FEASIBLE
! FLAG_DIFFUSION         FLAG INDICATING IF DIFFUSION SOLVER IS CALLED
! RADIUS                 RADIUS OF NUCLEI
! THETA_STAR             CRITICAL STORED ENERGY FOR RECRYSTALLIZATION NUCLEATION
! NUCLEATION_TYPE        HOMOGENEOUS OR HETERGENEOUS NUCLEATION
! INPUT_FLAG             FLAG INDICATING READ INITIAL MICROSTRUCTURE FROM A FILE OR CONSTRUCTED USING TESSLEATION.
! MICRO_OUT              FLAG INDICATING IF THE MICROSTRUCTURE IS OUTPUT AS *.VTI FILES
! FOLDERNAME             THE NAME OF FOLDER USED TO STORE OUTPUT FILES
! COMMAND                STRING USED TO MANIPULATE OUTPUT FILENAME
! ENERGY_FILE            THE NAME OF OUTPUT FILE STORING THE STORED ENERGY OF EACH GRAIN, USED TO CHECK IF THE STORED ENERGY DISTRIBUTION IS CONSISTENT WITH THE ONE IN STORE.TXT
! FACTOR1,2,3            NUMERICAL PARAMETERS USED TO IMPROVE COMPUTATIONAL EFFICIENCY
! VAL_RATE               3D MATRIX, USED TO STORE THE TRANSIENT INTERFACE VELOCITY
! V_MAX                  TRANS-INTERFACE DIFFUSIVITY/INTERFACE THICKNESS, USED IN THE SOLUTE DRAG MODEL DEVELOPED BY HILLERT.
! V_MAX0                 PRE-FACTOR OF V_MAX
! Q_SD                   ACTIVATION OF V_MAX
! BINDING_ENERGY         BINDING ENERGY OF MN
! M_NPLE                 MN NEGATIVE SPIKE IN FERRITE DURING AUSTENITE FORMATION UNDER NPLE
! N_NPLE                 MN POSITIVE SPIKE IN AUSTENITE DURING FERRITE FORMATION UNDER NPLE
! RATIO_K_PARA           RATIO CONSTANT BETWEEN DRIVING FORCE AND CARBON CONCENTRATION
! CORRECT_FACTOR         NUMERICAL FACTOR IN THE SOLUTE DRAG MODEL
! DRAG0                  SOLUTE DRAG FORCE FOR STATIONARY INTERFACES
! UX0                    CONCENTRATION OF SUBSTITUTIONAL ELEMENTS
! EX0                    BINDING ENERGY/RT
! DE_HEAT                SUBSTITUTIONAL POTENTIAL DIFERENCE BETWEEN AUSTENITE AND FERRITE DURING AUSTENITE FORMATION UNDER NPLE
! DE_COOL                SUBSTITUTIONAL POTENTIAL DIFFERENCE BETWEEN AUSTENITE AND FERRITE DURING FERRITE FORMATION UNDER NPLE
! DOMAIN_CARBON          AVERAGE CARBON CONCENTRATION IN THE DOMAIN_

      Integer Ix, Iy, Iz, Domain_size (3)
      Real I_intf, Width, Phasemin, Eta
!
      Integer Npp, Nngb, Nmax, Npp0, N_iteration, Dms, Tpp, E_level, &
     & N_phase
      Real Pi, Naf, R
      Parameter (Pi=3.141592654, Npp=6, Nmax=15, Naf=-1.E30, R=8.314, &
     & E_level=19, N_phase=3)
!     
      Pointer Val_ph1 (:, :, :, :), Val_ph2 (:, :, :, :), Val_tmp (:, &
     & :, :, :), Ind_ph1 (:, :, :, :), Ind_ph2 (:, :, :, :), Ind_tmp &
     & (:, :, :, :), Val_con1 (:, :, :, :), Val_con2 (:, :, :, :), &
     & Val_phase (:, :, :, :)
      Integer Ind_ph1, Ind_ph2, Ind_tmp
      Real Val_ph1, Val_ph2, Val_con1, Val_con2, Val_phase, Val_tmp
      
      Logical, Pointer :: If_potential_site (:, :, :), If_interface (:, &
     & :, :), If_interface2 (:, :, :), If_intf_tmp (:, :, :)      
      Real, Pointer :: Flag (:, :, :)     
      
      Real,Pointer:: Rex_energy (:)
      Real Matrix_conc, Stoich_conc, Auste_conc
!
      Real, Save :: Del_t, Time, Final_t, Write_t, &
     & Temperature, Temperature0, Dt_diff, T_path (2, 0:50)
      Integer Time_ratio
      Integer, Pointer :: Coherent_grain (:, :)!!!There Is Maximum 10000 Austenite Grains
      Integer Stoich_n, Npp1
      Real, Save :: Diff (N_phase), Diff0 (N_phase), Diff_q (N_phase), &
     & Mobl (N_phase, N_phase), Mobl0 (N_phase, N_phase), Mobl_q &
     & (N_phase, N_phase), Intf_e (N_phase, N_phase)
     
      Real, Save :: Ref_c_para (N_phase, N_phase), Slope_para (N_phase, N_phase), Partition_k_para (N_phase, &
     & N_phase), Partition_c_para (N_phase, N_phase),Entropy_para (N_phase, N_phase), Ref_c_nple (N_phase, N_phase), Ref_c_nple_r (N_phase, N_phase)
        
      Real, Save :: Gamma_nuc_t, Aus_nuc_ferrite_t, Aus_nuc_ferrite_rate
      Integer Alpha_nuc_n, Gamma_nuc_n
      Logical Flag_aus_nuc_pearlite, Flag_aus_nuc_ferrite 
      Logical Flag_diffusion

      Real Radius, Theta_star
      Integer Nucleation_type, Input_flag
      Logical Micro_out
!
      Character (Len=80) :: Foldername, Command, Energy_file

      Real Factor1, Factor2, Factor3
      

      Real, Pointer :: Val_rate (:, :, :)
      Real  V_max, V_max0, Q_sd, Binding_energy 
      Real M_nple (3), N_nple (3), Ratio_k_para (N_phase, N_phase), &
     & Correct_factor (2), Drag0 (2), Ux0 (3), Ex0, De_heat (3), &
     & De_cool (3)


      Real Domain_carbon
End Module