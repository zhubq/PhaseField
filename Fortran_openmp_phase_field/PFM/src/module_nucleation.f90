Module Module_nucleation
!!!!! THIS MODULE DEALS WITH NUCLEATION FOR FERRITE RECRYSTALLIZATION (INSTANTANEOUS) AND AUSTENITE (INSTANTANEOUS AT PEARITE/FERRITE INTERFACES AND FOLLOW
!!!!! CLASSICAL NUCLEATION MODEL AT FERRITE GRAIN BOUNDARIES)
!!!!!
!!!!!SUBROUTINE NOTATIONS
!!!!!
!!!!! Nucleation_rex_alpha(Nuc_n)		INSTATANEOUS NUCLEATION OF NUC_N NUCLEI FOR FERRITE RECRYSTALLIZATION
!!!!! HOMO(NUC_N)				INSTATANEOUS HOMGENEOUS NUCLEATION OF NUC_N NUCLEI FOR FERRITE RECRYSTALLIZATION
!!!!! CRITICAL(NUC_N)				INSTATANEOUS HETEROGENEOUS NUCLEATION OF NUC_N NUCLEI FOR FERRITE RECRYSTALLIZATION
!!!!! NUCLEATION_AUSTENITE()			AUSTENITE NUCLEATION COORDINATOR
!!!!! AUS_NUC_PEARLITE(NUC_N)			INSTANTANEOUS AUSTENITE NUCLEATION AT FERRITE/PEARLITE INTERFACES
!!!!! AUS_NUC_FERRITE(NUC_RATE0)		CONTINOUS NUCLEATION OF AUSTENITE AT FERRITE GRAIN BOUNDARIES
!!!!! AUS_NUC_FERRITE2(NUC_RATE0)		PREVIOUS VERSION
!!!!! SEEDING2(RAD,II,JJ,KK,PHASE)		PUT A NUCLEUS IN THE REGION CENTERED OF POINT(II,JJ,KK), RADIUS:RAD, NUCLEI PHASE:PHASE
!!!!! SEEDING(RAD,II,JJ,KK,PHASE1,PHASE2)	SIMILAR TO SEEDING2,BUT AVOID OCCUPING THE POINTS CORRESPONDING TO PHASE2

  Use Util_func
  Implicit None

 !!! TO SHIELD GRID POINTS AROUND A NELEUS, I.E. EXCLUDE THOSE POINTS FOR NEXT NUCLEATION EVENT
  Logical,Allocatable::Shield(:,:,:)
  


 Contains

Subroutine Nucleation_rex_alpha(Nuc_n)
  Integer Nuc_n
  If(Nuc_n==0)Return
  Write(*,*)"Nucleation Of New Ferrite Starts......."
  Call Random_number(Flag)
  Where(Ind_ph1(1,:,:,:)<=Stoich_n)Flag=Naf
  If(Nucleation_type==0)Then
        Call Homo(Nuc_n)
    Elseif(Nucleation_type==1)Then
        Call Critical(Nuc_n)
    Else
        Write(*,*)"Input Error For Nucleation Type"
        Stop
  Endif
  Write(*,*)"Ferrite Nucleation Of Rex"
  

End Subroutine


 Subroutine Homo(Nuc_n)
   Integer N,Xyz(3),I,J,K,Nuc_n
   
   Do N=1,Nuc_n
      Xyz=Maxloc(Flag)
      Call Seeding(Radius,Xyz(1),Xyz(2),Xyz(3),1,3)
   

   Enddo

 End Subroutine
 
 Subroutine Critical(Nuc_n)
 Integer N,Xyz(3),I,J,K,Nuc_n
 Write(*,*)"Critical_n Nucleating..."
 Theta_star=Theta_star*1.E-12
 !$omp Parallel Do &
 !$omp Shared(Ind_ph1,Ind_ph2,Val_ph1,Val_ph2,Flag) &
 !$omp Private(I,J,K)
 Do J=1,Iy
    Do K=1,Iz
         Do I=1,Ix
             
             If(Rex_energy(Ind_ph1(1,I,J,K))<Theta_star.Or.Any(Phase_index(Ind_ph1(:,I,J,K))==3))Flag(I,J,K)=Naf
         Enddo
    Enddo
 Enddo 
  !$omp End Parallel Do
  Do N=1,Nuc_n
     Xyz=Maxloc(Flag,Flag>0.)
     Call Seeding(Radius,Xyz(1),Xyz(2),Xyz(3),1,3)
!      Forall(I=1:Ix,J=1:Iy,K=1:Iz,&
!     Dis(I,J,K,(Xyz(1)-1)*Width,(Xyz(2)-1)*Width,(Xyz(3)-1)*Width)<=Max(2*Radius,2*Eta))Flag(I,J,K)=Naf
   Enddo

 End Subroutine

!!!!!!!!!!!!!***************************************
!!!!!!!!!!!!!***************************************

Subroutine Nucleation_austenite()

 
!!!!!Site-Saturated Nucleation In Pearlite-Ferrite Interface   
!!!!! Also Diffusion Starts To Be Considered    
        If(Flag_aus_nuc_pearlite.And.Temperature>Gamma_nuc_t)Then
                       Radius=Width*I_intf
                      ! Call Con_ini  !!!Initialize The Carbon Concentration Before Nucleation

                       Call Aus_nuc_pearlite(Gamma_nuc_n)
                       Flag_aus_nuc_pearlite=.False.
                       Flag_diffusion=.True.
                      

!                        Write(*,*)'Boundary Sites Number:',Ferrite_boundary_area0
        Endif

          If(Flag_aus_nuc_ferrite.And.Temperature>Aus_nuc_ferrite_t) Call Aus_nuc_ferrite2( Aus_nuc_ferrite_rate)   !!!! Continuous Nucleation At Ferrite Boundaries
        
        
End Subroutine



Subroutine Aus_nuc_pearlite(Nuc_n)
 Integer Substrate,Ordter
 Integer Nuc_n,I,J,K,N
 Integer Xyz(3),Ngb(3,0:Nngb),Grain_index(0:Nngb)
  Call Random_number(Flag)
!  !$omp Parallel Do&
!  !$omp Shared(Ind_ph1,Ind_ph2,Val_ph1,Val_ph2,Flag) &
!  !$omp Private(I,J,K)
   Do J=1,Iy
        Do K=1,Iz
           Do I=1,Ix

!                  If(Phase_index(Ind_ph1(1,I,J,K))/=3.Or.Ind_ph1(2,I,J,K)<Phasemin) Flag(I,J,K)=Naf



               If(Ind_ph1(1,I,J,K)>Stoich_n)Then
                   Flag(I,J,K)=Naf
                   
               Else
                 Call Get_ngb1(I,J,K,Ngb)
                 Do N=0,Nngb
                   Grain_index(N)=Maxloc(Val_phase(:,Ngb(1,N),Ngb(2,N),Ngb(3,N)),1)
                 Enddo
!                 If(All(Grain_index>Stoich_n))Flag(I,J,K)=Naf !!!Only Pearlite-Ferrite Interface On Ferrite Side
                 If(All(Grain_index==3))Flag(I,J,K)=Naf !!!Both Pearlite-Pearlite And Pearlite-Ferrite Interface
               Endif
           Enddo
       Enddo
   Enddo
!  !$omp End Parallel Do

 Do N=1,Nuc_n
    If(All(Flag<0.))Exit
    Xyz=Maxloc(Flag)
    
    Call Seeding(3.*Width,Xyz(1),Xyz(2),Xyz(3),2,1)
    
     
    
!     Forall(I=1:Ix,J=1:Iy,K=1:Iz,&
!     Dis(I,J,K,(Xyz(1)-1)*Width,(Xyz(2)-1)*Width,(Xyz(3)-1)*Width)<=Max(2*Radius,2*Eta))Flag(I,J,K)=Naf


 Enddo
 Write(*,*)N-1,"Austenite Nucleated In Pearlite Instantaneously"
 
End Subroutine

Subroutine Aus_nuc_ferrite(Nuc_rate0)
 
 Real Nuc_rate0,Nuc_rate,Randomnum
 Real,Save::Cumu_num=0.
 Integer,Save::Aus_num_fer=0
! Logical,Save::Flag_aus_nuc_ferrite=.True.
 
 Integer Ferrite_boundary_area,Num,I,J,K,N,Xyz(3)
  
If_potential_site=.True.

If(.Not.Flag_aus_nuc_ferrite)Return
  
  Ferrite_boundary_area=Count( All(Phase_index(Ind_ph1(1:2,:,:,:))==1,1))!!!Potential Area
  Nuc_rate=Nuc_rate0
 !Nuc_rate=Nuc_rate0*Exp(3.02e4*(1./Aus_nuc_ferrite_temp-1./Temperature))*Temperature/Aus_nuc_ferrite_temp
  Cumu_num=Cumu_num+Nuc_rate!*Ferrite_boundary_area /Ferrite_boundary_area0*Del_t
  Num=Floor(Cumu_num)!!! Nuclei Number
  
  If (Num==0)Return  !!!No Nuclei Is Introduced
 
!  !$omp Parallel Do&
!  !$omp Shared(Ind_ph1,Ind_ph2,Val_ph1,Val_ph2,Flag) &
!  !$omp Private(I,J,K,N)
  Do J=1,Iy
     Do K=1,Iz
        Do I=1,Ix
            
                  N=Count(Val_ph1(:,I,J,K)>Phasemin)
                  If( If_potential_site(I,J,K).And. N>1 .And.(All(Phase_index(Ind_ph1(1:N,I,J,K))==1)).And.(All(Ind_ph1(1:N,I,J,K)<=Npp0).Or.All(Ind_ph1(1:N,I,J,K)>Npp0)) )Then
                            Call Random_number(Randomnum)
                            Flag(I,J,K)=N*(Randomnum+Count(Ind_ph1(1:N,I,J,K)<=Npp0)/6.)
                  Else
                            Flag(I,J,K)=Naf
                  Endif
            
        Enddo
     Enddo
  Enddo
!  !$omp End Parallel Do
              
     If(All(Flag<0.))Then
          Flag_aus_nuc_ferrite=.False.
          Deallocate(If_potential_site)
          Write(*,*)"Austenite Nucleation At Ferrite Interface Is Done"
          Return
     Endif



  Do N=1,Num
    
     Xyz=Maxloc(Flag,If_potential_site)
     Aus_num_fer=Aus_num_fer+1
     
     Call Seeding2(Radius,Xyz(1),Xyz(2),Xyz(3),2)
!      Do K=1,Iz
!         Do J=1,Iy
!            Do I=1,Ix
!            
! !            If(Ind_ph1(1,I,J,K)==Tpp)Val_con1(1,I,J,K)=0.6  !!!!Initial Nucleus Concentration; By Defalt, It Is The Matrix Concentration Conforming To Mass Conservation
!            If(Dis(I,J,K,(Xyz(1)-1)*Width,(Xyz(2)-1)*Width,(Xyz(3)-1)*Width)<Max(2*Radius,2*Eta)) If_potential_site(I,J,K)=.False.
!            Enddo
!         Enddo
!      Enddo
   Enddo
  Cumu_num=Cumu_num-N+1
  Flag_aus_nuc_ferrite=.False.
End Subroutine

Subroutine Aus_nuc_ferrite2(Nuc_rate0)
 
 Real Nuc_rate0,Phi,Dg,Conc_fer
 Real,Save::Cumu_num=0.
 Integer Ferrite_boundary_area,Xyz(3),I,J,K,Nonzero,Iix,Iiy,Iiz,Ii,Jj,Kk,Shield_range
 Logical Nuc_flag
 
                        If(.Not.Allocated(Shield))Then
                                 Allocate(Shield(Ix,Iy,Iz))
                                 Shield=.True.
                       
                       Endif
!  If(.Not.Flag_aus_nuc_ferrite)Return

 Nuc_rate0=28.
  Conc_fer=Sum( Val_con1(1,:,:,:),Val_phase(1,:,:,:)>1.-Phasemin)/Sum(Val_phase(1,:,:,:),Val_phase(1,:,:,:)>1.-Phasemin)
  Dg=Ratio_k_para(1,2)*Partition_k_para(2,1)*(Conc_fer-Ref_c_para(1,2))+2.0e-12*Count(Ind_ph1(1,:,:,:)<=Npp0.And.Phase_index(Ind_ph1(1,:,:,:))==1)/Count(Phase_index(Ind_ph1(1,:,:,:))==1)
 If(Dg<0.)Return
  Phi=6.02e23*4.E-42
!  If_potential_site=All(Phase_index(Ind_ph1(1:2,:,:,:))==1).And.(Val_phase(2,:,:,:)<Phasemin).And.(All(Ind_ph1(:,:,:,:))
    Do K=1,Iz
      Do J=1,Iy
         Do I=1,Ix
            If (Shield(I,J,K).And.Val_ph1(2,I,J,K)>Phasemin.And.Val_phase(2,I,J,K)<Phasemin)Then  !!!! Interface Without Austenite
                 Nonzero=Count(Ind_ph1(:,I,J,K)>0)
                 If_potential_site(I,J,K)=All(Ind_ph1(1:Nonzero,I,J,K)<=Npp0).Or.All(Ind_ph1(1:Nonzero,I,J,K)>Npp0)
            Else
                 If_potential_site(I,J,K)=.False.
                 
            Endif
         Enddo
      Enddo
   Enddo
  Ferrite_boundary_area=Count(If_potential_site)
  If(Ferrite_boundary_area==0)Flag_aus_nuc_ferrite=.False.
 
!   If(Dg<Drag0(1))Return
 
  If(Dg>0.)Cumu_num=Cumu_num+Width**(Dms-1)*Ferrite_boundary_area*Exp(-(250000.+Phi/(Dg*Dg))/(8.314*Temperature) +Nuc_rate0)*Del_t
  
  
!    If  (Abs(Time-Write_t*Nint(Time/Write_t))<=Del_t/2.)Write(*,*)Cumu_num
Do While (.True.)
  If(Cumu_num<1.0)Return
!   If(Cumu_num>=2.0)Stop "Error: Nucleation Is Delayed"
 
  Call Random_number(Flag)
!   Where(Ind_ph1(1,:,:,:)<=Npp0)Flag=Flag+0.5!+Val_ph1(3,:,:,:)*10.
  
  Where(.Not.Shield)Flag=Naf
  Shield_range=Nint(1.5/Width)!2*I_intf
 
 
 !!!! Around Austenite There Is No Nucleation
  Do I=1,500
        Nuc_flag=.True.
        Xyz=Maxloc(Flag,If_potential_site)    
Outer:        Do Kk=-Min((Iz-1)/2,Shield_range),Min((Iz-1)/2,Shield_range)
           Iiz=Myngb(K,Iz,Kk)
           Do Jj=-Min((Iy-1)/2,Shield_range),Min((Iy-1)/2,Shield_range)
               Iiy=Myngb(J,Iy,Jj)
               Do Ii=-Min((Ix-1)/2,Shield_range),Min((Ix-1)/2,Shield_range)
                    Iix=Myngb(I,Ix,Ii)
                    If(Val_phase(2,Iix,Iiy,Iiz)>Phasemin) Then
                          Nuc_flag=.False.
                          Flag(Xyz(1),Xyz(2),Xyz(3))=-1.
                          Exit Outer
                    Endif                          
               Enddo
          Enddo
       Enddo Outer
       If(Nuc_flag)Exit
  Enddo
                    
                        
                    
  If(Nuc_flag) Call Seeding2(Eta,Xyz(1),Xyz(2),Xyz(3),2)
  Cumu_num=Cumu_num-1.0
 Enddo
End Subroutine

 Subroutine Seeding(Rad,Ii,Jj,Kk,Phase1,Phase2)!Radius:Nuclei Radius; X,Y,Z:Location; 

   Real Rad,Rex(Tpp),Totalc
   Integer I,J,K,Ii,Jj,Kk,Phase1,Phase2

   
   
   Tpp=Tpp+1

     !$omp Parallel Do &
    !$omp Shared(Ind_ph1,Ind_ph2,Val_ph1,Val_ph2,Val_con1,Val_con2,Val_phase,Flag) &
    !$omp Private(I,J,K) 
   Do K=1,Iz
      Do J=1,Iy
         Do I=1,Ix
           If(Dis(I,J,K,(Ii-1)*Width,(Jj-1)*Width,(Kk-1)*Width)<=Rad)Then  
                If(Any(Phase_index(Ind_ph1(:,I,J,K))==Phase2))Then
!                      If(Ind_ph1(1,I,J,K)>Npp0)  Coherent_grain(Tpp,1)=Ind_ph1(1,I,J,K)
                     Cycle 
                Endif
                Ind_ph1(1,I,J,K)=Tpp
                Ind_ph1(2:Npp,I,J,K)=0
                Val_ph1(1,I,J,K)=1.
                Val_ph1(2:Npp,I,J,K)=0.
                
                Val_con1(Phase1,I,J,K)=Dot_product(Val_con1(:,I,J,K),Val_phase(:,I,J,K))
                Val_phase(:,I,J,K)=0.
                Val_phase(Phase1,I,J,K)=1.
            Endif
            
            If(Dis(I,J,K,(Ii-1)*Width,(Jj-1)*Width,(Kk-1)*Width)<=Max(2*Radius,2*Eta))Then
               Flag(I,J,K)=Naf
               If_interface(I,J,K)=.True.
            Endif
          Enddo
        Enddo
     Enddo
     !$omp End Parallel Do
  
 
    
End Subroutine

Subroutine Seeding2(Rad,Ii,Jj,Kk,Phase)
!!!!Nucleation Within Boundary, Non-Spherical Shape, K-S Or N-W Relationship
   Real Rad,Rex(Tpp)
   Integer I,J,K,Ii,Jj,Kk,Index_vec(Npp),Phase
   integer Or_num !!Judge How Many Neighbor Grains Have K-S Or With Austenite: 50% With One And 50% With Two--Triple Junctions And Faces
   
   
   Index_vec=Ind_ph1(:,Ii,Jj,Kk)
   Tpp=Tpp+1
   Or_num=Min(2,Count(Phase_index(Index_vec)==1)-1 )
   Do I=1,Or_num
      If(Maxval(Index_vec)>0) Then
          Coherent_grain(Tpp,I)=Maxval(Index_vec)
          Index_vec(Maxloc(Index_vec))=0
      Endif
   Enddo

  
   Do K=1,Iz
      Do J=1,Iy
         Do I=1,Ix
           If(Dis(I,J,K,(Ii-1)*Width,(Jj-1)*Width,(Kk-1)*Width)<=Rad .And.(Ind_ph1(2,I,J,K)>0.Or.All(Coherent_grain(Tpp,:)/=Ind_ph1(1,I,J,K)))) Then  !!!.Or.All(Coherent_grain(Tpp,:)/=Ind_ph1(1,I,J,K))
                
                Ind_ph1(1,I,J,K)=Tpp
                Ind_ph1(2:Npp,I,J,K)=0
                Val_ph1(1,I,J,K)=1.
                Val_ph1(2:Npp,I,J,K)=0.

                Val_con1(Phase,I,J,K)=Dot_product(Val_con1(:,I,J,K),Val_phase(:,I,J,K))!Partition_k(Phase,1)*Dot_product(Val_con1(:,I,J,K),Val_phase(:,I,J,K))
                Val_phase(:,I,J,K)=0.
                Val_phase(Phase,I,J,K)=1.   
                
            Endif
            
           If(Dis(I,J,K,(Ii-1)*Width,(Jj-1)*Width,(Kk-1)*Width)<Max(1.5,2.*Eta)) Then
              Shield(I,J,K)=.False.
              If_interface(I,J,K)=.True.
           Endif

          Enddo
        Enddo
     Enddo
  
!     If(Maxval(Index_vec)>Npp0) Then
!         Coherent_grain(Tpp,1)=Maxval(Index_vec)
!         Index_vec(Maxloc(Index_vec))=0
!     Endif
  
!     Call Random_number(Or_num)
!    If(Or_num<1./Heat_rate2.And.Maxval(Index_vec)>Npp0)Coherent_grain(Tpp,2)=Maxval(Index_vec)!!!
End Subroutine

End Module
