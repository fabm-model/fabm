! The following program provides an example implimentation of the routines.

   PROGRAM Run_CO2_Dynamics

   Real*8 :: Temp, Sal, Depth, DIC, pco2w, pco2a, TA, ph
   Real*8 :: cco2, carba,bicarb,carb,henry,om_cal,om_arg,TCO2,dcf
   Real*8 :: wnd, flux

! set up inputs
   Temp  = 5.0		!C
   Sal   = 35.0		!psu
   Depth = 0.0		!m
   DIC   = 2200.0	!mmol/m3
   Wnd   = 10.0		!m/s
   pCO2a = 390.0 	!uatm
!!!! note you need to define total alkalinity in subroutine co2_dynamics (see lines 71-86) !!!

   Call CO2_dynamics (Temp, Sal, Depth, DIC,pco2w,TA,ph,carba,bicarb,carb,henry,om_cal,om_arg,TCO2,dcf)

   Call Air_sea_exchange (Temp, Wnd, pCO2w, pCO2a, Henry, dcf, flux)

! write inputs to screen
   WRITE(*,*) " "
   WRITE(*,'(A28)') "    .........Inputs........."
   WRITE(*,'(A24,F10.3)') "    Temperature   (T) = ",Temp
   WRITE(*,'(A24,F10.3)') "    Salinity    (psu) = ",Sal
   WRITE(*,'(A24,F10.3)') "    Depth         (m) = ",Depth
   WRITE(*,'(A24,F10.3)') "    DIC     (mmol/m3) = ",DIC
   WRITE(*,'(A24,F10.3)') "    Wind speed  (m/s) = ",Wnd
   WRITE(*,'(A24,F10.3)') "    pCO2 atmos (uatm) = ",pCO2a
 
! write outputs to screen
   WRITE(*,*) " "
   WRITE(*,'(A32,F10.3)') "    ..........Outputs..........."
   WRITE(*,'(A27,F10.3)') "    pH                 (pH) = ",pH
   WRITE(*,'(A27,F10.3)') "    DIC           (umol/kg) = ",TCO2
   WRITE(*,'(A27,F10.3)') "    TA            (umol/kg) = ",TA
   WRITE(*,'(A27,F10.3)') "    pco2w            (uatm) = ",pco2w
   WRITE(*,'(A27,F10.3)') "    carbonic acid (mmol/m3) = ",carba
   WRITE(*,'(A27,F10.3)') "    bicarbonate   (mmol/m3) = ",bicarb
   WRITE(*,'(A27,F10.3)') "    carbonate     (mmol/m3) = ",carb
   WRITE(*,'(A27,F10.3)') "    Omega calcite       (~) = ",om_cal
   WRITE(*,'(A27,F10.3)') "    Omega aragonite     (~) = ",om_arg
   WRITE(*,'(A27,F10.3)') "    air sea flux(mmol/m2/d) = ",flux
   WRITE(*,*) " "


   STOP
   END PROGRAM Run_CO2_Dynamics
