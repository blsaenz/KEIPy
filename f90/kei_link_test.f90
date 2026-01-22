PROGRAM KEI_link_test


  USE kei_parameters
  USE kei_common
  USE kei_icecommon
  USE kei_ice
  USE kei_ecocommon
  USE kei_ocncommon
  USE kei_eco
  USE kei_hacks
  USE kei_sw
  USE link

  INTEGER :: iii

  REAL :: &
    U_local(NZ,NVEL),  &    ! momentum
    X_local(NZ,NSCLR),  &   ! tracers
    Fcomp(1000,19), &   ! forcing
    dm_local(NZP1), &
    hm_local(NZP1), &
    zm_local(NZP1)

  DOUBLE PRECISION :: doy, sw_out(n_outputs)


  open(12, file="../test_data/kf_200_100_2000_savetxt.txt")
  read(12,*) Fcomp(:,1:16)
  close(12)
  Fcomp(:,msl_f_ind) = Fcomp(:,msl_f_ind)*0.01  ! Pa -> mbar -- need to sort this out in python

  ! seaweed forcing vars were not originally included in that output file above
  Fcomp(:,swh_f_ind) = 0.25   ! m
  Fcomp(:,mwp_f_ind) = 10.0   ! s
  Fcomp(:,cmag_f_ind) = 0.25  ! m/s


  open(12, file="../test_data/U_savetxt.txt")
  read(12,*) U_local
  close(12)

  open(12, file="../test_data/X_savetxt.txt")
  read(12,*) X_local
  close(12)

  open(12, file="../test_data/dm_savetxt.txt")
  read(12,*) dm_local
  close(12)

  open(12, file="../test_data/hm_savetxt.txt")
  read(12,*) hm_local
  close(12)

  open(12, file="../test_data/zm_savetxt.txt")
  read(12,*) zm_local
  close(12)

  CALL set_param_int("nend",999)

  CALL kei_param_init()

  CALL set_grid(dm_local,hm_local,zm_local)

  CALL set_tracers(U_local,X_local)

  CALL kei_compute_init()

  doy = 60.D0
  DO nt = 1,999

    print *,'Forcing:'
    DO iii=1,16
        print *,Fcomp(nt,iii)
    ENDDO
    CALL set_forcing(Fcomp(nt,:))
    CALL kei_compute_step(nt,doy)

    CALL get_tracers(U_local,X_local)
    print *,nt
    DO iii=1,60
        print *,X_local(iii,1),X_local(iii,2)+Sref ! temperature,salinity
    ENDDO
    CALL get_sw_data(sw_out)
    print *,'Seaweed Biomass:',sw_out(1)
    print *,'Seaweed QN:',sw_out(2)
    print *,'Seaweed d_NO3:',sw_out(16)
    print *,'Seaweed d_DOC:',sw_out(22)

    doy = doy + dtday
    if (doy >= 367.D0) then
      doy = doy - 367.D0
    endif

    ! write some functions to get macmods data and print it?

    ! Do a nutrient mass conservation test?

  ENDDO


END PROGRAM KEI_link_test