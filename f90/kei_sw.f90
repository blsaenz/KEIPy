
MODULE kei_sw

    use kei_parameters
    use kei_common
    use macmods_kinds_mod
    use macmods_param_mod
    use macmods_mod

    use kei_ecocommon, only :  po4_ind,  & ! dissolved inorganic phosphate
      no3_ind     ,  & ! dissolved inorganic nitrate
      sio3_ind    ,  & ! dissolved inorganic silicate
      nh4_ind     ,  & ! dissolved ammonia
      fe_ind      ,  & ! dissolved inorganic iron
      o2_ind      ,  & ! dissolved oxygen
      dic_ind     ,  & ! dissolved inorganic carbon
      alk_ind     ,  & ! alkalinity
      doc_ind     ,  & ! dissolved organic carbon
      spC_ind     ,  & ! small phytoplankton carbon
      spChl_ind   ,  & ! small phytoplankton chlorophyll
      spCaCO3_ind ,  & ! small phytoplankton caco3
      diatC_ind   ,  & ! diatom carbon
      diatChl_ind ,  & ! diatom chlorophyll
      zooC_ind    ,  & ! zooplankton carbon
      spFe_ind    ,  & ! small phytoplankton iron
      diatSi_ind  ,  & ! diatom silicon
      diatFe_ind  ,  & ! diatom iron
      diazC_ind   ,  & ! diazotroph carbon
      diazChl_ind ,  & ! diazotroph Chlorophyll
      diazFe_ind  ,  & ! diazotroph iron
      don_ind     ,  & ! dissolved organic nitrogen
      dofe_ind    ,  & ! dissolved organic iron
      dop_ind          ! dissolved organic phosphorus


    IMPLICIT NONE

    PUBLIC

    integer, parameter :: n_sw_output = 36

    ! Keep these around, so we can write them out
    REAL(KIND=r8), DIMENSION(IMT,JMT), SAVE :: sw_lat,sw_lon,sw_sst, &
      sw_chl,sw_par,sw_swh,sw_mwp,sw_cmag,sw_no3_,sw_nh4_,sw_po4_,sw_fe_ !underscores bc of f2py link conflict

    INTEGER(kind=i4), SAVE :: sw_i ! current index of seaweed into the 1D vertical grid

    REAL(KIND=r8), SAVE :: sw_absorp


  CONTAINS

    SUBROUTINE sw_init()

        CALL macmods_init(.false.) ! read the param file
        sw_i = -1

    END SUBROUTINE sw_init

    function sw_depth_index(doy) result(i)
        ! using dm, the grid cell boundary depths, to find the best location of seaweed in the grid
        ! i should be the 1-based index of the cell, because dm is zero-based

        real(kind=r8), intent(in) :: doy
        real(kind=r8) :: depth
        integer :: i

        call macmods_get_depth(doy,depth)
        do i=0,NZ
            if (depth < dm(i)) then
                exit
            endif
        enddo

    end function sw_depth_index


    SUBROUTINE kei_sw_step(U,X,doy,dtday,nt,par_phyto,swh_in,mwp_in,cmag_in)

        real :: U(NZP1,NVEL), X(NZP1,NSCLR) ! nvel=2 typically
        real(kind=r8), intent(in) :: doy,swh_in,mwp_in,cmag_in,dtday
        integer(kind=i4), intent(in) :: nt
        REAL(KIND=r8), DIMENSION(IMT,JMT) :: f_values
        real(KIND=r8), intent(in), dimension(NZ) :: par_phyto ! shortwave irradiance, transmitted through ice, attenudated by phytos

        integer(kind=i4) :: Xoffset = 2  ! ecosys tracers are offset by 2 ... sigh...
        real(KIND=r8), dimension(1,1) :: tendency


        ! This is important. zero stuff, since we will do accounting and store/write the block
        ! at every step.
        ! -------------------------------------------------------
        CALL macmods_zero_outputs()

        ! find (and store) current depth/grid index of seaweed
        ! for very long kelps, this could be modified to maybe use the max of
        ! a range of grid cells, or the average of grid cells?
        sw_i = sw_depth_index(doy)
        if (sw_i > NZ) then
            print *,'MACMODS seaweed depth is lower than grid (sw_i,NZ)!',sw_i,NZ
            CALL exit(-1)
        endif

        sw_lat = dlat
        CALL macmods_set_forcing('lat',sw_lat) ! degrees latitude
        sw_lon = dlon
        CALL macmods_set_forcing('lon',sw_lon) ! degrees longitude

        !f_values = X(sw_i,2) + Sref ! don't forget that salinity is not salinity w/out Sref...
        !CALL macmods_set_forcing('sss',f_values)
        sw_sst = X(sw_i,1)
        CALL macmods_set_forcing('sst',sw_sst)  ! sea (surface) temperature - legacy name, will be at whatever depth assigned

        !sw_chl = X(sw_i, Xoffset + diatChl_ind) &
        !           + X(sw_i, Xoffset + diazChl_ind) &
        !           + X(sw_i, Xoffset + spChl_ind)
        ! instead of above, we are using direct par_phyto from ecosys, for consistency
        sw_chl = -1.0 ! set to negative, and macmods_calc will assume par is at depth
        CALL macmods_set_forcing('chl',sw_chl)  ! phytoplankton chloropyll-a, used for attenuating light at depth, if enabled
        sw_par = par_phyto(sw_i)
        CALL macmods_set_forcing('par',sw_par)  ! photosyntheticaly active radiation (W)

        ! these are going to have to come from forcing or invented by KEI based on winds over time?
        sw_swh = swh_in
        CALL macmods_set_forcing('swh',sw_swh)  ! swell height (m)
        sw_mwp = mwp_in
        CALL macmods_set_forcing('mwp',sw_mwp)  ! mean wave period (s)
        sw_cmag = cmag_in
        CALL macmods_set_forcing('cmag',sw_cmag)  ! current magnitude (m/s)

        ! don't use nflux
        !nflux, &  ! maximum vertical NO3-flux at 100m, if used

        ! nutrients from ecosys tracers
        sw_no3_ = X(sw_i, Xoffset + no3_ind)
        CALL macmods_set_forcing('no3',sw_no3_)  ! ambient Nitrate (mmol/m3)
        sw_nh4_ = X(sw_i, Xoffset + nh4_ind)
        CALL macmods_set_forcing('nh4',sw_nh4_)  ! ambient ammonia (mmol/m3)
        sw_po4_ = X(sw_i, Xoffset + po4_ind)
        CALL macmods_set_forcing('po4',sw_po4_)  ! ambient Nitrate (mmol/m3)
        sw_fe_ = X(sw_i, Xoffset + fe_ind)
        CALL macmods_set_forcing('fe',sw_fe_)   ! ambient bioavailable iron (mmol/m3)

        ! COMPUTE
        CALL macmods_tendency_compute(dtday,doy)

        ! a guess at absorption from sw biomass?
        ! Currently heat and PAR are calculated differently ... start by
        ! passing this into ecosystem model, but this and chl-a should be incorporated
        ! info sw_frac in kei_ocn soon...  Maybe some fraction should be reflected/disappeared?
        sw_absorp = macmods_get_output_1(mo_d_B,1,1) / 50._r8

        ! APPLY TENDENCIES
        X(sw_i, Xoffset + dic_ind) = X(sw_i, Xoffset + dic_ind) + macmods_get_output_1(mo_d_DIC,1,1)
        X(sw_i, Xoffset + nh4_ind) = X(sw_i, Xoffset + nh4_ind) + macmods_get_output_1(mo_d_NH4,1,1)
        X(sw_i, Xoffset + no3_ind) = X(sw_i, Xoffset + no3_ind) + macmods_get_output_1(mo_d_NO3,1,1)
        X(sw_i, Xoffset + po4_ind) = X(sw_i, Xoffset + po4_ind) + macmods_get_output_1(mo_d_PO4,1,1)
        X(sw_i, Xoffset + fe_ind) = X(sw_i, Xoffset + fe_ind) + macmods_get_output_1(mo_d_Fe,1,1)
        !X(sw_i, Xoffset + poc_ind) = X(sw_i, Xoffset + poc_ind) + tend_POC = c0 ! need POC input to ecosys
        X(sw_i, Xoffset + doc_ind) = X(sw_i, Xoffset + doc_ind) + macmods_get_output_1(mo_d_DOC,1,1)
        X(sw_i, Xoffset + don_ind) = X(sw_i, Xoffset + don_ind) + macmods_get_output_1(mo_d_DON,1,1)
        X(sw_i, Xoffset + dofe_ind) = X(sw_i, Xoffset + dofe_ind) + macmods_get_output_1(mo_d_DOFe,1,1)
        X(sw_i, Xoffset + o2_ind) = X(sw_i, Xoffset + o2_ind) + macmods_get_output_1(mo_d_O2,1,1)

    END SUBROUTINE kei_sw_step


    ! copy all the tracers into a single array here.  This facilitates passing
    ! back to python in the the kei_link module.
    SUBROUTINE sw_get_outputs(output_array)

      real(kind=r8), dimension(n_sw_output), intent(inout) :: output_array

      ! moved to block copy using get-subroutines
      CALL macmods_get_output_block_1(output_array,1,1)

    END SUBROUTINE sw_get_outputs


END MODULE kei_sw