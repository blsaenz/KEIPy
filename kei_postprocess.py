import os,h5py,netCDF4,time,datetime
import numpy as np
import PyCO2SYS as pyco2
import matplotlib.pyplot as plt
from matplotlib.dates import num2date,date2num
import matplotlib.dates as mdates
import cmocean
import xarray as xr
import pandas as pd
import scipy.io as spio

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def kei_dn(k,start_yr):
    td = datetime.timedelta(hours=1)
    dt = datetime.datetime(start_yr,1,1) + datetime.timedelta(days=k['day'][0])
    dtt = [dt+td*i for i in range(len(k['day']))]
    return date2num(dtt)

def kei_read_legacy(filename):

    matfilename = filename[:-2] + "mat"
    sal_ref = 34.6;

    if os.path.exists(matfilename):
        k = loadmat(matfilename)['k']
        return k

    else:
        print('Reading KEI netcdf ...')
        k = {}
        stime=time.time()
        # with h5py.File(matfilename,'r') as f:
        #     for gname in f.keys():
        #         k[gname] = f[gname][...]
        ds = netCDF4.Dataset(filename,"r")
        for v in get_read_vars():
            print('reading',v,'...')
            k[v] = ds.variables[v][...]
        print('read: %s seconds'%(time.time()-stime))

        nt,nz = np.shape(k['S'])

        # types of carbonate ssytem params: 1=ALK [micromol/kg], 2=DIC [micromol/kg], 3=pH, 4=pCO2 [microATM], 6=carbonate ion [micromol/kg], 7=bicarbonate iom [micromol/kg]
        carb = pyco2.sys(k['ALK'], k['DIC'], 1, 2, temperature=k['T'], salinity=k['S'], pressure=-1*k['zgrid'],
                         total_phosphate=k['PO4'], total_silicate=k['SiO3'])

        k['pH'] = carb['pH_total']

        # pull interesting fluxes

        k['deep_ice_production'] = k['focn'][:,8]  # this is included in other terms below
        k['surface_ice_production'] = -k['focn'][:,7]*334000  #(convert kg to W)
        k['total_ice_production'] = k['fio'][:,7]*334000   # < 0 (W)
        # behavior of ice on ocn in response to everything - usually negative
        # b/c melting.  Positive if melted completely and didn't use up ML heat
        # actual melt flux
        k['ice_to_ocn_flux'] = k['fio'][:,4]  # (total, already scaled by ice fraction)
        # effectively the mixed layer melt potential
        k['ocn_to_ice_flux'] = k['fio'][:.8] # not scaled by ice fraction

        # freeze potential
        k['ocn_freeze_potential'] = k['fio'][:,7]*334000 # kg -> Watts

        #k['atm_flux_to_ice_surface']=double(ncread(filename,'atm_flux_to_ice_surface'));  % already scaled by fraction ice
        # k.atm_flux_to_ocn_surface=double(ncread(filename,'atm_flux_to_ocn_surface'));  % already scaled by fraction ocean
        #
        # k.total_ice_freeze=double(ncread(filename,'total_ice_freeze'));
        # k.total_ice_melt=double(ncread(filename,'total_ice_melt'));
        # k.frazil_ice_volume=double(ncread(filename,'frazil_ice_volume'));
        # k.congelation_ice_volume=double(ncread(filename,'congelation_ice_volume'));
        # k.snow_ice_volume=double(ncread(filename,'snow_ice_volume'));
        # k.snow_precip_mass=double(ncread(filename,'snow_precip_mass'));
        # k.ice_ocean_bottom_flux_potential=double(ncread(filename,'ice_ocean_bottom_flux_potential'));
        # k.ice_ocean_bottom_flux=double(ncread(filename,'ice_ocean_bottom_flux'));
        #

        # pressure = sw_pres(-k['zgrid'],-67)

        # calc extra bulk properties
        # for i=1:n_steps
        # for in range(nt):
        #     [~,k1] = min(k.T(:,i));
        #     k1=max([70,k1]);...
        #     k.deep_pyc(i) = mixedl(k.T(k1:end,i),k.T(k1,i),-1.*(k1-1:401-0.5));
        #     dp_i = round(abs(k.deep_pyc(i)));  % only works cause 1m grid!
        #     k.deep_pyc_end(i) = upside_down_mixedl(k.T(:,i),k.T(drange-1,i),-1.*(0:drange-0.5));
        #     k.deep_pyc_slope(i) = (k.T(dp_i+20,i) - k.T(dp_i,i))/20;
        #     k.deep_pyc_t(i) = mean(k.T(dp_i:dp_i+20,i));
        #     k.cp(:,i) = CPSW(k.S(:,i),k.T(:,i),-k.zgrid);
        #     k.sigma(:,i) = sw_dens(k.S(:,i),k.T(:,i),pressure);
        #     Tf = sw_fp(k.S(:,i),pressure);
        #     ml_idx = floor(-k.zml(i));
        #     k.heat(:,i) = k.cp(:,i).*k.sigma(:,i).*(k.T(:,i)-Tf);
        #     k.heat_60(i) = sum(k.heat(1:60,i));
        #     k.heat_ml(i) = sum(k.heat(1:ml_idx,i));
        #     k.heat_gr_60(i) = sum(k.heat(61:drange-1,i));
        #     sal_ref_i = find(k.S(60:end,i) > sal_ref,1) + 59;
        #     if isempty(sal_ref_i)
        #         sal_ref_i = drange-1;
        #     end
        #     kc = (sal_ref - k.S(1:sal_ref_i,i))./sal_ref;
        #     k.fresh(i) = sum(kc);
        #     kc = (sal_ref - k.S(1:ml_idx,i))./sal_ref;
        #     k.fresh_ml(i) = sum(kc);
        #
        #     k.cp_prev(:,i) = CPSW(k.Sprev(:,i),k.Tprev(:,i),-k.zgrid);
        #     sigma(:,i) = sw_dens(k.Sprev(:,i),k.Tprev(:,i),pressure);
        #     Tf = sw_fp(k.Sprev(:,i),pressure);
        #     k.heat_prev(:,i) = k.cp_prev(:,i).*sigma(:,i).*(k.Tprev(:,i)-Tf);
        #
        #     if (i > 1)
        #         %boundarylayer = round(k.hmx(i));
        #         %boundarylayer_1 = round(k.hmx(i-1));
        #         %k.pbl_flux(i) = sum(k.heat(1:boundarylayer,i) - k.heat_prev(1:boundarylayer,i))/(60*60); %dheat/timestep = flux
        #         k.pbl_flux(i) = (sum(k.heat(1:60,i)) - sum(k.heat(1:60,i-1)))/(60*60); %dheat/timestep = flux
        #         k.pyc_flux(i) = k.pbl_flux(i) - k.atm_flux_to_ocn_surface(i) - k.ice_to_ocn_flux(i) + k.total_ice_production(i);
        #     end
        # end
        #
        #
        # days = round(length(k.atm_flux_to_ice_surface)/24)-1;
        # for i=1:days
        #     k.fice_daily(i) = mean(k.fice(i*24-23:i*24));
        #     k.fice_mean24(i*24-23:i*24) = k.fice_daily(i);
        #     k.surface_ice_production_daily(i) = mean(k.surface_ice_production(i*24-23:i*24));
        #     k.surface_ice_production_mean24(i*24-23:i*24) = k.surface_ice_production_daily(i);
        #     k.total_ice_production_daily(i) = mean(k.total_ice_production(i*24-23:i*24));
        #     k.total_ice_production_mean24(i*24-23:i*24) = k.total_ice_production_daily(i);
        #     k.ice_to_ocn_flux_daily(i) = mean(k.ice_to_ocn_flux(i*24-23:i*24));
        #     k.ice_to_ocn_flux_mean24(i*24-23:i*24) = k.ice_to_ocn_flux_daily(i);
        #     k.ocn_to_ice_flux_daily(i) = mean(k.ocn_to_ice_flux(i*24-23:i*24));
        #     k.ocn_to_ice_flux_mean24(i*24-23:i*24) = k.ocn_to_ice_flux_daily(i);
        #     k.atm_flux_to_ice_surface_daily(i) = mean(k.atm_flux_to_ice_surface(i*24-23:i*24));
        #     k.atm_flux_to_ice_surface_mean24(i*24-23:i*24) = k.atm_flux_to_ice_surface_daily(i);
        #     k.atm_flux_to_ocn_surface_daily(i) = mean(k.atm_flux_to_ocn_surface(i*24-23:i*24));
        #     k.atm_flux_to_ocn_surface_mean24(i*24-23:i*24) = k.atm_flux_to_ocn_surface_daily(i);
        #
        #     k.ocn_freeze_potenial_daily(i) = mean(k.ocn_freeze_potenial(i*24-23:i*24));
        #     k.ocn_freeze_potenial_mean24(i*24-23:i*24) = k.ocn_freeze_potenial_daily(i);
        #
        #     k.ice_ocean_bottom_flux_daily(i)=mean(k.ice_ocean_bottom_flux(i*24-23:i*24));
        #
        #     k.hmx_daily(i) = mean(k.hmx(i*24-23:i*24));
        #     k.zml_daily(i) = mean(k.zml(i*24-23:i*24));
        #     k.hi_daily(i) = mean(k.hi(i*24-23:i*24));
        #     k.hsn_daily(i) = mean(k.hsn(i*24-23:i*24));
        #
        #     k.deep_pyc_daily(i) = mean(k.deep_pyc(i*24-23:i*24));
        #     k.deep_pyc_end_daily(i) = mean(k.deep_pyc_end(i*24-23:i*24));
        #     k.deep_pyc_slope_daily(i) = mean(k.deep_pyc_slope(i*24-23:i*24));
        #     k.deep_pyc_t_daily(i) = mean(k.deep_pyc_t(i*24-23:i*24));
        #
        #     k.fresh_daily(i) = mean(k.fresh(i*24-23:i*24));
        #     k.fresh_ml_daily(i) = mean(k.fresh_ml(i*24-23:i*24));
        #
        #     k.heat_60_daily(i) = mean(k.heat_60(i*24-23:i*24));
        #     k.heat_ml_daily(i) = mean(k.heat_ml(i*24-23:i*24));
        #     k.heat_gr_60_daily(i) = mean(k.heat_gr_60(i*24-23:i*24));
        #
        #     k.total_ice_freeze_daily(i) = sum(k.total_ice_freeze(i*24-23:i*24));
        #
        #     k.frazil_ice_volume_daily(i)=sum(k.frazil_ice_volume(i*24-23:i*24));
        #     k.congelation_ice_volume_daily(i)=sum(k.congelation_ice_volume(i*24-23:i*24));
        #     k.snow_ice_volume_daily(i)=sum(k.snow_ice_volume(i*24-23:i*24));
        #     k.total_ice_volume_daily(i)=k.frazil_ice_volume_daily(i) + k.congelation_ice_volume_daily(i) + k.snow_ice_volume_daily(i);
        #
        #     k.turb_daily_150(i)=sum(sum(k.wT(1:150,i*24-23:i*24)));
        #
        #
        #     k.pbl_flux_daily(i) = mean(k.pbl_flux(i*24-23:i*24));
        #     k.pyc_flux_daily(i) = mean(k.pyc_flux(i*24-23:i*24));
        #
        #     % find zero-degree depth
        #     k.T_daily(:,i) = mean(k.T(:,i*24-23:i*24),2);
        #     k.d0(i) = -1;
        #     j = 200;
        #     while (k.d0(i) == -1 ) && (j > 0)
        #         if k.T_daily(j,i) < 0
        #             k.d0(i) = j;
        #         end
        #         j=j-1;
        #     end
        #
        #     if j == 0
        #         %disp( 'zero degree water not found in kei_read_plotting!!!')
        #         %disp(sprintf('filename: %s, day: %i',filename,i))
        #         [~,k.d0(i)] = min(k.T_daily(50:end,i));
        #         k.d0(i) = k.d0(i)+49;
        #         %return
        #     end
        #
        # end
        #
        # k.d0 = k.d0*(-1);
        #
        # % save for later
        # save(matfilename,'k')
        #
        #
        # return




def get_read_vars():
    return [
    'zgrid',
    'd',
    'step',
    'time',
    'day',
    'hour',
    'wT',
    'wS',
    'T',
    'S',
    'hmx',
    'zml',
    'km',
    'ks',
    'fatm',
    'fao',
    'fai',
    'fio',
    'focn',
    'flx',
    'hi',
    'hsn',
    'fice',
    'hfsnow',
    'ns',
    'ni',
    'dzi',
    'dzs',
    'ti',
    'si',
    'ts',
    'qsi',
    'tas',
    'ta',
    'ks',
    'kt',
    'tot_prod',
    'diatC',
    'diatChl',
    'spC',
    'spChl',
    'Fe',
    'zooC',
    'NO3',
    'NH4',
    'PO4',
    'SiO3',
    'DIC',
    'ALK',
    'Tprev',
    'Sprev',
    'atm_flux_to_ice_surface',
    'atm_flux_to_ocn_surface',
    'total_ice_freeze',
    'total_ice_melt',
    'frazil_ice_volume',
    'congelation_ice_volume',
    'snow_ice_volume',
    'snow_precip_mass',
    'ice_ocean_bottom_flux_potential',
    'ice_ocean_bottom_flux',
    ]

def CPSW(S,T1,P0):
    # ******************************************************************
    # UNITS:
    #       PRESSURE        P0       DECIBARS
    #       TEMPERATURE     T        DEG CELSIUS (IPTS-68)
    #       SALINITY        S        (IPSS-78)
    #       SPECIFIC HEAT   CPSW     J/(KG DEG C)
    # ***
    # REF: MILLERO ET AL,1973,JGR,78,4499-4507
    #       MILLERO ET AL, UNESCO REPORT NO. 38 1981 PP. 99-188.
    # PRESSURE VARIATION FROM LEAST SQUARES POLYNOMIAL
    # DEVELOPED BY FOFONOFF 1980.
    # ***
    # CHECK VALUE: CPSW = 3849.500 J/(KG DEG. C) FOR S = 40 (IPSS-78),
    # T = 40 DEG C, P0= 10000 DECIBARS

    # check that temperature is above -2
    T = T1
    T[T < -2.0] = -2.0

    # SCALE PRESSURE TO BARS
    P=P0/10.0

    # SQRT SALINITY FOR FRACTIONAL TERMS
    SR = np.sqrt(np.abs(S))
    # SPECIFIC HEAT CP0 FOR P=0 (MILLERO ET AL ,UNESCO 1981)
    A = (-1.38385E-3*T+0.1072763)*T-7.643575
    B = (5.148E-5*T-4.07718E-3)*T+0.1770383
    C = (((2.093236E-5*T-2.654387E-3)*T+0.1412855)*T -3.720283)*T+4217.4
    CP0 = (B*SR + A)*S + C
    # CP1 PRESSURE AND TEMPERATURE TERMS FOR S = 0
    A = (((1.7168E-8*T+2.0357E-6)*T-3.13885E-4)*T+1.45747E-2)*T -0.49592
    B = (((2.2956E-11*T-4.0027E-9)*T+2.87533E-7)*T-1.08645E-5)*T +2.4931E-4
    C = ((6.136E-13*T-6.5637E-11)*T+2.6380E-9)*T-5.422E-8
    CP1 = ((C*P+B)*P+A)*P
    # CP2 PRESSURE AND TEMPERATURE TERMS FOR S > 0
    A = (((-2.9179E-10*T+2.5941E-8)*T+9.802E-7)*T-1.28315E-4)*T +4.9247E-3
    B = (3.122E-8*T-1.517E-6)*T-1.2331E-4
    A = (A+B*SR)*S
    B = ((1.8448E-11*T-2.3905E-9)*T+1.17054E-7)*T-2.9558E-6
    B = (B+9.971E-8*SR)*S
    C = (3.513E-13*T-1.7682E-11)*T+5.540E-10
    C = (C-1.4300E-12*T*SR)*S
    CP2 = ((C*P+B)*P+A)*P
    # SPECIFIC HEAT RETURN
    return CP0 + CP1 + CP2


def TfrzC(S,Db):
    """Freezing point of water in degrees C at salinity S in PSU
    and pressure Db in decibars.
    TfrzC = (-0.0575 +1.710523e-3 *sqrt(S) -2.154996e-4 *S) *S - 7.53e-4 *Db
    """
    return -0.054*S - 7.53e-4*Db

def wd(arr,IO,all_pos_neg='all'):
    """ find weekly mean, or not"""
    r = len(arr)%7
    if IO:
        this_arr = np.mean(arr[:-r].reshape(-1, 7), axis=1)
    else:
        this_arr = np.copy(arr)

    if all_pos_neg == 'pos':
        this_arr[this_arr < 0] = 0
        this_arr[this_arr==0] = np.nan
    elif all_pos_neg == 'neg':
        this_arr[this_arr > 0] = 0
        this_arr[this_arr==0] = np.nan

    return this_arr

def daily_fluxes(k,sday,eday):
    f = {}
    # atm -> pbl
    f['atm2pbl'] = np.copy(k['atm_flux_to_ocn_surface_daily'])
    f['atm2pbl'][f['atm2pbl']<0] = 0.0
    # pbl -> atm
    f['pbl2atm'] = np.copy(k['atm_flux_to_ocn_surface_daily']) * -1.0
    f['pbl2atm'][f['pbl2atm']<0] = 0.0
    # atm -> ice
    f['atm2ice'] = np.copy(k['atm_flux_to_ice_surface_daily'])
    f['atm2ice'][f['atm2ice'] < 0] = 0.0
    # ice -> atm
    f['ice2atm'] = np.copy(k['atm_flux_to_ice_surface_daily']) * -1.0
    f['ice2atm'][f['ice2atm'] < 0] = 0.0
    # ice -> pbl -- never happens, ice can't heat water
    # pbl -> ice
    f['pbl2ice'] = np.copy(k['ice_ocean_bottom_flux_daily'])
    f['pbl2ice'][k['fice_daily'] < 0.01] = 0.0
    # pyc -> pbl
    f['pyc2pbl'] = np.copy(k['pyc_flux_daily'])
    f['pyc2pbl'][f['pyc2pbl'] < 0] = 0.0
    # pbl -> pyc
    f['pbl2pyc'] = np.copy(k['pyc_flux_daily']) * -1.0
    f['pbl2pyc'][f['pbl2pyc'] < 0] = 0.0

    for k,v in f.items():
        f[k] = v[sday:eday]
    return f


def stackplot_plus_minus(axx,xvalues,yvalues,y_order):

    neg_tot = np.zeros(np.shape(xvalues))
    pos_tot = np.zeros(np.shape(xvalues))
    zorder = len(y_order) + 1
    p = []
    for y_name in y_order:
        y_pos = np.copy(yvalues[y_name])
        y_pos[y_pos<0] = 0
        y_neg = np.copy(yvalues[y_name])
        y_neg[y_neg>0] = 0

        pos_tot += y_pos
        neg_tot += y_neg

        p.append(axx.stackplot(xvalues,pos_tot+neg_tot,zorder=zorder))

        zorder-=1

    return p

def stackbar_plus_minus(axx,xvalues,yvalues,y_order):

    neg_tot = np.zeros(np.shape(xvalues))
    pos_tot = np.zeros(np.shape(xvalues))
    zorder = len(y_order) + 1
    p = []
    for y_name in y_order:
        y_pos = np.copy(yvalues[y_name])
        y_pos[y_pos<0] = 0
        y_neg = np.copy(yvalues[y_name])
        y_neg[y_neg>0] = 0

        p.append(axx.bar(xvalues,y_pos,width=4,bottom=pos_tot))
        axx.bar(xvalues,y_neg,width=4,bottom=neg_tot,color=p[-1][-1].get_facecolor(),label='_nolegend_')
        pos_tot += y_pos
        neg_tot += y_neg

        zorder-=1

    return p


def stacked_flux_plot(nc1,nc2,yr,sday,eday,wk=True,iceup=None,icedur=None,runname1='IC',runname2='ITC',
                      legend_loc='upper right',legend_ax=0,title=""):
    k = kei_read(nc1)
    k2 = kei_read(nc2)
    dn = kei_dn(k,yr)
    dn_day = dn[0::24]

    months = mdates.MonthLocator()
    monthsFmt = mdates.DateFormatter('%b')

    f1 = daily_fluxes(k,sday,eday)
    f2 = daily_fluxes(k2,sday,eday)

    fdiff = {}
    for k in f1.keys():
        fdiff[k] = f2[k] - f1[k]

    fig,axx = plt.subplots(3,1,sharex=True,figsize=(6,10))
    axx=axx.ravel()
    for i,f in enumerate([f1,f2,fdiff]):
        axx[i].stackplot(wd(dn_day[sday:eday],wk), wd(f['pyc2pbl'],wk),wd(f['pbl2pyc'],wk),wd(f['pbl2ice'],wk),
                     wd(f['atm2ice'],wk),wd(f['ice2atm'],wk),wd(f['atm2pbl'],wk),wd(f['pbl2atm'],wk),
                     labels = ['Pycnocline -> PBL','PBL -> Pycnocline','PBL -> Sea-ice','Atmosphere -> Sea-ice',
                            'Sea-ice -> Atmosphere','Atmosphere -> PBL','PBL -> Atmosphere']
                     )
    axx[2].xaxis.set_major_locator(months)
    axx[2].xaxis.set_major_formatter(monthsFmt)
    axx[0].set_ylim([0,270])
    axx[1].set_ylim([0,270])
    axx[1].set_ylabel('Directional Flux [$W  m^{-2}$]')
    axx[legend_ax].legend(loc=legend_loc,fontsize = 'small')

    axx[0].text(0.1,0.87,'%i'%yr+': %s'%runname1,fontsize='large',transform=axx[0].transAxes)
    axx[1].text(0.1,0.87,'%i'%yr+': %s'%runname2,fontsize='large',transform=axx[1].transAxes)
    axx[2].text(0.1,0.87,'%i'%yr+' Difference: %s - %s'%(runname2,runname1),fontsize='large',transform=axx[2].transAxes)

    if iceup is not None:
        dniceup = date2num(iceup)
        dnmelt = date2num(iceup + icedur)
        for ax in axx:
            ylims = ax.get_ylim()
            ax.plot([dniceup,dniceup],ylims,'k:')
            ax.plot([dnmelt,dnmelt],ylims,'k:')

    return fig,axx



def stacked_flux_plot_w_negatives(nc1,nc2,yr,sday,eday,wk=True,iceup=None,icedur=None,runname1='IC',runname2='ITC',
                      legend_loc='upper right',legend_ax=0,title=""):
    k = kei_read(nc1)
    k2 = kei_read(nc2)
    dn = kei_dn(k,yr)
    dn_day = dn[0::24]

    months = mdates.MonthLocator()
    monthsFmt = mdates.DateFormatter('%b')

    f1 = daily_fluxes(k,sday,eday)
    f2 = daily_fluxes(k2,sday,eday)

    fdiff = {}
    for k in f1.keys():
        fdiff[k] = f2[k] - f1[k]

    fig,axx = plt.subplots(3,1,sharex=True,figsize=(6,10))
    axx=axx.ravel()

    labels = ['Pycnocline -> PBL','PBL -> Pycnocline','PBL -> Sea-ice','Atmosphere -> Sea-ice',
                            'Sea-ice -> Atmosphere','Atmosphere -> PBL','PBL -> Atmosphere']
    for i,f in enumerate([f1,f2,fdiff]):

        # p = axx[i].stackplot(wd(dn_day[sday:eday],wk,'all'), wd(f['pyc2pbl'],wk,'pos'),wd(f['pbl2pyc'],wk,'pos'),
        #              wd(f['pbl2ice'],wk,'pos'),wd(f['atm2ice'],wk,'pos'),wd(f['ice2atm'],wk,'pos'),
        #              wd(f['atm2pbl'],wk,'pos'),wd(f['pbl2atm'],wk,'pos'),
        #              labels = labels
        #              )
        # colors = [p1.get_facecolor()[0] for p1 in p]
        # axx[i].stackplot(wd(dn_day[sday:eday],wk,'all'), wd(f['pyc2pbl'],wk,'neg'),wd(f['pbl2pyc'],wk,'neg'),
        #              wd(f['pbl2ice'],wk,'neg'),wd(f['atm2ice'],wk,'neg'),wd(f['ice2atm'],wk,'neg'),
        #              wd(f['atm2pbl'],wk,'neg'),wd(f['pbl2atm'],wk,'neg'),
        #              labels = ['_nolegend_','_nolegend_','_nolegend_','_nolegend_',
        #                     '_nolegend_','_nolegend_','_nolegend_'],
        #              colors = colors
        #              )

        #if i==2:
        #    fig2 = plt.figure()
        #    plt.plot(wd(dn_day[sday:eday],wk), np.array([wd(f['pyc2pbl'],wk),wd(f['pbl2pyc'],wk),wd(f['pbl2ice'],wk),
        #             wd(f['atm2ice'],wk),wd(f['ice2atm'],wk),wd(f['atm2pbl'],wk),wd(f['pbl2atm'],wk)]).transpose())
        #    plt.savefig('2007_fluxes_FREE_TC_2022_03_31_lineplot.png', dpi=1200, bbox_inches='tight')

        f_wk = {k:wd(v,wk) for k,v in f.items()}
        p = stackbar_plus_minus(axx[i],wd(dn_day[sday:eday],wk,'all'),f_wk,
                                 ['pyc2pbl','pbl2pyc','pbl2ice','atm2ice','ice2atm','atm2pbl','pbl2atm'])

    axx[2].xaxis.set_major_locator(months)
    axx[2].xaxis.set_major_formatter(monthsFmt)
    axx[0].set_ylim([0,270])
    axx[1].set_ylim([0,270])
    #axx[2].set_ylim([-270,270])
    axx[1].set_ylabel('Directional Flux [$W  m^{-2}$]')
    axx[legend_ax].legend(labels,loc=legend_loc,fontsize = 'small')

    axx[0].text(0.1,0.87,'%i'%yr+': %s'%runname1,fontsize='large',transform=axx[0].transAxes)
    axx[1].text(0.1,0.87,'%i'%yr+': %s'%runname2,fontsize='large',transform=axx[1].transAxes)
    axx[2].text(0.1,0.87,'%i'%yr+' Difference: %s - %s'%(runname2,runname1),fontsize='large',transform=axx[2].transAxes)

    if iceup is not None:
        dniceup = date2num(iceup)
        dnmelt = date2num(iceup + icedur)
        for ax in axx:
            ylims = ax.get_ylim()
            ax.plot([dniceup,dniceup],ylims,'k:')
            ax.plot([dnmelt,dnmelt],ylims,'k:')

    return fig,axx



def mooring_T_plot(forcing_start_dt=datetime.datetime(2007,1,1),
                   forcing_netcdf_fname=r'/Users/blsaenz/KEI_run/DATA_moorings/kf_300_07_ECMWF.nc'):
    nc = netCDF4.Dataset(forcing_netcdf_fname,'r')
    mooring_t = nc.variables['wct'][...].filled(np.nan)
    ntime,nz = np.shape(mooring_t)
    hour = datetime.timedelta(hours=1)
    time_dn = date2num([forcing_start_dt+hour*i for i in range(ntime)])

    mooring_t[mooring_t==-999.] = np.nan

    fig,axx=plt.subplots(3,1,figsize=(6,7))
    pcs = []

    for i,year in enumerate([2007,2008,2011]):
        yr_b_1 = date2num(datetime.datetime(year, 1, 1))
        yr_b_2 = date2num(datetime.datetime(year + 1, 1, 1))

        t1 = np.argmin(np.abs(time_dn - yr_b_1))
        t2 = np.argmin(np.abs(time_dn - yr_b_2))

        pc = axx[i].pcolormesh(time_dn[t1:t2],np.arange(nz)*(-1.0),np.transpose(mooring_t[t1:t2,:]),cmap=cmocean.cm.thermal,
                          shading='nearest')
        pcs.append(pc)
        months = mdates.MonthLocator()
        monthsFmt = mdates.DateFormatter('%b')
        axx[i].xaxis.set_major_locator(months)
        axx[i].xaxis.set_major_formatter(monthsFmt)

        if year==2011:
            iceup = datetime.datetime(2011, 6, 25)
            icedur = datetime.timedelta(days=158)
        elif year==2008:
            iceup = datetime.datetime(2008, 7, 6)
            icedur = datetime.timedelta(days=98)
        elif year==2007:
            iceup = datetime.datetime(2007, 6, 29)
            icedur = datetime.timedelta(days=115)
        x1 = date2num(iceup)
        x2 = date2num(iceup+icedur)
        axx[i].plot([x1,x1],axx[i].get_ylim(),'k:')
        axx[i].plot([x2,x2],axx[i].get_ylim(),'k:')

        axx[i].text(time_dn[t1]+8, -50, '%i'%year, fontsize=14)

    pos = axx[0].get_position()
    cax = fig.add_axes([pos.x0,0.92,pos.width,0.035])
    plt.colorbar(pcs[0],cax=cax,orientation='horizontal')
    fig.text(0.08, 0.93, r'$\degree C$', fontsize=10)
    axx[0].set_xticklabels([])
    axx[1].set_xticklabels([])
    axx[1].set_ylabel('Depth (m)')
    plt.savefig('Mooring_temps_2021_10_19.png',dpi=1200,bbox_inches='tight')
    plt.show()



def T_diff_pcolor_plot(year,ITC_nc_fname,FREE_nc_fname):
    itc = kei_read(ITC_nc_fname)
    free = kei_read(FREE_nc_fname)

    #t =
    ntime,nz = np.shape(mooring_t)
    hour = datetime.timedelta(hours=1)
    time_dn = date2num([forcing_start_dt+hour*i for i in range(ntime)])

    mooring_t[mooring_t==-999.] = np.nan

    fig,axx=plt.subplots(3,1,figsize=(6,7))
    pcs = []

    for i,year in enumerate([2007,2008,2011]):
        yr_b_1 = date2num(datetime.datetime(year, 1, 1))
        yr_b_2 = date2num(datetime.datetime(year + 1, 1, 1))

        t1 = np.argmin(np.abs(time_dn - yr_b_1))
        t2 = np.argmin(np.abs(time_dn - yr_b_2))

        pc = axx[i].pcolormesh(time_dn[t1:t2],np.arange(nz)*(-1.0),np.transpose(mooring_t[t1:t2,:]),cmap=plt.get_cmap('viridis'),
                          shading='nearest')
        pcs.append(pc)
        months = mdates.MonthLocator()
        monthsFmt = mdates.DateFormatter('%b')
        axx[i].xaxis.set_major_locator(months)
        axx[i].xaxis.set_major_formatter(monthsFmt)

        if year==2011:
            iceup = datetime.datetime(2011, 6, 25)
            icedur = datetime.timedelta(days=158)
        elif year==2008:
            iceup = datetime.datetime(2008, 7, 6)
            icedur = datetime.timedelta(days=98)
        elif year==2007:
            iceup = datetime.datetime(2007, 6, 29)
            icedur = datetime.timedelta(days=115)
        x1 = date2num(iceup)
        x2 = date2num(iceup+icedur)
        axx[i].plot([x1,x1],axx[i].get_ylim(),'r:')
        axx[i].plot([x2,x2],axx[i].get_ylim(),'r:')

        axx[i].text(time_dn[t1]+8, -50, '%i'%year, fontsize=14)

    pos = axx[0].get_position()
    cax = fig.add_axes([pos.x0,0.92,pos.width,0.035])
    plt.colorbar(pcs[0],cax=cax,orientation='horizontal')
    fig.text(0.08, 0.93, r'$\degree C$', fontsize=10)
    axx[0].set_xticklabels([])
    axx[1].set_xticklabels([])
    axx[1].set_ylabel('Depth (m)')
    plt.savefig('Mooring_temps_2021_10_04.png',dpi=1200,bbox_inches='tight')
    plt.show()


def T_diff_pcolor_plot(ITC_nc_fname,IC_nc_fname):
    itc = kei_read(ITC_nc_fname)
    ic = kei_read(IC_nc_fname)

def kei_read_xr(kei_nc):
    return xr.open_dataset(kei_nc)

if __name__ == '__main__':

    ds = kei_read_xr("/Users/blsaenz/temp/keipy_output/keipytest_macmods_2025-05-13_082010/keipytest_macmods.nc")
    dude=1


    for y in [2007,2008,2011]:

        fname = r'/Users/blsaenz/data/KEI/output/v12/300.100_'

        k = kei_read_legacy(fname + '%i'%y + '_IC.nc')
        f, axx = plt.subplots(2, 1)
        axx[0].plot(k['day'], k['deep_pyc'])
        ax0 = axx[0].twinx()
        ax0.plot(k['day'], k['ice_ocean_bottom_flux'],'C1')
        axx[1].plot(k['day'],k['fice'])
        ax1 = axx[1].twinx()
        ax1.plot(k['day'], k['total_ice_freeze'],'C2')
        axx[0].set_title('%i IC'%y)
        axx[0].set_ylim([-70,-150])
        ax0.set_ylim([-10,210])
        ax1.set_ylim([-1000,5.e5+1000])

        k = kei_read_legacy(fname + '%i'%y + '_ITC.nc')
        f, axx = plt.subplots(2, 1)
        axx[0].plot(k['day'], k['deep_pyc'])
        ax0 = axx[0].twinx()
        ax0.plot(k['day'], k['ice_ocean_bottom_flux'],'C1')
        axx[1].plot(k['day'],k['fice'])
        ax1 = axx[1].twinx()
        ax1.plot(k['day'], k['total_ice_freeze'], 'C2')
        axx[0].set_title('%i ITC'%y)
        axx[0].set_ylim([-70,-150])
        ax0.set_ylim([-10,210])
        ax1.set_ylim([-1000,5.e5+1000])

    dude=1

    exit()


    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='IC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    exit()
    plt.savefig('2007_bar_fluxes_IC_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')



    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_TC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='FREE', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_bar_fluxes_FREE_TC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_TC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='FREE', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_bar_fluxes_FREE_TC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_TC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='FREE', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_bar_fluxes_FREE_TC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_IC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='FREE', runname2='IC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_bar_fluxes_FREE_IC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_IC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='FREE', runname2='IC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_bar_fluxes_FREE_IC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_IC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='FREE', runname2='IC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_bar_fluxes_FREE_IC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_TC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='TC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_bar_fluxes_TC_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_TC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='TC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_bar_fluxes_TC_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_TC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='TC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_bar_fluxes_TC_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_TC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='IC', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_bar_fluxes_IC_TC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_TC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='IC', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_bar_fluxes_IC_TC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_TC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='IC', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_bar_fluxes_IC_TC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='IC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_bar_fluxes_IC_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='IC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_bar_fluxes_IC_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='IC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_bar_fluxes_IC_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='FREE', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_bar_fluxes_FREE_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='FREE', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_bar_fluxes_FREE_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_ITC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='FREE', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_bar_fluxes_FREE_ITC_2022_03_31.png', dpi=1200, bbox_inches='tight')
    exit()


    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_TC.nc'
    fig, ax = stacked_flux_plot_w_negatives(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='FREE', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_stacked_fluxes_FREE_TC_2022_03_31.png', dpi=1200, bbox_inches='tight')
    fig, ax = stacked_flux_plot(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                            icedur=datetime.timedelta(days=127), runname1='FREE', runname2='TC',
                                            legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_stacked_fluxes_FREE_TC_2022_03_31_orig.png', dpi=1200, bbox_inches='tight')
    exit()


    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_TC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='FREE', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_stacked_fluxes_FREE_TC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_TC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='FREE', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_stacked_fluxes_FREE_TC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_TC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='FREE', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_stacked_fluxes_FREE_TC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_IC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='FREE', runname2='IC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_stacked_fluxes_FREE_IC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_IC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='FREE', runname2='IC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_stacked_fluxes_FREE_IC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_IC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='FREE', runname2='IC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_stacked_fluxes_FREE_IC_2021_11_25.png', dpi=1200, bbox_inches='tight')
    # exit()

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_TC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='TC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_stacked_fluxes_TC_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_TC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='TC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_stacked_fluxes_TC_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_TC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='TC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_stacked_fluxes_TC_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')
    # exit()

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_TC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='IC', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_stacked_fluxes_IC_TC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_TC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='IC', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_stacked_fluxes_IC_TC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_TC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='IC', runname2='TC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_stacked_fluxes_IC_TC_2021_11_25.png', dpi=1200, bbox_inches='tight')
    # exit()

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='IC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_stacked_fluxes_IC_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='IC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_stacked_fluxes_IC_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_IC.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='IC', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_stacked_fluxes_IC_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')
    # exit()

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2007_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2007, 47, 475, True, iceup=datetime.datetime(2007, 6, 29),
                                icedur=datetime.timedelta(days=127), runname1='FREE', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2007_stacked_fluxes_FREE_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')
    # exit()

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2008_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2008, 47, 475, True, iceup=datetime.datetime(2008, 7, 6),
                                icedur=datetime.timedelta(days=102), runname1='FREE', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2008_stacked_fluxes_FREE_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')
    # exit()

    fn = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_FREE.nc'
    fn2 = r'/Users/blsaenz/data/KEI/output/v12/300.100_2011_ITC.nc'
    fig, ax = stacked_flux_plot(fn, fn2, 2011, 47, 475, True, iceup=datetime.datetime(2011, 6, 25),
                                icedur=datetime.timedelta(days=158), runname1='FREE', runname2='ITC',
                                legend_loc='upper right', legend_ax=0)
    plt.savefig('2011_stacked_fluxes_FREE_ITC_2021_11_25.png', dpi=1200, bbox_inches='tight')

    exit()

    #mooring_T_plot()
    #exit()

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_TC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2007,47,475,True,iceup=datetime.datetime(2007,6,29),
                               icedur=datetime.timedelta(days=127),runname1='FREE',runname2='TC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2007_stacked_fluxes_FREE_TC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_TC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2008,47,475,True,iceup=datetime.datetime(2008,7,6),
                               icedur=datetime.timedelta(days=102),runname1='FREE',runname2='TC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2008_stacked_fluxes_FREE_TC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_TC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2011,47,475,True,iceup=datetime.datetime(2011,6,25),
                               icedur=datetime.timedelta(days=158),runname1='FREE',runname2='TC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2011_stacked_fluxes_FREE_TC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_IC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2007,47,475,True,iceup=datetime.datetime(2007,6,29),
                               icedur=datetime.timedelta(days=127),runname1='FREE',runname2='IC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2007_stacked_fluxes_FREE_IC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_IC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2008,47,475,True,iceup=datetime.datetime(2008,7,6),
                               icedur=datetime.timedelta(days=102),runname1='FREE',runname2='IC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2008_stacked_fluxes_FREE_IC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_IC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2011,47,475,True,iceup=datetime.datetime(2011,6,25),
                               icedur=datetime.timedelta(days=158),runname1='FREE',runname2='IC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2011_stacked_fluxes_FREE_IC_2021_11_09.png',dpi=1200,bbox_inches='tight')
    #exit()


    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_TC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2007,47,475,True,iceup=datetime.datetime(2007,6,29),
                               icedur=datetime.timedelta(days=127),runname1='TC',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2007_stacked_fluxes_TC_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_TC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2008,47,475,True,iceup=datetime.datetime(2008,7,6),
                               icedur=datetime.timedelta(days=102),runname1='TC',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2008_stacked_fluxes_TC_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_TC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2011,47,475,True,iceup=datetime.datetime(2011,6,25),
                               icedur=datetime.timedelta(days=158),runname1='TC',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2011_stacked_fluxes_TC_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')
    #exit()

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_TC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2007,47,475,True,iceup=datetime.datetime(2007,6,29),
                               icedur=datetime.timedelta(days=127),runname1='IC',runname2='TC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2007_stacked_fluxes_IC_TC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_TC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2008,47,475,True,iceup=datetime.datetime(2008,7,6),
                               icedur=datetime.timedelta(days=102),runname1='IC',runname2='TC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2008_stacked_fluxes_IC_TC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_TC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2011,47,475,True,iceup=datetime.datetime(2011,6,25),
                               icedur=datetime.timedelta(days=158),runname1='IC',runname2='TC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2011_stacked_fluxes_IC_TC_2021_11_09.png',dpi=1200,bbox_inches='tight')
    #exit()

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2007,47,475,True,iceup=datetime.datetime(2007,6,29),
                               icedur=datetime.timedelta(days=127),runname1='IC',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2007_stacked_fluxes_IC_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2008,47,475,True,iceup=datetime.datetime(2008,7,6),
                               icedur=datetime.timedelta(days=102),runname1='IC',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2008_stacked_fluxes_IC_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2011,47,475,True,iceup=datetime.datetime(2011,6,25),
                               icedur=datetime.timedelta(days=158),runname1='IC',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2011_stacked_fluxes_IC_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')
    #exit()


    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2007,47,475,True,iceup=datetime.datetime(2007,6,29),
                               icedur=datetime.timedelta(days=127),runname1='FREE',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2007_stacked_fluxes_FREE_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')
    #exit()

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2008,47,475,True,iceup=datetime.datetime(2008,7,6),
                               icedur=datetime.timedelta(days=102),runname1='FREE',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2008_stacked_fluxes_FREE_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')
    #exit()

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2011,47,475,True,iceup=datetime.datetime(2011,6,25),
                               icedur=datetime.timedelta(days=158),runname1='FREE',runname2='ITC',
                               legend_loc='upper right',legend_ax=0)
    plt.savefig('2011_stacked_fluxes_FREE_ITC_2021_11_09.png',dpi=1200,bbox_inches='tight')
    exit()


    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_FREE.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_ITC.nc'
    fig, ax = T_diff_pcolor_plot('2007', fn2, fn)



    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2011_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2011,47,475,True,iceup=datetime.datetime(2011,6,25),
                               icedur=datetime.timedelta(days=158))
    plt.savefig('2011_stacked_fluxes_2021_10_04.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2008_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2008,47,475,True,iceup=datetime.datetime(2008,7,6),
                               icedur=datetime.timedelta(days=98))
    plt.savefig('2008_stacked_fluxes_2021_10_04.png',dpi=1200,bbox_inches='tight')

    fn = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_IC.nc'
    fn2 = r'/Volumes/reservior/data/KEI/output/moorings_paper_v10/300.100_2007_ITC.nc'
    fig,ax = stacked_flux_plot(fn,fn2,2007,47,475,True,iceup=datetime.datetime(2007,6,29),
                               icedur=datetime.timedelta(days=115))
    plt.savefig('2007_stacked_fluxes_2021_10_04.png',dpi=1200,bbox_inches='tight')

    dude=1

