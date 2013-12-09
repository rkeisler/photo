import numpy as np
import matplotlib.pylab as pl
pl.ion()
import pdb
import cPickle as pickle
import ipdb

wise_datapath='/home/rkeisler/wise/'

def colmag_wise_boss(frac=0.0065):
    # frac=0.0065 corresponds to 2.0 square degrees.
    colcol_subset_boss(frac=frac, sample='lowz', nofirst=False, sym1='k.', sym2='b.')
    colcol_subset_boss(frac=frac, sample='cmass', nofirst=True, sym2='g.')
    pl.legend(['WISE','BOSS LOWZ','BOSS CMASS'], 'upper left')

def colcol_subset_boss(frac=0.0065, sample='cmass', nofirst=False, sym1='b.', sym2='r.'):
    label, boss_match_ind, min_dist_sq, x = match_wise_boss(wise_part=44, sample=sample, quick=True)
    w1 = x[:, 2]
    w2 = x[:, 3]
    w3 = x[:, 4]
    w4 = x[:, 5]
    wh_boss = np.where(label==1)[0]
    nshow1 = int(frac*len(w1))
    nshow2 = int(frac*len(wh_boss))
    if not(nofirst): colmag_plot(w1,w2,w3,w4,nshow=nshow1,symbol=sym1)
    colmag_plot(w1[wh_boss],w2[wh_boss],w3[wh_boss],w4[wh_boss],nshow=nshow2, symbol=sym2, noclf=True)

    # explore this color cut
    cA = w1
    cB = w2-w4-w1
    wh_cut = np.where((cA>=14.8) & (cA<=16.0) & (cB>=-10.0) & (cB<=-8.8))[0]
    wh_lots = np.where((cA>=14.8) & (cA<=16.0) & (cB>=-10.0) & (cB<=-8.8) & (label==0))[0]
    wh_few = np.where((cA>=14.8) & (cA<=16.0) & (cB>=-10.0) & (cB<=-8.8) & (label==1))[0]
    print 1.*len(wh_lots)/len(wh_few),len(wh_cut)



def colcol_subset(nshow=10000):
    x = np.load(wise_datapath+'wise-allsky-cat-part44.npy')
    w1 = x[:, 2]
    w2 = x[:, 3]
    w3 = x[:, 4]
    w4 = x[:, 5]    
    colcol_plot(w1,w2,w3,w4,nshow=nshow)


def colmag_plot(w1_,w2_,w3_,w4_,nshow=10000,symbol='b.',noclf=False):
    ngals = len(w1_)
    ind = np.arange(ngals)
    np.random.shuffle(ind)
    ind = ind[0:nshow]
    w1 = w1_[ind]
    w2 = w2_[ind]
    w3 = w3_[ind]
    w4 = w4_[ind]
    d1 = {'color':w1, 'name':'w1', 'range':[8,20]}    
    d12 = {'color':w1-w2, 'name':'w1-w2', 'range':[-2,3]}
    d23 = {'color':w2-w4-w1, 'name':'w2-w4-w1', 'range':[-14,-6]}
#    d23 = {'color':w2-w4, 'name':'w2-w4', 'range':[-2,10]}    
    d34 = {'color':w3-w4, 'name':'w3-w4', 'range':[-1,5]}
    d14 = {'color':w1-w3, 'name':'w1-w3', 'range':[-1,11]}
    d = [d1, d12, d23, d34, d14]
    if not(noclf): pl.clf()
    pl.figure(1,figsize=(15,8))
    count=0
    fs=14
    from itertools import combinations
    for a,b in combinations(d,2):
        count+=1
        pl.subplot(2,5,count)
        pl.plot(a['color'], b['color'], symbol)
        pl.xlabel(a['name'], fontsize=fs)
        pl.ylabel(b['name'], fontsize=fs)
        pl.xlim(a['range'])
        pl.ylim(b['range'])        

    
def colcol_plot(w1_,w2_,w3_,w4_,nshow=10000,symbol='b.',noclf=False):
    ngals = len(w1_)
    ind = np.arange(ngals)
    np.random.shuffle(ind)
    ind = ind[0:nshow]
    w1 = w1_[ind]
    w2 = w2_[ind]
    w3 = w3_[ind]
    w4 = w4_[ind]
    d12 = {'color':w1-w2, 'name':'w1-w2', 'range':[-2,3]}
    d23 = {'color':w2-w3, 'name':'w2-w3', 'range':[-2,8]}
    d34 = {'color':w3-w4, 'name':'w3-w4', 'range':[-1,5]}
    d14 = {'color':w1-w4, 'name':'w1-w4', 'range':[-1,11]}
    d = [d12, d23, d34, d14]
    if not(noclf): pl.clf()
    pl.figure(1,figsize=(15,8))
    count=0
    fs=14
    from itertools import combinations
    for a,b in combinations(d,2):
        count+=1
        pl.subplot(2,3,count)
        pl.plot(a['color'], b['color'], symbol)
        pl.xlabel(a['name'], fontsize=fs)
        pl.ylabel(b['name'], fontsize=fs)
        pl.xlim(a['range'])
        pl.ylim(b['range'])        
               

def dumb(sample='cmass'):
    min_dist=study_subset(sample=sample)
    dd=np.linspace(0.5,15,30)
    frac = np.zeros_like(dd)
    for i,ddd in enumerate(dd): frac[i]=1.*len(np.where(np.sqrt(min_dist)*3600.<=ddd)[0])/1000.
    pl.clf()
    pl.plot(dd,frac)
    pl.plot(dd,frac,'bo')
    pl.ylabel('Fraction of '+sample.upper()+' galaxies matched to WISE objects',fontsize=14)
    pl.xlabel('Match Radius (arcsec)',fontsize=14)
    pl.ylim(0.8,1.0)


def match_wise_boss(wise_part=44, sample='cmass',match_rad_arcsec=2.0, quick=False):
    savename = wise_datapath+'match_wise%i_boss'%wise_part+sample.upper()+'_%0.1frad.pkl'%match_rad_arcsec
    if quick:
        label, boss_match_ind, min_dist_sq, X = pickle.load(open(savename,'r'))
        return label, boss_match_ind, min_dist_sq, X

    from time import time
    # get BOSS positions
    ra_boss, dec_boss = parse_boss(sample=sample)

    # get WISE positions
    X = get_wise_part(wise_part)
    ra_wise_tmp = X[:,0]
    dec_wise_tmp = X[:,1]

    # the WISE positions are by definition bounded in a tight range of
    # declination.  Let's only try to match BOSS sourecs that fall in
    # that dec range.  Similarly, only try to match WISE sources that
    # fall in the RA range of the subset of BOSS sources in this DEC
    # range.
    wh_dec = np.where((dec_boss>=np.min(dec_wise_tmp)) & (dec_boss<=np.max(dec_wise_tmp)))[0]
    ra_boss = ra_boss[wh_dec]
    dec_boss = dec_boss[wh_dec]
    nboss = len(ra_boss)
    wh_ra = np.where((ra_wise_tmp>=np.min(ra_boss)) & (ra_wise_tmp<=np.max(ra_boss)))[0]
    X = X[wh_ra, :]
    ra_wise = X[:,0]
    dec_wise = X[:,1]
    nwise = len(ra_wise)

    # loop over the BOSS sources and find matches.
    min_dist_sq = []
    label = np.zeros(nwise,dtype=int)
    boss_match_ind = np.zeros(nwise,dtype=int)
    match_deg_sq = (match_rad_arcsec/3600.)**2.
    count=-1
    timea = time()
    for iboss in range(nboss):
#    for iboss in range(1000):    
        count+=1
        if (count % 100)==0:
            print count,nboss
            print 'took ',time()-timea
            timea=time()
        this_ra = ra_boss[iboss]
        this_dec = dec_boss[iboss]
        dist_sq = (this_dec-dec_wise)**2. + ((this_ra-ra_wise)*np.cos(this_dec*np.pi/180.))**2.
        this_min_dist_sq = np.min(dist_sq)
        if (this_min_dist_sq <= match_deg_sq):
            this_wise_ind = np.argmin(dist_sq)
            label[this_wise_ind] = 1
            boss_match_ind[this_wise_ind] = iboss
            #print np.sqrt(this_min_dist_sq)*3600.
        min_dist_sq.append(np.min(dist_sq))
    pickle.dump((label, boss_match_ind, min_dist_sq, X), open(savename,'w'))
    return label, boss_match_ind, min_dist_sq, X


def parse_boss(sample='cmass'):
    import pyfits
    x=pyfits.open(wise_datapath+'bosstile-final-collated-boss2-boss32.fits')
    tbdata = x[1].data 
    ra=tbdata.field(12)
    dec=tbdata.field(13)
    wh_lowz=np.where(((tbdata.field(44)&1)>0) & (tbdata.field(44)>0))[0]
    wh_cmass=np.where(((tbdata.field(44)&2)>0) & (tbdata.field(44)>0))[0]
    if sample=='lowz': wh_use=wh_lowz
    if sample=='cmass': wh_use=wh_cmass
    return ra[wh_use], dec[wh_use]


def get_wise_part(part=44):
    savename=wise_datapath+'wise-allsky-cat-part44.npy'
    return np.load(open(savename,'r'))


def boil_down_wise(filename=wise_datapath+'wise-allsky-cat-part44', 
                   savepath=wise_datapath, npersub=100000):
    savename = savepath+filename.split('/')[-1]
    f=open(filename)
    default_mag = 30.
    w_ind = wise_cat_ind()
    np_ind = np_array_ind()
    nsave = len(np_array_ind())
    count = 0
    subcount = 0
    subname = 0
    X = []
    for line in f:
        count+=1
        subcount+=1
        if (count % 100000)==0: print count
        tmp = line.split('|')
        thisX = np.zeros(nsave)
        for k in w_ind:
            valtmp = tmp[w_ind[k]]
            if k=='cc_flags':
                if valtmp=='': valtmp=1
                else: valtmp = not(valtmp=='0000')
            else:
                if valtmp=='': valtmp=default_mag
                else: valtmp=float(valtmp)
            thisX[np_ind[k]]=valtmp
        X.append(thisX)        
        if (subcount>=npersub):
            subcount=0
            np.save(savename+'_%i'%subname, np.vstack(X))            
            subname += 1
            X=[]
    np.save(savename+'_%i'%subname, np.vstack(X))            


def wise_cat_ind():
    return {'ra':1, 'dec':2,
            'w1':16, 'w2':20, 'w3':24, 'w4':28,
            'j':287, 'h':289, 'k':291,
            'cc_flags':55}

def np_array_ind():
    tmp={'ra':0, 'dec':1,
         'w1':2, 'w2':3, 'w3':4, 'w4':5,
         'j':6, 'h':7, 'k':8,
         'cc_flags':9}
    if len(tmp)!=max(tmp.values())+1:
        print 'huh?'
        pdb.set_trace()
    return tmp


def process_wise_dec_strip(part):
    from os import system
    from time import time
    timeo = time()

    savename_prefix = 'wise-allwise-cat-part%02d'%part
    savepath = '/home/rkeisler/wise/tmp/'
    bz_name = savepath+savename_prefix+'.bz2'
    unzip_name = savepath+savename_prefix

    print '...downloading...'
    this_address = 'http://irsadist.ipac.caltech.edu/wise-allwise/'+savename_prefix+'.bz2'
    cmd = 'curl -o '+bz_name+' '+this_address
    print cmd
    timea = time()
    system(cmd)
    timeb = time()
    print 'that took ',timeb-timea

    print '...unzipping...'
    cmd = 'bunzip2 '+bz_name
    print cmd
    timea = time()
    system(cmd)
    timeb = time()
    print 'that took ',timeb-timea

    print '...reducing...'
    print "boil_down_wise(filename=unzip_name, savepath='/home/rkeisler/wise/allwise/')"
    print unzip_name
    timea = time()
    boil_down_wise(filename=unzip_name, savepath='/home/rkeisler/wise/allwise/')
    timeb = time()
    print 'that took ',timeb-timea

    print '...removing raw data...'
    cmd = 'rm '+unzip_name
    print cmd
    system(cmd)
    print '************************************************'
    print 'one part, EVERYTHING took ',time()-timeo
    print '************************************************'


def process_many(imin_incl=1, imax_incl=48):
    for i in range(imin_incl, imax_incl+1): process_wise_dec_strip(i)


def read_cosmos(cosmosname='/Users/rkeisler/cosmos_zphot_mag25.tbl.tbl'):
    f=open(cosmosname,'r')
    ra = []
    dec = []
    z = []
    type = []
    for line in f:
        if '|' in line: continue
        if '\\' in line: continue
        tmp=line.split()
        if tmp[2]=='null': continue
        if tmp[4]=='null':
            # these are in masked areas (bright stars, etc.).
            tmp[4]=0.
            tmp[5]=-99
        ra.append(float(tmp[2]))
        dec.append(float(tmp[3]))
        z.append(float(tmp[4]))
        type.append(int(tmp[5]))
    f.close()
    ra = np.array(ra)
    dec = np.array(dec)
    z = np.array(z)
    type = np.array(type, dtype=int)
    return ra, dec, z, type


def read_wise_cat_cosmos(filename='/Users/rkeisler/wise_cat_COSMOS.txt'):
    f=open(filename,'r')
    ra=[]
    dec=[]
    w1=[]
    w2=[]
    w3=[]
    w4=[]
    w1sn=[]
    w1cov=[]
    j=[]
    h=[]
    k=[]
    arrived_labels=False
    default_mag=30.
    for line in f:
        if ('|' in line):
            arrived_labels=True
            continue
        else:
            if not(arrived_labels): continue
        tmp=line.split()
        ra.append(float(tmp[1]))
        dec.append(float(tmp[2]))
        w1.append(float(tmp[5]))
        w1sn.append(float(tmp[6]))
        w2.append(float(tmp[7]))
        w3.append(float(tmp[8]))
        w4.append(float(tmp[9]))
        w1cov.append(float(tmp[10]))
        if tmp[11]=='null': j.append(default_mag)
        else: j.append(float(tmp[11]))
        if tmp[12]=='null': h.append(default_mag)
        else: h.append(float(tmp[12]))
        if tmp[13]=='null': k.append(default_mag)
        else: k.append(float(tmp[13]))
    ra=np.array(ra)
    dec=np.array(dec)
    wh_cosmos=np.where((ra>=149.41143) & (ra<=150.82684) & (dec>=1.498815) & (dec<=2.91273))[0]
    ra = ra[wh_cosmos]
    dec = dec[wh_cosmos]
    w1=np.array(w1)[wh_cosmos]
    w2=np.array(w2)[wh_cosmos]
    w3=np.array(w3)[wh_cosmos]
    w4=np.array(w4)[wh_cosmos]
    w1sn=np.array(w1sn)[wh_cosmos]
    w1cov=np.array(w1cov)[wh_cosmos]
    j=np.array(j)[wh_cosmos]
    h=np.array(h)[wh_cosmos]
    k=np.array(k)[wh_cosmos]
    return ra,dec,w1,w2,w3,w4,j,h,k


def study_match_rad():
    ra_c, dec_c, z_c, type_c = read_cosmos()
    ra,dec,w1,w2,w3,w4,j,h,k = read_wise_cat_cosmos()
    count=0
    n_wise = len(ra)
    min_dist=[]
    for this_ra, this_dec in zip(ra,dec):
        count+=1
        if (count % 1000)==0: print count,n_wise
        this_dist = np.sqrt((this_dec-dec_c)**2. + ((this_ra-ra_c)*np.cos(this_dec*np.pi/180.))**2.)*3600.
        this_min_dist = np.min(this_dist)
        this_cosmos_ind = np.argmin(this_dist)
        if (type_c[this_cosmos_ind]!=0): continue
        min_dist.append(this_min_dist)
    min_dist = np.array(min_dist)
    rad=np.linspace(0.2,5.0,49)
    frac=np.zeros_like(rad)
    for i,this_rad in enumerate(rad):
        frac[i] = 1.*len(np.where(min_dist<=this_rad)[0])/len(min_dist)
    pl.plot(rad,frac)
    pl.plot(rad,frac,'bo')
    pl.xlabel('Match Radius (arcsec)',fontsize=14)
    pl.ylabel('Fraction Matched',fontsize=14)
    pl.title('WISE and COSMOS matching',fontsize=16)
    pdb.set_trace()


def match_wise_cosmos(quick=False):
    savename='match_wise_cosmos.pkl'
    if quick: 
        X,z,type=pickle.load(open(savename,'r'))
        return X,z,type
    ra_c, dec_c, z_c, type_c = read_cosmos()
    ra,dec,w1,w2,w3,w4,j,h,k = read_wise_cat_cosmos()
    count=0
    n_wise = len(ra)
    z = np.zeros(n_wise)
    type = np.zeros(n_wise)
    for i in range(n_wise):
        count+=1
        if (count % 1000)==0: print count,n_wise
        this_ra = ra[i]
        this_dec = dec[i]
        this_dist = np.sqrt((this_dec-dec_c)**2. + ((this_ra-ra_c)*np.cos(this_dec*np.pi/180.))**2.)*3600.
        this_cosmos_ind = np.argmin(this_dist)
        z[i] = z_c[this_cosmos_ind]
        type[i] = type_c[this_cosmos_ind]
    X = np.vstack((w1,w2,w3,w4,j,h,k)).T
    whuse = np.where(type!=-99)[0]
    X = X[whuse,:]
    z = z[whuse]
    type = type[whuse]
    pickle.dump((X,z,type), open(savename,'w'))
    return X,z,type
    
        
def study_dndz(cA_cen=15.4, cA_width=1.2,
               cB_cen=-9.0, cB_width=1.2):
    x, z, type = match_wise_cosmos(quick=True)
    w1=x[:,0]
    w2=x[:,1]
    w3=x[:,2]
    w4=x[:,3]
    j=x[:,4]
    h=x[:,5]
    k=x[:,6]

    cA = w1
    cB = w2-w4-w1

    #cA_cen = 15.4
    #cA_width = 1.2
    #cB_cen = -9.4
    #cB_width = 1.2

    cA_min = cA_cen-cA_width/2.
    cA_max = cA_cen+cA_width/2.
    cB_min = cB_cen-cB_width/2.
    cB_max = cB_cen+cB_width/2.
    wh_cut = np.where((cA>=cA_min) & (cA<=cA_max) & (cB>=cB_min) & (cB<=cB_max) & (type==0))[0]
    n_cut = len(wh_cut)
    fcolor='green'

    pl.figure(1)
    pl.clf()
    pl.subplot(2,1,1)
    pl.plot(cA,cB,'k.')
    pl.plot(cA[wh_cut], cB[wh_cut],'r.')
    meanz = np.mean(z[wh_cut])
    medz = np.median(z[wh_cut])
    pl.title('%i gals, <z>=%0.2f, |z|=%0.2f'%(n_cut, meanz, medz))

    pl.subplot(2,1,2)
    zbins = np.linspace(0.1,2.0,15)
    n, bins, patches = pl.hist(z[wh_cut], zbins, normed=0, facecolor=fcolor, alpha=0.5)
    zcen = 0.5*(zbins[0:-1]+zbins[1:])



def study_dndz_interactive():
    x, z, type = match_wise_cosmos(quick=True)
    w1=x[:,0]
    w2=x[:,1]
    w3=x[:,2]
    w4=x[:,3]
    j=x[:,4]
    h=x[:,5]
    k=x[:,6]

    cA = w1
    cB = w2-w4-0.7*w1
    xlab = 'W1'
    ylab = 'W2-W4-0.7W1'
    fs = 16

    from matplotlib import rc
    rc('xtick', labelsize=17) 
    rc('ytick', labelsize=17) 
    fcolor='green'
    fcolor='blue'
    pl.figure(1, figsize=(12,8))
    pl.clf()
    pl.subplot(2,1,1)
    pl.plot(cA,cB,'k.')
    pl.xlabel(xlab,fontsize=fs);pl.ylabel(ylab,fontsize=fs)

    keep_going=True
    from matplotlib.pylab import ginput
    import cosmolopy as cp
    while keep_going:
        pt = ginput(2,timeout=100000)
        pl.clf()
        pl.subplot(2,1,1)
        pl.plot(cA,cB,'k.')
        pl.xlabel(xlab,fontsize=fs);pl.ylabel(ylab,fontsize=fs)
        if pt[0]==pt[1]:
            keep_going=False
        else:
            x0 = pt[0][0]
            y0 = pt[0][1]
            x1 = pt[1][0]
            y1 = pt[1][1]
            cA_min = np.min([x0,x1])
            cA_max = np.max([x0,x1])
            cB_min = np.min([y0,y1])
            cB_max = np.max([y0,y1])
            wh_cut = np.where((cA>=cA_min) & (cA<=cA_max) & (cB>=cB_min) & (cB<=cB_max) & (type==0))[0]
            n_cut = len(wh_cut)

            pl.subplot(2,1,1)
            pl.plot(cA[wh_cut], cB[wh_cut],'r.')
            pl.xlabel(xlab,fontsize=fs);pl.ylabel(ylab,fontsize=fs)
            meanz = np.mean(z[wh_cut])
            medz = np.median(z[wh_cut])
            pl.title('[%0.2f, %0.2f], [%0.2f, %0.2f], %i gals, <z>=%0.2f, |z|=%0.2f, sig=%0.2f'%(cA_min, cA_max, cB_min, cB_max, n_cut, meanz, medz,np.std(z[wh_cut])), fontsize=16)
#            pl.legend(['%0.2f, %0.2f'%(cA_min, cA_max), '%0.2f, %0.2f'%(cB_min, cB_max)], 'upper left')

            pl.subplot(2,1,2)
            zbins = np.linspace(0.1,2.0,15)
            n, bins, patches = pl.hist(z[wh_cut], zbins, normed=0, facecolor=fcolor, alpha=0.5)
            zcen = 0.5*(zbins[0:-1]+zbins[1:])
            rsound_com = 153.
            rsound_phs = rsound_com/(1.+zcen)
            dang = cp.distance.angular_diameter_distance(zcen,**cp.fidcosmo)
            ang_sound = rsound_phs/dang
            ang_sound /= np.median(ang_sound)
            pl.xlabel('Redshift (z)',fontsize=20)
            pl.ylabel('N',fontsize=20)
            #pl.plot(zcen, ang_sound*np.max(n)/4.5, 'b', linewidth=2)
            #print 30./1e3/dang*180./np.pi*3600. # galaxy size

def get_hpix(nside=2**9, w1_min = 16.0, w1_max = 16.7,
             col_min = -999.0, col_max = 999.0, 
             coord='G', quick=False):
    savename = '/home/rkeisler/wise/hpix_'+coord+'_nside%i_w1_%3.1f_%3.1f'%(nside, w1_min, w1_max)+'.pkl'
    if quick: return pickle.load(open(savename,'r'))

    import healpy as hp
    from glob import glob
    files = glob('/home/rkeisler/wise/allwise/*npy')
    npix = hp.nside2npix(nside)
    nmap = np.zeros(npix)
    count=0
    nfiles = len(files)
    for file in files:
        count+=1
        print count,nfiles
        x = np.load(file)
        ra=x[:,0]
        dec=x[:,1]
        w1=x[:,2]
        w2=x[:,3]
        w3=x[:,4]
        w4=x[:,5]
        cc_flags=x[:,9]
        col = w2-w4-0.7*w1
        wh_cut = np.where((w1>=w1_min) & (w1<=w1_max) & (cc_flags==0))[0]
        #wh_cut = np.where((w1>=w1_min) & (w1<=w1_max) & (col>=col_min) & (col<=col_max))[0]        
        phi_cut = ra[wh_cut]*np.pi/180.
        th_cut = np.pi/2.-dec[wh_cut]*np.pi/180.

        # if desired, convert from celestial coordinates.
        if (coord!='C'):
            th_cut_old = th_cut.copy()
            phi_cut_old = phi_cut.copy()
            r = hp.Rotator(coord=['C',coord])
            th_cut, phi_cut = r(th_cut_old, phi_cut_old)
        
        ind = hp.ang2pix(nside, th_cut, phi_cut)
        for thing in ind: nmap[thing]+=1
    pickle.dump(nmap, open(savename,'w'))
    return nmap

def chunk_dict():
    dd={'a':{'ca_min':14.04, 'ca_max':14.71, 'cb_min':-4.94, 'cb_max':-4.15},
        'b':{'ca_min':14.69, 'ca_max':15.15, 'cb_min':-4.92, 'cb_max':-4.25},
        'c':{'ca_min':15.17, 'ca_max':15.58, 'cb_min':-5.00, 'cb_max':-4.09},
        'd':{'ca_min':15.58, 'ca_max':15.91, 'cb_min':-4.94, 'cb_max':-3.76},
        'e':{'ca_min':15.92, 'ca_max':16.34, 'cb_min':-4.94, 'cb_max':-3.58},
        'f':{'ca_min':16.37, 'ca_max':16.85, 'cb_min':-5.15, 'cb_max':-3.47},
        'g':{'ca_min':16.85, 'ca_max':17.25, 'cb_min':-5.23, 'cb_max':-3.68},
        'h':{'ca_min':17.20, 'ca_max':17.48, 'cb_min':-5.20, 'cb_max':-3.99}}
    return dd


def make_many_hpix():
    nside=2**9
    dicts = [{'w1_min':15.7, 'w1_max':16.7},
             {'w1_min':16.0, 'w1_max':16.7},
             {'w1_min':15.7, 'w1_max':18.0},
             {'w1_min':16.7, 'w1_max':18.0},
             {'w1_min':16.0, 'w1_max':18.0}]             

    for d in dicts:
        this_map=get_hpix(nside=nside, 
                          w1_min=d['w1_min'], w1_max=d['w1_max'],
                          quick=False)
        

def mask_from_map(nmap, fwhm_deg=3.0, final_fwhm_deg=5.0, 
                  thresh_min=0.9, thresh_max=3.0,
                  lat_gal_cut=20., ecl_pole_rad=20.,
                  coord='G'):

    import healpy as hp
    mask = np.ones_like(nmap)

    # mask out low galactic latitudes
    npix = len(nmap)
    nside = hp.npix2nside(npix)
    ipix = np.arange(npix)
    theta_c, phi_c = hp.pix2ang(nside, ipix, nest=False)
    r = hp.Rotator(coord=[coord,'G'])  # transforms to galactic
    theta_gal, phi_gal = r(theta_c, phi_c)
    lat_gal_deg = (np.pi/2.-theta_gal)*180./np.pi
    wh_gal = np.where(np.abs(lat_gal_deg)<lat_gal_cut)[0]
    mask[wh_gal]=0.

    # mask out the ecliptic poles
    # actually, isn't the coverage *better* at the poles??
    r = hp.Rotator(coord=[coord,'E'])  # transforms to ecliptic
    theta_ecl, phi_ecl = r(theta_c, phi_c)
    lat_ecl_deg = (np.pi/2.-theta_ecl)*180./np.pi
    wh_ecl_pole = np.where(np.abs(np.abs(lat_ecl_deg)-90.)<ecl_pole_rad)[0]
    mask[wh_ecl_pole]=0.

    # mask out Magellanic clouds.
    # tmpp

    # mask according to <nmap>, which is a proxy for coverage.
    sm = hp.smoothing(nmap*mask, fwhm=fwhm_deg*np.pi/180.)
    wh=np.where((sm<(thresh_min*np.median(sm))))[0]
    mask[wh]=0.
    wh=np.where((sm>(thresh_max*np.median(sm))))[0]
    mask[wh]=0.

    # apply a  final smoothing/threshold/smoothing
    mask = hp.smoothing(mask, fwhm=final_fwhm_deg*np.pi/180.)
    #hp.mollview(mask);ipdb.set_trace()
    mask = hp.smoothing(mask>0.95, fwhm=0.5*final_fwhm_deg*np.pi/180.)
    if False:
        hp.mollview(mask*nmap)
        pdb.set_trace()
    return mask


def many_spectra():
    import healpy as hp
    chunks = ['f','g','h']
    thresh = {'f':1.0, 'g':1.0, 'h':1.8}
    fwhm_deg = {'f':7.0, 'g':7.0, 'h':5.0}
    final_fwhm_deg = {'f':7.0, 'g':7.0, 'h':7.0}
    cl={}
    pl.figure(1)
    pl.clf()
    xlim = [10,700]
    ylim = [1e-7, 3e-4]
    for k in chunks:
        print k
        nmap = get_hpix(quick=True,name=k)
        mask = mask_from_map(nmap, fwhm_deg=fwhm_deg[k], final_fwhm_deg=final_fwhm_deg[k], thresh=thresh[k])
#        mean_nmap = np.mean(nmap[np.where(mask!=0)[0]])
        mean_nmap = np.mean(nmap[np.where(mask>0.5)[0]])
        delta = (nmap-mean_nmap)/mean_nmap
        hp.mollview(hp.smoothing(delta*mask,fwhm=1.*np.pi/180.),title=k,min=-0.3,max=0.3)
        continue
        this_cl = hp.anafast(delta*mask)/np.mean(mask**2.)
        l=np.arange(len(this_cl))
        # smooth in l*cl
        lcl = this_cl*l
        sm_lcl = np.zeros_like(lcl)
        rbox = 15
        for i in l:
            imin = np.max([0,i-rbox])
            imax = np.min([np.max(l), i+rbox])
            sm_lcl[i]=np.mean(lcl[imin:imax])
#        cl[k:this_cl]
#        pl.loglog(l,this_cl,linewidth=2)
        pl.loglog(l,sm_lcl,linewidth=2)
        pdb.set_trace()
    pl.xlim(xlim)
#    pl.ylim(ylim)
    pl.legend(chunks)
    pl.xlabel('L')
    pl.ylabel('L*CL')
#    pl.ylabel('CL')


def get_wise_in_cfhtls(field='f1'):
    from glob import glob

    if field=='f1':
        ra_min = 30.2777134
        ra_max = 38.7223038
        dec_min = -11.128138
        dec_max = -3.8051642
        parts = ['22','23','24','25']

    ra = []
    dec = []
    w1 = []
    w2 = []
    w3 = []
    w4 = []
    for part in parts:
        print part
        files = glob(wise_datapath+'/allwise/wise-allwise-cat-part'+part+'*.npy')
        for file in files:
            x = np.load(file)
            ra_tmp = x[:,0]
            dec_tmp = x[:,1]
            wh = np.where((ra_tmp>ra_min) & (ra_tmp<ra_max) & (dec_tmp>dec_min) & (dec_tmp<dec_max))[0]
            if len(wh)==0: continue
            ra.append(ra_tmp[wh])
            dec.append(dec_tmp[wh])
            w1.append(x[wh,2])
            w2.append(x[wh,3])
            w3.append(x[wh,4])
            w4.append(x[wh,5])            
            
    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    w1 = np.concatenate(w1)
    w2 = np.concatenate(w2)
    w3 = np.concatenate(w3)
    w4 = np.concatenate(w4)

    pickle.dump((ra, dec, w1, w2, w3, w4), open(wise_datapath+'/wise_in_cfhtls_'+field+'.pkl','w'))

    
    
def study_cfhtls_match_rad():
    import cfhtls
    print '...loading CFHTLS...'
    ra_c, dec_c, z_c, i_c = cfhtls.parse_cfhtls(1, quick=True)
    print '...loading WISE in CFHTLS...'
    ra,dec,w1,w2,w3,w4 = pickle.load(open('/Users/rkeisler/wise_in_cfhtls_f1.pkl','r'))
    from scipy.spatial import KDTree
    print '...building KDTree...'
    tree = KDTree(np.vstack((ra_c, dec_c)).T)
    print '...querying tree...'
    x = tree.query(np.vstack((ra, dec)).T)
    dist_arcsec = x[0]*3600.
    ind_c = x[1]

    pickle.dump((tree, ind_c, dist_arcsec), open('tree_ind_dist.pkl','w'))

    my_wise_color = w2-w4-0.7*w1
    wh_wise_cut = np.where((w1>16)&(w1<17)&(my_wise_color>-5.0)&(my_wise_color<-3.5))[0]
    ind_wise_cut = ind_c[wh_wise_cut]
    ipdb.set_trace()

        
def correlate_with_planck_lensing():
    import pyfits
    import healpy as hp
    # get wise stuff
    print '...loading wise...'
    nmap = get_hpix(nside=2**9, w1_min = 16.0, w1_max = 16.7,
                    coord='G', quick=True)
    mask_wise = mask_from_map(nmap, coord='G',thresh_min=0.9, thresh_max=2.3)
    whok = np.where(mask_wise>0.5)[0]
    nbar = np.median(nmap[whok])
    delta = (1.*nmap-nbar)/nbar

    # get planck lensing stuff
    print '...loading planck...'
    pdatadir='/home/rkeisler/cib_delens/planck_data/'
    phi = pyfits.getdata(pdatadir+'COM_CompMap_Lensing_2048_R1.10.fits',1)['PHIBAR']
    phi = hp.reorder(phi, out='RING', inp='NESTED')
    mask_lens = pyfits.getdata(pdatadir+'COM_CompMap_Lensing_2048_R1.10.fits')['MASK']
    mask_lens = hp.reorder(mask_lens, out='RING', inp='NESTED')
    # downgrade to wise resolution
    print '...down-grading...'
    nside_wise = hp.npix2nside(len(delta))
    nside_planck = hp.npix2nside(len(phi))
    phi_dg = hp.ud_grade(phi, nside_wise)
    mask_lens_dg = hp.ud_grade(mask_lens, nside_wise)

    mask = mask_wise*mask_lens_dg
    print '...anafast...'
    cl_phi_delta = hp.anafast(mask*delta, map2=mask*phi_dg)/np.mean(mask**2.)
    # correct for phi transfer function
    tf_phi = pyfits.getdata(pdatadir+'COM_CompMap_Lensing_2048_R1.10.fits',2)['RLPP']
    cl_phi_delta /= tf_phi[0:len(cl_phi_delta)]
    nl = len(cl_phi_delta)
    ll = np.arange(nl)

    lmin=100.
    lmax=1500.
    dl=100.    
    nbins = np.ceil((lmax-lmin)/dl)
    llo = np.linspace(lmin, lmax, nbins)
    lcen = llo+0.5*dl
    delta_l_factor = 0.5
    cl_phi_delta_theory = np.load('/data/rkeisler/cl_phi_delta_mz0p65_sz0p25_b1p0.npy')
    cl_phi_delta_theory = cl_phi_delta_theory[0:len(cl_phi_delta)]
    y = cl_phi_delta/cl_phi_delta_theory
    #y = cl_phi_delta*ll**4.
    ybin = np.zeros(nbins)
    yerr = np.zeros(nbins)
    for i in np.arange(nbins):
        wh=np.where((ll>llo[i])&(ll<=(llo[i]+dl)))[0]
        ybin[i] = np.mean(y[wh])
        nl_tmp = 1.*len(wh)*delta_l_factor
        yerr[i] = np.std(y[wh])/np.sqrt(nl_tmp)

    pl.clf(); pl.errorbar(lcen, ybin, yerr=yerr)
    fs=17
    pl.xlabel('L',fontsize=fs)
    pl.ylabel(r'$C^{\phi \delta}_{\rm{meas}} / C^{\phi \delta}_{\rm{theory}}$',fontsize=fs)
    pl.title(r'Assuming dN/dz=N(0.65,0.25), b=1.0',fontsize=fs-1)
    pl.plot([0,1600],[0,0],'k--')
    pl.plot([0,1600],[1,1],'k--')
    pl.xlim(0,max(lcen)+dl/2.)
    print np.sqrt(np.sum((ybin/yerr)**2.))
    
    
    
    ipdb.set_trace()
