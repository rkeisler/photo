import numpy as np
import matplotlib.pylab as pl
pl.ion()
import pdb
import pickle

wise_datapath='/data/rkeisler/wise/'


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
    #from matplotlib.pyplot import tight_layout
    #tight_layout()
    #pl.legend(['WISE','BOSS ('+sample.upper()+')'])


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
    d23 = {'color':w2-w4, 'name':'w2-w4', 'range':[-2,10]}
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


def boil_down_wise(filename=wise_datapath+'wise-allsky-cat-part44', savepath=wise_datapath):
    savename = savepath+filename.split('/')[-1]
    f=open(filename)
    default_mag = 30.
    w_ind = wise_cat_ind()
    np_ind = np_array_ind()
    nsave = len(np_array_ind())
    count = -1
    X = []
    for line in f:
        count+=1
        if (count % 100000)==0: print count
        tmp = line.split('|')
        thisX = np.zeros(nsave)
        for k in w_ind:
            valtmp = tmp[w_ind[k]]
            if valtmp=='': valtmp=default_mag
            else: valtmp=float(valtmp)
            thisX[np_ind[k]]=valtmp
        X.append(thisX)        
    np.save(savename, np.vstack(X))


def wise_cat_ind():
    return {'ra':1, 'dec':2,
            'w1':16, 'w2':20, 'w3':24, 'w4':28,
            'j':273, 'h':275, 'k':277,
            'w1sn':18, 'w1cov':54}


def np_array_ind():
    tmp={'ra':0, 'dec':1,
         'w1':2, 'w2':3, 'w3':4, 'w4':5,
         'j':6, 'h':7, 'k':8,
         'w1sn':9, 'w1cov':10}
    if len(tmp)!=max(tmp.values())+1:
        print 'huh?'
        pdb.set_trace()
    return tmp


def process_wise_dec_strip(part):
    from os import system
    print '...downloading...'
    this_address = 'http://irsadist.ipac.caltech.edu/wise-allsky/wise-allsky-cat-part%02d.bz2'%part
    print this_address
    system('curl -O '+this_address)
    print '...unzipping...'
    system('bunzip2 '+savename)
    print '...reducing...'
    boil_down_wise(filename=savename, savepath='/data/wise/allsky/')
    print '...removing raw data...'
    system('rm '+savename)

def process_many(imin_incl=1, imax_excl=51):
    for i in range(imin(incl, imax_excl): process_wise_dec_strip(i)


    
