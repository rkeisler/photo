import numpy as np
import matplotlib.pylab as pl
pl.ion()
import pdb

def parse_cmass():
    import pyfits
    x=pyfits.open('bosstile-final-collated-boss2-boss32.fits')
    tbdata = x[1].data 
    ra=tbdata.field(12)
    dec=tbdata.field(13)
    wh_cmass=np.where(((tbdata.field(44)&2)>0) & (tbdata.field(44)>0))[0]
    wh_cmass_decrange=np.where(((tbdata.field(44)&2)>0) & (tbdata.field(44)>0) & (dec>47.385300) & (dec<=50.678000))[0]
    pdb.set_trace()


def boil_down_wise2(filename='wise-allsky-cat-part44'):
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
    np.save(filename, np.vstack(X))


def wise_cat_ind():
    return {'ra':1, 'dec':2,
            'w1':16, 'w2':20, 'w3':24, 'w4':28,
            'j':273, 'h':275, 'k':277}


def np_array_ind():
    tmp={'ra':0, 'dec':1,
         'w1':2, 'w2':3, 'w3':4, 'w4':5,
         'j':6, 'h':7, 'k':8}
    if len(tmp)!=max(tmp.values())+1:
        print 'huh?'
        pdb.set_trace()
    return tmp

