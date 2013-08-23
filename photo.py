import pickle
import numpy as np
import matplotlib.pylab as pl
import ipdb
from os import system
from os.path import exists
datapath = '/data/photo_data/'    


def make_irsa_table(filename=datapath+'atlas_0p05_to_1p0.csv', ngals=None, shuffle=True):
    # get SDSS data
    d = read_sdss(filename,keys2grab=['ra','dec'])
    id = d.keys()
    if shuffle:
        np.random.seed(3)
        np.random.shuffle(id)        
    ngals_total = len(id)
    if ngals==None: ngals=ngals_total
    tmp = filename.split('.')[0]
    tmp = tmp.split('_')
    zrange = tmp[-3]+'_to_'+tmp[-1]
    savename = datapath + 'irsa_table_'+zrange+'_%igals.txt'%ngals
    fileout = open(savename,'w')
    fileout.write('|cntr    |ra        |dec            |\n')
    for j in range(ngals):
        this_d = d[id[j]]
        fileout.write(' %7i   %9.5f  %9.5f\n'%(id[j], this_d['ra'], this_d['dec']))
    fileout.close()


def read_sdss(filename, keys2grab = ['redshift','objid','u','g','r','i','z','petroR50_r','petroR90_r','expAB_r']):
    dtype = {k:float for k in keys2grab}
    dtype['objid']=int
    print '...reading '+filename+'...'
    import csv
    f = open(filename,'r')
    reader = csv.reader(f)
    keys = reader.next()
    dkey = {}
    for itmp,key in enumerate(keys): dkey[key]=itmp
    ind2use = [dkey[key] for key in keys2grab]
    d = {}
    this_id = -1
    for line in reader:
        this_id += 1
        this_d = {}
        for key,ind in zip(keys2grab,ind2use): this_d[key] = map(dtype[key],[line[ind]])[0]
        d[this_id] = this_d
    return d



def read_wise_and_2mass(filename, quick=False):
    savename = filename.split('.')[0]+'.pkl'
    if quick: return pickle.load(open(savename,'r'))

    default_mag=50.0
    f=open(filename,'r')
    d = {}
    for line in f:
        if 'ra_01' in line:
            tmp = line.replace('|',' ')
            keys = tmp.split()
            for key in keys: d[key]=[]
        if 'double' in line:
            tmp = line.replace('|',' ')
            dtypes = tmp.split()
        if ('\\' in line) | ('|' in line): continue
        tmp = line.split()
        for ikey,thing in enumerate(tmp):
            key = keys[ikey]
            dtype = dtypes[ikey]
            if dtype=='char': d[key].append(thing)
            if dtype=='int': d[key].append(int(thing))
            if dtype=='double': 
                if thing=='null': thing=default_mag
                d[key].append(float(thing))
    f.close()

    # now create the output dictionary.  the keys are the cntr_01, and
    # the values are the four WISE magnitudes.  for each unique
    # cntr_01, find all matches and take the closest one.
    cntrs = np.array(d['cntr_01'])
    output = {}
    count=0
    nuniq = len(set(cntrs))
    for this_cntr in set(cntrs):
        count += 1
        if (count % 10000)==0: print count,nuniq
        wh=np.where(cntrs==this_cntr)[0]
        if len(wh)==1:
            output[this_cntr] = {'w1':d['w1mpro'][wh], 
                                 'w2':d['w2mpro'][wh],
                                 'w3':d['w3mpro'][wh], 
                                 'w4':d['w4mpro'][wh],
                                 'j':d['j_m_2mass'][wh],
                                 'h':d['h_m_2mass'][wh],
                                 'k':d['k_m_2mass'][wh]}
        else:
            c_ra = d['ra_01'][wh[0]]
            c_dec = d['dec_01'][wh[0]]
            cos_tmp = np.cos(c_dec*np.pi/180.)
            dmin=9e9
            for ind in wh:
                this_d = d['dist_x'][ind]
                #this_d = ((c_ra-d['ra'][ind])*cos_tmp)**2. + (c_dec-d['dec'][ind])**2.
                if this_d<dmin:
                    dmin=this_d
                    ind2use=ind
            output[this_cntr] = {'w1':d['w1mpro'][ind2use], 
                                 'w2':d['w2mpro'][ind2use],
                                 'w3':d['w3mpro'][ind2use], 
                                 'w4':d['w4mpro'][ind2use],
                                 'j':d['j_m_2mass'][ind2use],
                                 'h':d['h_m_2mass'][ind2use],
                                 'k':d['k_m_2mass'][ind2use]}

    # finally, replace the default null value with the minimum flux
    # seen per band.
    print '...fixing default...'
    kk=output.keys()
    bands = output[kk[0]].keys()
    default_wm_dict = {}
    for band in bands:
        tmp = np.array([vv[band] for vv in output.values()])
        true_max = np.max(tmp[np.where(tmp<default_mag)[0]])
        default_wm_dict[band]=true_max
        for k,v in output.iteritems():
            if (v[band]>=true_max): output[k][band]=true_max

    print '...saving...'
    pickle.dump((output, default_wm_dict), open(savename,'w'))
    return output, default_wm_dict


def main(filename=datapath+'atlas_0p05_to_1p0.csv', 
         doPlot=False, n_estimators=30, n_jobs=1, 
         surveys=['sdss','wise','2mass']):

    from sklearn.ensemble import ExtraTreesRegressor as rfRegressor
    from sklearn.cross_validation import train_test_split

    # get data
    d = read_sdss(filename)

    # get wise and 2mass data
    wm, default_wm_dict = read_wise_and_2mass(datapath+'wise_2mass_100000gals_0p05_to_1p0.txt',quick=True)

    # fold one dataset into another.  NB: need to be careful here, as
    # currently i'm effectively requiring a WISE match.
    d = fold_sdss_into_wm(d, wm)

    # create feature vector and labels.
    X = get_feautres(d, surveys=surveys)
    y = np.array([v['redshift'] for v in d.values()])

    #pickle.dump((X,y),open(datapath+'Xy.pkl','w'))

    # split data into train and test.
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    
    # define and fit a RF.  use cross-validation to determine hyper-parameters.
    rf = rfRegressor(n_jobs=n_jobs, n_estimators=n_estimators)

    print '...training...'
    rf.fit(X_train, y_train)
    print '...predicting...'
    predict_test = rf.predict(X_test)
    diff = predict_test-y_test
    this_sigma = np.std(diff)
    this_sigma68 = get_sigma68(diff)

    # summary vs redshift, and vs r-magnitude.
    zbin, zbias, zbias2, zsigma, zsigma68 = summary_vs_redshift(y_test, diff)
    #r_ind = 0
    #r_test = X_test[:,r_ind]
    #rbin, rbias, rbias2, rsigma, rsigma68 = summary_vs_r(r_test, diff)
    
    # if desired, plot some stuff.
    if doPlot: plot_summary(y_test, predict_test, y_train, diff, 
                            zbin, zbias, zbias2, zsigma, zsigma68)
#                            rbin, rbias, rbias2, rsigma, rsigma68)
    ipdb.set_trace()




def remove_bias(z_predict, z_true):
    z = np.linspace(np.min(z_predict), np.max(z_predict), 10)
    for i in np.arange(len(z)-1):
        wh=np.where((z_predict>=z[i]) & (z_predict<z[i+1]))[0]
        z_predict[wh] -= np.mean((z_predict-z_true)[wh])
    return z_predict

def get_sigma68(diff):
    return np.percentile(np.abs(diff),68.)

def plot_ngals_vs_match_radius():
    f = open(datapath+'num_sources_v_match_radius.txt','r')
    f.next()
    rad = []
    n = []
    for line in f:
        tmp = line.split()
        rad.append(float(tmp[0]))
        n.append(int(tmp[1]))
    pl.clf()
    pl.loglog(rad,n,'bo')
    pl.loglog(rad,n,'b')
    pl.loglog(rad, 550.*(np.array(rad)/15.)**2. + 990.,'r')
    pl.ylim(800,np.max(n))
    f.close()

def summary_vs_redshift(redshift, diff):
    zlo = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75])
    zhi = np.array([0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85])
    zbin = 0.5*(zlo+zhi)
    nbin = len(zbin)
    bias = np.zeros(nbin)
    bias2 = np.zeros(nbin)
    sigma = np.zeros(nbin)
    sigma68 = np.zeros(nbin)
    for i in np.arange(nbin):
        wh=np.where((redshift>=zlo[i]) & (redshift<zhi[i]))[0]
        nwh = len(wh)
        bias[i] = diff[wh].mean()
        bias2[i] = np.median(diff[wh])
        sigma[i] = diff[wh].std()
        sigma68[i] = get_sigma68(diff[wh])
    return zbin, bias, bias2, sigma, sigma68


def summary_vs_r(r_test, diff):    
    rlo = np.array([16,17,18,19,20])
    rhi = np.array([17,18,19,20,21])
    rbin = 0.5*(rlo+rhi)
    nrbin = len(rbin)
    rbias = np.zeros(nrbin)
    rbias2 = np.zeros(nrbin)
    rsigma = np.zeros(nrbin)
    rsigma68 = np.zeros(nrbin)
    for i in np.arange(nrbin):
        wh=np.where((r_test>=rlo[i]) & (r_test<rhi[i]))[0]
        nwh = len(wh)
        rbias[i] = diff[wh].mean()
        rbias2[i] = np.median(diff[wh])
        rsigma[i] = diff[wh].std()
        rsigma68[i] = get_sigma68(diff[wh])
    return rbin, rbias, rbias2, rsigma, rsigma68

def plot_summary(y_test, predict_test, y_train, diff, 
                 zbin, zbias, zbias2, zsigma, zsigma68):
#                 rbin, rbias, rbias2, rsigma, rsigma68):
    pl.ion()
    pl.figure(1, figsize=(15.5, 9.5))
    pl.clf()
    pl.subplot(321)
    pl.plot(y_test, predict_test,'.')
    pl.plot([0,2],[0,2],'k--',linewidth=5)
    pl.xlim(y_train.min(), y_train.max())
    pl.ylim(y_train.min(), y_train.max())
    pl.xlabel('specZ')
    pl.ylabel('photoZ')

    pl.subplot(322)
    pl.plot(y_test, diff,'.')
    pl.plot([0,2],[0,0],'k--',linewidth=5)
    pl.xlim(y_train.min(), y_train.max())
    pl.ylim(np.array([-1,1])*0.1)
    pl.xlabel('specZ')
    pl.ylabel('photoZ-specZ')

    pl.subplot(323)
    pl.plot(zbin, zbias,'ro')
    pl.plot(zbin, zbias2,'go')
    pl.legend(['mean','median'])
    pl.plot(zbin, zbias,'r')
    pl.plot(zbin, zbias2,'g')
    pl.plot([0,2],[0,0],'k--',linewidth=5)
    pl.xlim(y_train.min(), y_train.max())
    pl.ylim(np.array([-1,1])*0.1)
    pl.xlabel('specZ')
    pl.ylabel('bias')

    pl.subplot(325)
    pl.plot(zbin, zsigma,'ro')
    pl.plot(zbin, zsigma68,'go')
    pl.legend(['RMS','68%'])
    pl.plot(zbin, zsigma,'r')
    pl.plot(zbin, zsigma68,'g')
    pl.xlim(y_train.min(), y_train.max())
    pl.ylim(0.0, 0.1)
    pl.xlabel('specZ')
    pl.ylabel('sigma')
    
    if False:

        pl.subplot(324)
        pl.plot(rbin, rbias,'ro')
        pl.plot(rbin, rbias2,'go')
        pl.legend(['mean','median'])
        pl.plot(rbin, rbias,'r')
        pl.plot(rbin, rbias2,'g')
        pl.plot([0,30],[0,0],'k--',linewidth=5)
        pl.xlim(np.min(rbin), np.max(rbin))
        pl.ylim(np.array([-1,1])*0.05)
        pl.xlabel('r')
        pl.ylabel('bias')

        pl.subplot(326)
        pl.plot(rbin, rsigma,'ro')
        pl.plot(rbin, rsigma68,'go')
        pl.legend(['RMS','68%'])
        pl.plot(rbin, rsigma,'r')
        pl.plot(rbin, rsigma68,'g')
        pl.xlim(np.min(rbin), np.max(rbin))
        pl.ylim(0.0, 0.05)
        pl.xlabel('r')
        pl.ylabel('sigma')




def get_feautres(d, surveys=['sdss','wise','2mass']):
    ngals = len(d)
    X = None

    if 'sdss' in surveys:
        u = np.array([v['u'] for v in d.values()])
        g = np.array([v['g'] for v in d.values()])
        r = np.array([v['r'] for v in d.values()])
        i = np.array([v['i'] for v in d.values()])
        z = np.array([v['z'] for v in d.values()])
        petroR50_r = np.array([v['petroR50_r'] for v in d.values()])
        petroR90_r = np.array([v['petroR90_r'] for v in d.values()])
        expAB_r = np.array([v['expAB_r'] for v in d.values()])
        # EDIT FEATURES HERE
        tmp = np.vstack((r, u-g, g-r, r-i, i-z, petroR50_r, petroR90_r, expAB_r))
        if X==None: X=tmp
        else: X=np.vstack((X,tmp))

    if 'wise' in surveys:
        w1 = np.array([v['w1'] for v in d.values()])
        w2 = np.array([v['w2'] for v in d.values()])
        w3 = np.array([v['w3'] for v in d.values()])
        w4 = np.array([v['w4'] for v in d.values()])
        # EDIT FEATURES HERE
        tmp = np.vstack((w2, w1-w2, w2-w3, w4))
        if X==None: X=tmp
        else: X=np.vstack((X,tmp))

    if '2mass' in surveys:
        j = np.array([v['j'] for v in d.values()])
        h = np.array([v['h'] for v in d.values()])
        k = np.array([v['k'] for v in d.values()])
        # EDIT FEATURES HERE
        tmp = np.vstack((j, j-h, h-k))
        if X==None: X=tmp
        else: X=np.vstack((X,tmp))

    if ('sdss' in surveys) & ('wise' in surveys):
        X=np.vstack((X, i-w2))

    if ('sdss' in surveys) & ('2mass' in surveys):
        X=np.vstack((X, i-j))

    if ('wise' in surveys) & ('2mass' in surveys):
        X=np.vstack((X, j-w2))

    return X.T

def fold_wm_into_sdss(d, wm, default_wm_dict):
    for k in d.keys():
        if k in wm: 
            d[k] = dict(d[k].items()+wm[k].items())
        else:
            d[k] = dict(d[k].items()+default_wm_dict.items())
    return d


def fold_sdss_into_wm(d, wm):
    for k in wm.keys(): wm[k] = dict(d[k].items()+wm[k].items())
    return wm
    
