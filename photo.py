import pickle
import numpy as np
import matplotlib.pylab as pl
import ipdb
from os import system
from os.path import exists
datapath = '/data/photo_data/'    


def make_irsa_table(filename=datapath+'atlas_0p05_to_1p0.csv', ngals=None, shuffle=True):
    # get data
    d = read_sdss(filename)
    ngals_total = len(d['ra'])
    if ngals==None: ngals=ngals_total
    # grab a random (sub)set
    ind = np.arange(ngals_total, dtype=int)
    if shuffle:
        np.random.seed(3)
        np.random.shuffle(ind)
    ind = ind[0:ngals]
    tmp = filename.split('.')[0]
    tmp = tmp.split('_')
    zrange = tmp[1]+'_to_'+tmp[3]
    savename = datapath + 'irsa_table_'+zrange+'_%igals.txt'%ngals
    fileout = open(savename,'w')
    fileout.write('|cntr     |ra        |dec            |\n')
    for j in ind:
        fileout.write(' %7i   %9.5f  %9.5f\n'%(j,d['ra'][j],d['dec'][j]))
    fileout.close()


def read_sdss(filename):
    print '...reading '+filename+'...'
    import csv
    f = open(filename,'r')
    reader = csv.reader(f)
    keys = reader.next()
    d = {key:[] for key in keys}
    for line in reader:
        for key,val in zip(keys,line):
            d[key].append(np.float(val))
    for k,v in d.iteritems(): d[k] = np.array(v,dtype='float')
    id = np.arange(len(d[d.keys()[0]]), dtype=int)
    d['id'] = id
    return d

   
def read_wise(filename, quick=True):
    savename = filename.split('.')[0]+'.pkl'
    if quick: return pickle.load(open(savename,'r'))
    f=open(filename,'r')
    d = {}
    for line in f:
        if 'ra_01' in line:
            tmp = line.replace('|','')
            keys = tmp.split()
            for key in keys: d[key]=[]
        if 'double' in line:
            tmp = line.replace('|','')
            dtypes = tmp.split()
        if ('\\' in line) | ('|' in line): continue
        tmp = line.split()
        for ikey,thing in enumerate(tmp):
            key = keys[ikey]
            dtype = dtypes[ikey]
            if dtype=='char': d[key].append(thing)
            if dtype=='int': d[key].append(int(thing))
            if dtype=='double': 
                if thing=='null': thing=20.
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
                                 'w4':d['w4mpro'][wh]}
        else:
            c_ra = d['ra_01'][wh[0]]
            c_dec = d['dec_01'][wh[0]]
            cos_tmp = np.cos(c_dec*np.pi/180.)
            dmin=9e9
            for ind in wh:
                this_d = ((c_ra-d['ra'][ind])*cos_tmp)**2. + (c_dec-d['dec'][ind])**2.
                if this_d<dmin:
                    dmin=this_d
                    ind2use=ind
            output[this_cntr] = {'w1':d['w1mpro'][ind2use], 
                                 'w2':d['w2mpro'][ind2use],
                                 'w3':d['w3mpro'][ind2use], 
                                 'w4':d['w4mpro'][ind2use]}
    print '...saving...'
    pickle.dump(output,open(savename,'w'))
    return output
                    


def read_2mass(filename, quick=True):
    savename = filename.split('.')[0]+'.pkl'
    if quick: return pickle.load(open(savename,'r'))

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
                if thing=='null': thing=20.
                if thing=='-': thing=-1.
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
        wh=np.where(cntrs==this_cntr)[0][0]
        output[this_cntr] = {'j':d['j_m'][wh], 
                             'h':d['h_m'][wh],
                             'k':d['k_m'][wh]}
    print '...saving...'
    pickle.dump(output,open(savename,'w'))
    return output



def main(filename=datapath+'atlas_0p05_to_1p0.csv', 
         doPlot=False, n_estimators=30, n_jobs=1):

    from sklearn.ensemble import ExtraTreesRegressor as rfRegressor
    from sklearn.cross_validation import train_test_split

    # get data
    d = read_sdss(filename)

    # get wise data
    wise=read_wise(datapath+'wise_100000gals_0p05_to_1p0.txt',quick=True)
    dd={}
    for k in d:
        dd[k] = d[k][wise.keys()]
    for k in wise.values()[0]:
        dd[k] = np.array([thing[k] for thing in wise.values()])
    d=dd


    # get 2mass data
    #mass=read_2mass(datapath+'2mass_100000gals_0p05_to_1p0.txt',quick=True)
    #keys2use = ['j','h','k']
    #defaults = {'j':20., 'h':20., 'k':20.}
    #for key in keys2use: d[key]=np.zeros(len(d['redshift']))+defaults[key]
    #for key in keys2use:
    #    for id,v in mass.iteritems(): d[key][id] = v[key]
    #ipdb.set_trace()

    # create feature vector and labels.
    X = np.vstack((
            d['u'],
            d['g'],
            d['r'],
            d['i'],
            d['z'],
            d['u']-d['g'], 
            d['g']-d['r'], 
            d['r']-d['i'], 
            d['i']-d['z'], 
            d['g']-d['i'], 
            d['g']-d['z'], 
            d['expAB_u'], 
            d['expAB_g'], 
            d['expAB_r'], 
            d['expAB_i'], 
            d['expAB_z'], 
            d['petroR50_u'], 
            d['petroR90_u'], 
            d['petroR50_g'], 
            d['petroR90_g'], 
            d['petroR50_r'], 
            d['petroR90_r'], 
            d['petroR50_i'], 
            d['petroR90_i'], 
            d['petroR50_z'], 
            d['petroR90_z'], 
            d['petroR50_u']/d['petroR90_u'],
            d['petroR50_g']/d['petroR90_g'],
            d['petroR50_r']/d['petroR90_r'],
            d['petroR50_i']/d['petroR90_i'],
            d['petroR50_z']/d['petroR90_z'],
            d['w1']-d['u'],
            d['w1']-d['g'],
            d['w1']-d['r'],
            d['w1']-d['i'],
            d['w1']-d['z'],
            d['w1']-d['w2'],
            d['w1']-d['w3'],
            d['w1']-d['w4'],
            d['w2']-d['u'],
            d['w2']-d['g'],
            d['w2']-d['r'],
            d['w2']-d['i'],
            d['w2']-d['z'],
            d['w2']-d['w1'],
            d['w2']-d['w3'],
            d['w2']-d['w4'],
            d['w3']-d['u'],
            d['w3']-d['g'],
            d['w3']-d['r'],
            d['w3']-d['i'],
            d['w3']-d['z'],
            d['w3']-d['w2'],
            d['w3']-d['w1'],
            d['w3']-d['w4'],
            d['w4']-d['u'],
            d['w4']-d['g'],
            d['w4']-d['r'],
            d['w4']-d['i'],
            d['w4']-d['z'],
            d['w4']-d['w2'],
            d['w4']-d['w3'],
            d['w4']-d['w1'],
            d['w2'],
            d['w1'],
            )).T
    y = d['redshift']

    pickle.dump((X,y),open(datapath+'Xy.pkl','w'))

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
    r_ind = 0
    r_test = X_test[:,r_ind]
    rbin, rbias, rbias2, rsigma, rsigma68 = summary_vs_r(r_test, diff)
    
    # if desired, plot some stuff.
    if doPlot: plot_summary(y_test, predict_test, y_train, diff, 
                            zbin, zbias, zbias2, zsigma, zsigma68, 
                            rbin, rbias, rbias2, rsigma, rsigma68)
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
                 zbin, zbias, zbias2, zsigma, zsigma68, 
                 rbin, rbias, rbias2, rsigma, rsigma68):
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
    pl.ylim(np.array([-1,1])*0.05)
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
    pl.ylim(np.array([-1,1])*0.05)
    pl.xlabel('specZ')
    pl.ylabel('bias')


    pl.subplot(325)
    pl.plot(zbin, zsigma,'ro')
    pl.plot(zbin, zsigma68,'go')
    pl.legend(['RMS','68%'])
    pl.plot(zbin, zsigma,'r')
    pl.plot(zbin, zsigma68,'g')
    pl.xlim(y_train.min(), y_train.max())
    pl.ylim(0.0, 0.05)
    pl.xlabel('specZ')
    pl.ylabel('sigma')



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

