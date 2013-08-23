import numpy as np
import matplotlib.pylab as pl
import ipdb
from os import system
from os.path import exists
#import sdsspy
#pl.ion()

def cutouts_to_batches():
    import pickle
    #import cPickle as pickle
    from glob import glob

    nbatches = 100
    ngals_per_batch = 500

    system('rm batches/data_batch*')
    objid_to_redshift = get_objid_redshift_dict()

    files = glob('cutouts/*npz')

    # read the first file to get shape info.
    tmp = np.load(files[0])
    this_data = tmp['arr_0']
    
    nbands, nx, ny = this_data.shape
    nfeatures = nbands*nx*ny
   
    global_counter = -1
    sum_data = np.zeros(nfeatures)
    for ibatch in np.arange(nbatches):
        data = np.zeros((nfeatures, ngals_per_batch),dtype=np.float16)
        labels = np.zeros(ngals_per_batch)
        for igal in np.arange(ngals_per_batch):
            print global_counter
            global_counter += 1
            this_file = files[global_counter]
            this_objid = this_file.split('/')[-1].split('.')[0]
            this_redshift = np.float(objid_to_redshift[this_objid])
            labels[igal] = this_redshift

            tmp = np.load(this_file)
            this_data = tmp['arr_0'].ravel()
            data[:,igal] = this_data
            sum_data += this_data
        print '...dumping...'
        output = {'data':data, 'labels':labels}
        pickle.dump(output, open('batches/data_batch_%i'%ibatch, 'w'))

    mean_data = 1.*sum_data/(global_counter+1.)
    meta = {'num_cases_per_batch':ngals_per_batch, 
            'num_vis':nfeatures, 
            'data_mean':mean_data}
    pickle.dump(meta, open('batches/batches.meta', 'w'))


        
def get_objid_redshift_dict(filename='/data/atlas_0p05_to_1p0.csv'):
    f = open(filename,'r')
    line = f.readline()
    header = line.rstrip().split(',')
    i=-1
    output = {}
    while True:
        i += 1
        if ((i%100000)==0): print i
        line = f.readline()
        if not line:
            break
        tmp = line.rstrip().split(',')
        dic_tmp = {thing[0]:thing[1] for thing in zip(header,tmp)}
        output[dic_tmp['objid']]=dic_tmp['redshift']
    return output


def process_many_in_serial(filename, imin_inclusive, imax_exclusive, nlines=None):
    # get the indices of the lines to process.
    ind, nlines = ind_to_process(filename, imin_inclusive, imax_exclusive, nlines=nlines)
    f = open(filename,'r')
    line = f.readline()
    header = line.rstrip().split(',')
    i=-1
    while True:
        i += 1
        if ((i%100000)==0): print i,nlines,1.*i/nlines
        line = f.readline()
        if not line:
            break
        if not(i in ind): continue
        tmp = line.rstrip().split(',')
        dic_tmp = {thing[0]:thing[1] for thing in zip(header,tmp)}
        result = process_one(dic_tmp)
        if (result==-1): break

    f.close()

def process_one(dic):
    this_address = address_from_dict(dic)
    savename = '/data/cutouts/'+dic['objid']
    #if exists(savename): return -1
    

    print '...downloading...'
    print this_address
    system('curl -O '+this_address)
    atlas_filename = filename_from_address(this_address)

    print '...processing atlas...'
    process_atlas(atlas_filename, dic, savename=savename)

    print '...removing atlas data...'
    system('rm '+atlas_filename)
    return 0
    

def process_atlas(atlas_filename, dic, nside=64, verbose=False, savename=None):
    # read the atlas images
    obj = int(dic['obj'])
    if verbose: print 'obj',obj
    d=sdsspy.atlas.read_atlas(atlas_filename,obj,trim=False)
    
    # define the center of the image to the 
    # brightest pixel in the r+g image.
    rg=d['images'][2]+d['images'][3]
    xc,yc = np.unravel_index(rg.argmax(), rg.shape)

    nbands = 5
    output = np.zeros((nbands, nside, nside))
    cal_keys = ['nMgyPerCount_u','nMgyPerCount_g','nMgyPerCount_r','nMgyPerCount_i','nMgyPerCount_z']
    for i in np.arange(nbands):
        img = d['images'][i]
        if verbose: print img.shape

        # get the cutout
        nx, ny = img.shape
        xmin_incl = np.max([xc-nside/2, 0])
        xmax_excl = np.min([xc+nside/2, nx])
        ymin_incl = np.max([yc-nside/2, 0])
        ymax_excl = np.min([yc+nside/2, ny])
        cutout = 1.*img[xmin_incl:xmax_excl, ymin_incl:ymax_excl]

        # subtract the soft offset from the cutout.
        if verbose: print 'soft_bias: ',d['SOFT_BIAS']
        if verbose: print 'min/max cutout raw: ',np.min(cutout), np.max(cutout)
        cutout -= np.float(d['SOFT_BIAS'])
        
        # calibrate from counts to nanomaggies.
        this_nMgyPerCount = np.float(dic[cal_keys[i]])
        if verbose: print 'this_nMgyPerCount: ',this_nMgyPerCount
        cutout *= this_nMgyPerCount

        dx_lo = xc-xmin_incl
        dx_hi = xmax_excl-xc
        dy_lo = yc-ymin_incl
        dy_hi = ymax_excl-yc

        # put the cutout in the output array.
        output[i, nside/2-dx_lo:nside/2+dx_hi, nside/2-dy_lo:nside/2+dy_hi] = cutout
        if verbose: print cal_keys[i]
        if verbose: print output[i,:,:]


    #np.save(savename, output)
    np.savez_compressed(savename, output)



def address_from_dict(dic):
    this_run = int(dic['run'])
    this_camcol = int(dic['camcol'])
    this_field = int(dic['field'])
    this_address = 'http://data.sdss3.org/sas/dr9/boss/photo/redux/301/%i/objcs/%i/fpAtlas-%06i-%i-%04i.fit'%(this_run, this_camcol, this_run, this_camcol, this_field)
    return this_address


def ind_to_process(filename, imin_inclusive, imax_exclusive, nlines=None):

    # first get the number of lines in this file.
    if nlines==None:
        nlines = get_nlines(filename)
        # we don't want to count the header as a line to process.
        nlines -= 1
    
    # now create a randomly shuffled list of indices.
    # THE SEED (3) MUST ALWAYS BE THE SAME.
    ind = np.arange(nlines)
    np.random.seed(3)
    np.random.shuffle(ind)
    return ind[imin_inclusive:imax_exclusive], nlines



def get_nlines(filename):
    import subprocess
    return int(subprocess.check_output('wc -l '+filename,shell=True).split()[0])


def try_downloads(ndownload=10):
    #d = read_csv('atlas_top1e4_0p05_to0p1.csv')
    d = read_csv('atlas_top1e3_0p05_to0p1.csv')
    #d = read_csv('atlas_short.csv')
    run = d['run']
    camcol = d['camcol']
    obj = d['obj']
    objid = d['objid']
    field = d['field']
    nrow = len(obj)
    address = []
    for i in np.arange(nrow):
        this_address = 'http://data.sdss3.org/sas/dr9/boss/photo/redux/301/%i/objcs/%i/fpAtlas-%06i-%i-%04i.fit'%(run[i], camcol[i], run[i], camcol[i], field[i])
        #print this_address
        address.append(this_address)

    address = np.array(address)
    download_random(address,ndownload=ndownload)


def download_random(addresses, ndownload=1):
    ind = np.arange(len(addresses),dtype=int)
    np.random.shuffle(ind)
    these_addresses = addresses[ind[0:ndownload]]
    for adr in these_addresses:
        system('curl -O '+adr)
        #system('rm '+filename_from_address(adr))

def filename_from_address(address):
    return address.split('/')[-1]
    


def main(filename='atlas_0p05_to_1p0.csv'):
    from sklearn.grid_search import GridSearchCV
    from PyWiseRF import WiseRF as rfRegressor
    #from sklearn.ensemble import RandomForestRegressor as rfRegressor
    from sklearn.cross_validation import train_test_split
    from sklearn.svm import SVR
    from sklearn.tree import DecisionTreeRegressor
    from sklearn.linear_model import LogisticRegression

    # get data
    d = read_csv(filename)
    

    # normalize data
    # not necessary for random forest, but we might try non-RF regressors too.
    #for thing in [u,g,r,i,z]:
    #    thing -= np.mean(thing)
    #    thing /= np.std(thing)

    # create feature vector and labels.
    #X = np.vstack((u,g,r,i,z)).T
    # create features
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
            )).T
    y = d['redshift']

    # split data into train and test.
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    # furthermore, split train into A and B
    X_trainA, X_trainB, y_trainA, y_trainB = train_test_split(X_train, y_train, test_size=0.5)
    
    # define and fit a RF.  use cross-validation to determine hyper-parameters.
    n_estimators = 30

    print '...training...'
    if True:
        depths = np.arange(1,50,10)
        ndepths = len(depths)
        sigma = np.zeros(ndepths)
        sigma68 = np.zeros(ndepths)
        for i,depth in enumerate(depths):
        #grf = GridSearchCV(rf, dict(max_depth=depths), cv=3).fit(X_train, y_train

            rf = rfRegressor(n_jobs=1, n_estimators=n_estimators, max_depth=depth)
            grf = rf.fit(X_train, y_train)
            predict_test = grf.predict(X_test)
            diff = predict_test-y_test
            this_sigma = np.std(diff)
            this_sigma68 = get_sigma68(diff)
            sigma[i] = this_sigma
            sigma68[i] = this_sigma68
            print i,depth,this_sigma, this_sigma68
  

            pl.figure(1)
            pl.clf()
            pl.plot(depths, sigma, 'b')
            pl.plot(depths, sigma68, 'g')
            
    else:
        rf = rfRegressor(n_jobs=3, n_estimators=n_estimators, min_samples_split=2)
        rf.fit(X_trainA, y_trainA)
        rf_predict_trainB = rf.predict(X_trainB)

        print '...training second thing...'
        rfB = rfRegressor(n_jobs=3, n_estimators=n_estimators, min_samples_split=2)
        metaX_trainB = np.hstack((X_trainB, rf_predict_trainB[:,np.newaxis]))
        rfB.fit(metaX_trainB, y_trainB)


        #sv = SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.1, gamma=0.0, kernel='rbf', max_iter=-1, probability=False,shrinking=True, tol=0.001, verbose=False)
        #sv = rfRegressor(n_jobs=3, n_estimators=n_estimators, min_samples_split=2)
        #sv = DecisionTreeRegressor()
        #sv.fit(rf_predict_train[:,np.newaxis], y_train)
        #sv = LogisticRegression(penalty='l2', dual=False, tol=0.0001, C=1.0, fit_intercept=True, intercept_scaling=1, class_weight=None, random_state=None)
        #pdb.set_trace()
        #sv.fit(rf_predict_train[:,np.newaxis], y_train)
        
        
        rf_predict_test = rf.predict(X_test)
        if True:
            metaX_test = np.hstack((X_test, rf_predict_test[:,np.newaxis]))
            predict_test = rfB.predict(metaX_test)
        else:
            #predict_test = rf_predict_test
            #predict_test = remove_bias(predict_test, y_test)
            fit = np.polyfit(rf_predict_train, y_train, 2)
            predict_test = np.polyval(fit, rf_predict_test)


        diff = predict_test-y_test


        this_sigma = np.std(diff)
        this_sigma68 = get_sigma68(diff)
        print this_sigma, this_sigma68



    #zlo = np.array([0.05, 0.15, 0.25, 0.35, 0.45])
    #zhi = np.array([0.15, 0.25, 0.35, 0.45, 0.55])
    zlo = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75])
    zhi = np.array([0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85])
    zbin = 0.5*(zlo+zhi)
    nbin = len(zbin)
    bias = np.zeros(nbin)
    bias2 = np.zeros(nbin)
    sigma = np.zeros(nbin)
    sigma68 = np.zeros(nbin)
    for i in np.arange(nbin):
        wh=np.where((y_test>=zlo[i]) & (y_test<zhi[i]))[0]
        nwh = len(wh)
        bias[i] = diff[wh].mean()
        bias2[i] = np.median(diff[wh])
        sigma[i] = diff[wh].std()
        sigma68[i] = get_sigma68(diff[wh])



    
    rlo = np.array([14,15,16,17,18,19,20])
    rhi = np.array([15,16,17,18,19,20,21])
    rbin = 0.5*(rlo+rhi)
    nrbin = len(rbin)
    rbias = np.zeros(nrbin)
    rbias2 = np.zeros(nrbin)
    rsigma = np.zeros(nrbin)
    rsigma68 = np.zeros(nrbin)
    r_ind = 0
    r_test = X_test[:,r_ind]
    for i in np.arange(nrbin):
        wh=np.where((r_test>=rlo[i]) & (r_test<rhi[i]))[0]
        nwh = len(wh)
        rbias[i] = diff[wh].mean()
        rbias2[i] = np.median(diff[wh])
        rsigma[i] = diff[wh].std()
        rsigma68[i] = get_sigma68(diff[wh])
    


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
    pl.plot(zbin, bias,'ro')
    pl.plot(zbin, bias2,'go')
    pl.legend(['mean','median'])
    pl.plot(zbin, bias,'r')
    pl.plot(zbin, bias2,'g')
    pl.plot([0,2],[0,0],'k--',linewidth=5)
    pl.xlim(y_train.min(), y_train.max())
    pl.ylim(np.array([-1,1])*0.05)
    pl.xlabel('specZ')
    pl.ylabel('bias')


    pl.subplot(325)
    pl.plot(zbin, sigma,'ro')
    pl.plot(zbin, sigma68,'go')
    pl.legend(['RMS','68%'])
    pl.plot(zbin, sigma,'r')
    pl.plot(zbin, sigma68,'g')
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


    pdb.set_trace()


def remove_bias(z_predict, z_true):

    z = np.linspace(np.min(z_predict), np.max(z_predict), 10)
    for i in np.arange(len(z)-1):
        wh=np.where((z_predict>=z[i]) & (z_predict<z[i+1]))[0]
        z_predict[wh] -= np.mean((z_predict-z_true)[wh])
    return z_predict


def get_sigma68(diff):
    return np.percentile(np.abs(diff),68.)



def read_csv(filename):
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




def sparse_filter_images(ncodes=100, nbatch=10):
    import pickle
    import sf

    for i in np.arange(nbatch): 
        print i
        tmp = pickle.load(open('batches/data_batch_%i'%i, 'r'))
        nsamples = tmp['data'].shape[1]
        xtmp = tmp['data'].reshape(4,32,32,nsamples)[0:3,:,:].reshape(3*32*32,nsamples)
        if i==0:
            X=xtmp
            y=tmp['labels']
        else:
            X=np.hstack((X,xtmp))
            y=np.hstack((y,tmp['labels']))
    X = X.T
    w = sf.sparse_filter(X, n_codes=ncodes, activation='softabs', x0=None)
    pickle.dump(w,open('just_in_case.pkl','w'))
    pickle.dump(w, open('codes_3band32x32_to_%icodes.pkl'%ncodes,'w'))
    


def fit_sparse_filter(ncodes=100, nbatch=10):
    import pickle
    import sf
    w = pickle.load(open('codes_3band32x32_to_%icodes.pkl'%ncodes,'r'))
    for i in np.arange(nbatch): 
        print i
        tmp = pickle.load(open('batches/data_batch_%i'%i, 'r'))
        nsamples = tmp['data'].shape[1]
        xtmp = tmp['data'].reshape(4,32,32,nsamples)[0:3,:,:].reshape(3*32*32,nsamples)
        if i==0:
            X=xtmp
            y=tmp['labels']
        else:
            X=np.hstack((X,xtmp))
            y=np.hstack((y,tmp['labels']))
    X = X.T
    X = np.dot(X, w)
    
    from PyWiseRF import WiseRF as rfRegressor
    from sklearn.cross_validation import train_test_split
    from sklearn.svm import SVR
    rf = rfRegressor(n_jobs=3, n_estimators=100, compute_importances=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    rf.fit(X_train, y_train)
    #y_predict = rf.predict(X_test)
    #diff = y_predict-y_test

    wh=np.where(rf.feature_importances_>0.01)[0]
    from sklearn.grid_search import GridSearchCV

    svm = SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.1, gamma=0.0, kernel='rbf', max_iter=-1, probability=False,shrinking=True, tol=0.001, verbose=False)        
    svm.fit(X_train[:,wh], y_train)
    best_svm = GridSearchCV(svm, dict(epsilon=10.**np.arange(-4,3)), cv=3).fit(X_train, y_train)
    y_predict = svm.predict(X_test[:,wh])
    diff = y_predict-y_test
    print np.std(diff), get_sigma68(diff)


        
        
        
def mag_vs_z(nbatch=10):
    import pickle

    for i in np.arange(nbatch): 
        print i
        tmp = pickle.load(open('batches/data_batch_%i'%i, 'r'))
        nsamples = tmp['data'].shape[1]
        xtmp = tmp['data'].reshape(4,32,32,nsamples)  #[0:3,:,:].reshape(3*32*32,nsamples)
        xtmp = np.exp(xtmp)-1.5
        #for iband in np.arange(xtmp.shape[0]):
            #mag[iband,
        
        if i==0:
            X=xtmp
            y=tmp['labels']
        else:
            X=np.hstack((X,xtmp))
            y=np.hstack((y,tmp['labels']))

    ipdb.set_trace()




def sparse_filter_summary(ncodes=100,filename='atlas_0p05_to_1p0.csv'):
    import pickle
    import sf


    # get data
    d = read_csv(filename)

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
            )).T
    y = d['redshift']

    # normalize features
    X -= (X.mean(axis=0))
    for i in np.arange(X.shape[1]):
        X[:,i] /= get_sigma68(X[:,i])


    #ind = np.arange(X.shape[0])
    #np.random.seed(5)
    #np.random.shuffle(ind)
    #ind = ind[0:100000]
    #X = X[ind,:]

    #X = X.T

    w = sf.sparse_filter(X, n_codes=ncodes, activation='softabs', x0=None)
    pickle.dump(w,open('just_in_case.pkl','w'))
    pickle.dump(w, open('summary_to_%icodes.pkl'%ncodes,'w'))




def fit_sparse_filter_summary(filename='atlas_0p05_to_1p0.csv'):

    import pickle
    import sf


    # get data
    d = read_csv(filename)

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
            )).T
    y = d['redshift']

    # normalize features
    X -= (X.mean(axis=0))
    for i in np.arange(X.shape[1]):
        X[:,i] /= get_sigma68(X[:,i])

    w = np.load('sf_tmp.npy')

    w = w.reshape(50,len(w)/50)
    #X = np.dot(X, w.T)

    
    from PyWiseRF import WiseRF as rfRegressor
    from sklearn.cross_validation import train_test_split
    from sklearn.svm import SVR
    rf = rfRegressor(n_jobs=3, n_estimators=100, compute_importances=True)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    rf.fit(X_train, y_train)
    y_predict = rf.predict(X_test)
    diff = y_predict-y_test
    print np.std(diff), get_sigma68(diff)

    ipdb.set_trace()


def prep_for_NN(filename='atlas_0p05_to_1p0.csv', full='False',nbatches=6):

    import pickle
    import sf

    # get data
    d = read_csv(filename)

    # create feature vector and labels.
    if full:
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
                )).T
    else:
        X = np.vstack((
                d['u'],
                d['g'],
                d['r'],
                d['i'],
                d['z'])).T
    y = d['redshift']

    # normalize features
    for i in np.arange(X.shape[1]):
        
        med = np.median(X[:,i])
        s68 = get_sigma68(X[:,i])
        whbad = np.where(np.abs(X[:,i]-med)>(5.*s68))[0]
        X[whbad,i] = med
        X[:,i] -= med
        X[:,i] /= s68


    import pickle
    if full: savdir='full/' 
    else: savdir='lite/'
    from os import system
    system('rm '+savdir+'/*')
    ngals_total = X.shape[0]
    nfeatures = X.shape[1]
    ngals_per_batch = np.floor(1.*ngals_total/nbatches)

    ind = np.arange(ngals_total)
    np.random.seed(3)
    np.random.shuffle(ind)

    mean_data = X.mean(axis=0)[:,np.newaxis]
    meta = {'num_cases_per_batch':ngals_per_batch, 
            'num_vis':nfeatures, 
            'data_mean':mean_data}
    pickle.dump(meta, open(savdir+'batches.meta', 'w'))    

    for ibatch in np.arange(nbatches):
        print ibatch
        ind_tmp = ind[ibatch*ngals_per_batch: (ibatch+1)*ngals_per_batch]
        pickle.dump({'data':X.T[:,ind_tmp], 'labels':y[ind_tmp]}, open(savdir+'data_batch_%i'%ibatch,'w'))



def make_irsa_table(filename='atlas_0p05_to_1p0.csv', ngals=None, shuffle=True):


    # get data
    d = read_csv(filename)
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
    savename = 'irsa_table_'+zrange+'_%igals.txt'%ngals
    fileout = open(savename,'w')
    fileout.write('|cntr     |ra        |dec            |\n')
    for j in ind:
        fileout.write(' %7i   %9.5f  %9.5f\n'%(j,d['ra'][j],d['dec'][j]))
    fileout.close()


def plot_ngals_vs_match_radius():
    f = open('num_sources_v_match_radius.txt','r')
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

    
def read_wise(filename='wise_100000gals_0p05_to_1p0.txt', quick=True):
    import pickle
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
                    



def mainwise(filename='atlas_0p05_to_1p0.csv'):
    from sklearn.grid_search import GridSearchCV
    #from PyWiseRF import WiseRF as rfRegressor
    from sklearn.ensemble import RandomForestRegressor as rfRegressor
    from sklearn.cross_validation import train_test_split
    from sklearn.svm import SVR
    from sklearn.tree import DecisionTreeRegressor
    from sklearn.linear_model import LogisticRegression

    # get data
    d = read_csv(filename)

    # get wise data
    wise=read_wise(filename='wise_100000gals_0p05_to_1p0.txt',quick=True)
    dd={}
    for k in d:
        dd[k] = d[k][wise.keys()]
    for k in wise.values()[0]:
        dd[k] = np.array([thing[k] for thing in wise.values()])
    d=dd


    # get 2mass data
    mass=read_2mass(filename='2mass_100000gals_0p05_to_1p0.txt',quick=True)
    keys2use = ['j','h','k']
    defaults = {'j':20., 'h':20., 'k':20.}
    for key in keys2use: d[key]=np.zeros(len(d['redshift']))+defaults[key]
    for key in keys2use:
        for id,v in mass.iteritems(): d[key][id] = v[key]
    ipdb.set_trace()


    # normalize data
    # not necessary for random forest, but we might try non-RF regressors too.
    #for thing in [u,g,r,i,z]:
    #    thing -= np.mean(thing)
    #    thing /= np.std(thing)

    # create feature vector and labels.
    #X = np.vstack((u,g,r,i,z)).T
    # create features
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

    #pickle.dump((X,y),open('Xy.pkl','w'))


    # split data into train and test.
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    # furthermore, split train into A and B
    X_trainA, X_trainB, y_trainA, y_trainB = train_test_split(X_train, y_train, test_size=0.5)
    
    # define and fit a RF.  use cross-validation to determine hyper-parameters.
    n_estimators = 30

    print '...training...'
    if True:
        depths = np.arange(1,50,10)
        ndepths = len(depths)
        sigma = np.zeros(ndepths)
        sigma68 = np.zeros(ndepths)
        for i,depth in enumerate(depths):
        #grf = GridSearchCV(rf, dict(max_depth=depths), cv=3).fit(X_train, y_train

            rf = rfRegressor(n_jobs=1, n_estimators=n_estimators, max_depth=depth)
            grf = rf.fit(X_train, y_train)
            predict_test = grf.predict(X_test)
            diff = predict_test-y_test
            this_sigma = np.std(diff)
            this_sigma68 = get_sigma68(diff)
            sigma[i] = this_sigma
            sigma68[i] = this_sigma68
            print i,depth,this_sigma, this_sigma68
  

        pl.figure(1)
        pl.clf()
        pl.plot(depths, sigma, 'b')
        pl.plot(depths, sigma68, 'g')
            
    else:
        rf = rfRegressor(n_jobs=3, n_estimators=n_estimators, min_samples_split=2)
        rf.fit(X_trainA, y_trainA)
        rf_predict_trainB = rf.predict(X_trainB)

        print '...training second thing...'
        rfB = rfRegressor(n_jobs=3, n_estimators=n_estimators, min_samples_split=2)
        metaX_trainB = np.hstack((X_trainB, rf_predict_trainB[:,np.newaxis]))
        rfB.fit(metaX_trainB, y_trainB)
       
        rf_predict_test = rf.predict(X_test)
        if True:
            metaX_test = np.hstack((X_test, rf_predict_test[:,np.newaxis]))
            predict_test = rfB.predict(metaX_test)
        else:
            #predict_test = rf_predict_test
            #predict_test = remove_bias(predict_test, y_test)
            fit = np.polyfit(rf_predict_train, y_train, 2)
            predict_test = np.polyval(fit, rf_predict_test)


        diff = predict_test-y_test


        this_sigma = np.std(diff)
        this_sigma68 = get_sigma68(diff)
        print this_sigma, this_sigma68



    #zlo = np.array([0.05, 0.15, 0.25, 0.35, 0.45])
    #zhi = np.array([0.15, 0.25, 0.35, 0.45, 0.55])
    zlo = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75])
    zhi = np.array([0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85])
    zbin = 0.5*(zlo+zhi)
    nbin = len(zbin)
    bias = np.zeros(nbin)
    bias2 = np.zeros(nbin)
    sigma = np.zeros(nbin)
    sigma68 = np.zeros(nbin)
    for i in np.arange(nbin):
        wh=np.where((y_test>=zlo[i]) & (y_test<zhi[i]))[0]
        nwh = len(wh)
        bias[i] = diff[wh].mean()
        bias2[i] = np.median(diff[wh])
        sigma[i] = diff[wh].std()
        sigma68[i] = get_sigma68(diff[wh])



    
    rlo = np.array([14,15,16,17,18,19,20])
    rhi = np.array([15,16,17,18,19,20,21])
    rbin = 0.5*(rlo+rhi)
    nrbin = len(rbin)
    rbias = np.zeros(nrbin)
    rbias2 = np.zeros(nrbin)
    rsigma = np.zeros(nrbin)
    rsigma68 = np.zeros(nrbin)
    r_ind = 0
    r_test = X_test[:,r_ind]
    for i in np.arange(nrbin):
        wh=np.where((r_test>=rlo[i]) & (r_test<rhi[i]))[0]
        nwh = len(wh)
        rbias[i] = diff[wh].mean()
        rbias2[i] = np.median(diff[wh])
        rsigma[i] = diff[wh].std()
        rsigma68[i] = get_sigma68(diff[wh])
    


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
    pl.plot(zbin, bias,'ro')
    pl.plot(zbin, bias2,'go')
    pl.legend(['mean','median'])
    pl.plot(zbin, bias,'r')
    pl.plot(zbin, bias2,'g')
    pl.plot([0,2],[0,0],'k--',linewidth=5)
    pl.xlim(y_train.min(), y_train.max())
    pl.ylim(np.array([-1,1])*0.05)
    pl.xlabel('specZ')
    pl.ylabel('bias')


    pl.subplot(325)
    pl.plot(zbin, sigma,'ro')
    pl.plot(zbin, sigma68,'go')
    pl.legend(['RMS','68%'])
    pl.plot(zbin, sigma,'r')
    pl.plot(zbin, sigma68,'g')
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


    pdb.set_trace()


def read_2mass(filename='2MASS_100000gals_0p05_to_1p0.txt', quick=True):
    import pickle
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

