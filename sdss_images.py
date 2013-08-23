import numpy as np
import matplotlib.pylab as pl
import ipdb
from os import system
from os.path import exists


def cutouts_to_batches():
    import pickle
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


