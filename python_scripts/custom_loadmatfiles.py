import scipy.io as spio


###############################
# Loading MAT-files
###############################
# Adapted from:
# https://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries/7418519
# So much better than loadmat by scipy

def load_matlab_data(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    mat_data = spio.loadmat(filename, mat_dtype=True, struct_as_record=False, squeeze_me=True)
    return _check_keys(mat_data)

def _check_keys(mat_data):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    list_keys_discard = ['__header__','__globals__','__version__']
    out_data={}
    for key in mat_data:
        if key not in list_keys_discard :
            if isinstance(mat_data[key], spio.matlab.mio5_params.mat_struct):
                out_data[key] = _todict(mat_data[key])
            else:
                out_data[key] = mat_data[key]
    return out_data        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    out_data = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            out_data[strg] = _todict(elem)
        else:
            out_data[strg] = elem
    return out_data


