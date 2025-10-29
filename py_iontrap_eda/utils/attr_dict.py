class AttrDict(dict):
    '''
    Sets up dot.notation access to dictionary attributes
    dot.notation entails that dict.key = value
    '''
    
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__