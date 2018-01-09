import pickle

def save(file,var):
    fd=open(file,'wb')
    pickle.dump(var,fd)
    fd.close()

def restore(file):
    fd=open(file,'rb')
    var=pickle.load(fd)
    fd.close()
    return var