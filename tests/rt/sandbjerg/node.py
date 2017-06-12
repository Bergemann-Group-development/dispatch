'''A class representing a compute node'''

class thread:
    def __init__(self,id):
        self.id=id
        self.color=0.1*(1+id)

class node:
    nthreads=6
    def __init__(self):
        self.threads=[]
        for i in range(self.nthreads):
            self.threads.append(thread(i))
