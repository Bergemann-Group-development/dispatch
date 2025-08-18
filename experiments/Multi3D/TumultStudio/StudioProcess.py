from .auxfuncs import *
from . import ServerManager

class StudioProcess(object):
    def __init__(self, label):
        self.label = label

    def delete(self):
        pass

    def chosen(self, studio):
        return None