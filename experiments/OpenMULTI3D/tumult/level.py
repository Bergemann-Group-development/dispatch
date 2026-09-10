from .const import CC, EE, HH


class Level(object):
    def __init__(self, atom, i):
        self.i = i
        ecm, g, label, ion = atom._read(quotes="'")[:4]
        self.ecm   = float(ecm)
        self.g     = float(g)
        self.label = label
        self.ion   = int(ion)
        
        self.ev    = self.ecm*CC*HH/EE
        
        self.atom = atom
        self.trans = {}