from .ServerManager import update_source_data, clear_data
from .auxfuncs import *
from .widgets import HTMLDataFrame
from .LineNavigator import init_LineNavigator
import traitlets
from tumult.auxfuncs import tbdfits
from bokeh.models import ColumnDataSource as CDS
from .norm_form import norm_form
import pandas as pd

class StudioObs(traitlets.HasTraits):
    icon = 'satellite-dish'
    meta = None
    file = None
    normalised = False
    snr = 'unknown'

    def __init__(self, file, normalised=False, header=None, sep=r'\s+', **kwargs):
        super().__init__()
        self.normalised = normalised
        self.file = file

        if file.endswith('.fits'):
            out = tbdfits(file, **kwargs)
            if out is None:
                raise ValueError(f"Error reading file {file}")
            else:
                lam, flux, header = out
                source = CDS(dict(x=lam, y=flux))
                meta = pd.DataFrame(data=np.array(header.cards), columns=['key', 'value', 'comment'])
                self.meta = HTMLDataFrame(meta)

                # Check for SNR in the header
                snr_keys = [key for key in header if key.startswith("SNR")]
                if len(snr_keys) == 1:
                    self.snr = header[snr_keys[0]]

        else:
            df = pd.read_csv(file, header=header, sep=sep, **kwargs)
            source = CDS(dict(x=df.iloc[:,0].values, y=df.iloc[:,1].values))
        
        if not normalised:
            data = source.data.copy()
            xdata = source.data['x']
            ydata = source.data['y']
            mask = np.isfinite(ydata) * (ydata > 0)
            for key in data:
                data[key] = np.array(data[key])[mask]
            data['y'] /= data['y'].max()  # normalise y data
            source.data = data

        self.xdata = source.data['x'] # original data
        self.ydata = source.data['y'] # original data
        self.cntm = np.ones_like(self.xdata)  # normalisation spline
        self.source = source


    def chosen(self, studio):
        if not hasattr(self, 'navigator'):
            self.navigator = init_LineNavigator(self, studio)

        tabs = (self.navigator, self.info,)
        titles = ('Line Navigator', 'Info',)
        if not self.normalised:
            tabs += (norm_form.vbox,)
            titles += ('Norm Settings',)
        if hasattr(self, 'metadata'):
            tabs += (self.meta,)
            titles += ('Metadata',)

        instructions =  {'tabs':tabs,
                         'titles':titles,}
        return instructions

    
    @property
    def info(self):
        data = {'File': self.file,
                'Points': len(self.xdata),
                'xmin': self.xdata[0],
                'xmax': self.xdata[-1],
                'Normalised': self.normalised,
                'SNR': self.snr}
        df = pd.DataFrame(list(data.items()), columns=['key', 'value'])
        return HTMLDataFrame(df)
        