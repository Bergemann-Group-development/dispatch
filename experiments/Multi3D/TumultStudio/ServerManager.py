import os
import threading
import time
import asyncio
import socket
from .auxfuncs import *
from IPython.display import display, IFrame, HTML
import ipywidgets as widgets
from bokeh.server.server import Server
from bokeh.themes import Theme
from bokeh.models import *
from bokeh.io import curdoc
import json
import inspect
from contextlib import closing
import numpy as np
import os

theme_file = module_dir / 'theme.json'

# Global server manager instance+
def create_global():
    global server
    server = BokehServerManager()

def load_theme_from_file(filepath):
    with open(filepath) as f:
        return Theme(json=json.load(f))

class BokehServerManager:
    _instance = None
    _server = None
    _server_thread = None
    _current_port = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    def find_free_port(self, start_port=5008, max_attempts=50):
        """Find a free port starting from start_port"""
        for port in range(start_port, start_port + max_attempts):
            if self.is_port_free(port):
                return port
        raise RuntimeError(f"Could not find free port in range {start_port}-{start_port + max_attempts}")
    
    def is_port_free(self, port):
        """Check if a port is free"""
        try:
            with closing(socket.socket(socket.AF_INET, socket.SOCK_STREAM)) as sock:
                sock.settimeout(1)
                result = sock.connect_ex(('localhost', port))
                return result != 0  # Port is free if connection fails
        except:
            return False
    
    def stop_existing_server(self):
        """Properly stop existing server and thread"""
        if self._server:
            try:
                self._server.stop()
                if hasattr(self._server, 'io_loop'):
                    self._server.io_loop.stop()
            except:
                pass
            self._server = None
        
        self._server_thread = None
        self._current_port = None
    
    def start_server(self, layout, preferred_port=5008):
        """Start server on preferred port or find a free one"""
        # Stop any existing server first
        self.stop_existing_server()
        
        # Find a free port starting from preferred_port
        port = self.find_free_port(preferred_port)
        self._current_port = port
        
        def bokeh_app(doc):
            try:
                # Load theme if available
                if os.path.exists(theme_file):
                    doc.theme = load_theme_from_file(theme_file)
                else:
                    print(f"Warning: Theme file not found at {theme_file}")
                        
            except Exception as e:
                print(f"Error loading theme: {e}")
                import traceback
                traceback.print_exc()
            
            # Add your main layout
            self.doc = doc
            doc.add_root(layout)

        self._server = Server(
            {'/': bokeh_app},
            port=port,
            allow_websocket_origin=[f"localhost:{port}"],
            io_loop=None  # Use default event loop
        )
        
        # Start server
        self._server.start()
        self.loop = self._server.io_loop
        
        return port
    
    def get_current_port(self):
        """Get the port the server is currently running on"""
        return self._current_port

    def next_tick_update(self, update):
        sig = inspect.signature(update)
        nparams = len(sig.parameters.keys())
        params = [None,]*nparams

        if not hasattr(self, 'doc'):
            logging.info('server has no doc')
            update(*params)
        elif self.on_server():
            logging.info('immediate callback')
            update(*params)
        else:
            #logging.info(self.doc._pending_writes)
            @log_callback_errors
            def wrapper(*args, **kwargs):
                logging.info('tick callback')
                update(*params)
            self.doc.add_next_tick_callback(wrapper)

    def on_server(self):
        return self.doc._session_context()._session._pending_writes is not None

def iframe_server(layout, PORT=5008, width='100%', height='100%'):
    server.stop_existing_server()
    port = server.start_server(layout, PORT)
    # Wait for server to start
    #time.sleep(2)
    # Display iframe
    iframe = IFrame(src=f"http://localhost:{port}", width=width, height=height)
    # Create output widget to contain the iframe
    iframe_output = widgets.Output()
    # Display the iframe in the output widget
    with iframe_output:
        display(iframe)
    return iframe_output

def reset_xaxis(p, xx, margin=0):
    if len(xx) == 0: return
    xx = np.array(xx)
    xx = xx[~np.isnan(xx)]
    xmin = np.min(xx)
    xmax = np.max(xx)
    d = (xmax-xmin)*margin

    def update():
        p.x_range.start = xmin-d
        p.x_range.end   = xmax+d

    server.next_tick_update(update)
    logging.info(f'Resetting x-axis to {xmin-d} - {xmax+d}')
        
def reset_yaxis(p, yy, margin=0.05):
    if len(yy) == 0: return
    yy = np.array(yy)
    yy = yy[~np.isnan(yy)]
    ymin = np.min(yy)
    ymax = np.max(yy)
    d = (ymax-ymin)*margin

    def update():
        p.y_range.start = ymin-d
        p.y_range.end   = ymax+d

    server.next_tick_update(update)
    logging.info(f'Resetting y-axis to {ymin-d} - {ymax+d}')

def update_source_data(source):
    def decorator(update_func):
        @wraps(update_func)
        def wrapper():
            @log_callback_errors
            def update():
                current_data = source.data.copy()
                # Call the specific update function
                update_func(current_data)
                source.data = current_data
                        
            server.next_tick_update(update)
        return wrapper
    return decorator

def reorder_data(source, iold, inew):
    @update_source_data(source)
    def reorder(data):
        for key in data:
            item = data[key].pop(iold)
            data[key].insert(inew, item)
    reorder()

def clear_data(source):
    @update_source_data(source)
    def clear(data):
        for key in data:
            data[key] = []
    clear()


def spinning(spinner=None):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if spinner:
                def spin():
                    advance_spinner(spinner, 'value')
                spinner_value = spinner.value
                cid = server.doc.add_periodic_callback(spin, 100)

            result = func(*args, **kwargs)

            if spinner:
                server.doc.remove_periodic_callback(cid)
                spinner.value = spinner_value
            return result
        return wrapper
    return decorator
