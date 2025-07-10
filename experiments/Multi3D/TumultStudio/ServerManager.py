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
import atexit
from contextlib import closing
import numpy as np
import os

theme_file = module_dir + 'theme.json'

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
        
        #if self._server_thread and self._server_thread.is_alive():
        #    time.sleep(0.5)
        
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
        
        def start_server_thread():
            try:
                # Create new event loop for this thread
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                
                self._server = Server(
                    {'/': bokeh_app}, 
                    port=port, 
                    allow_websocket_origin=[f"localhost:{port}"]
                )
                self._server.start()
                self._server.io_loop.start()
            except Exception as e:
                print(f"Server error on port {port}: {e}")
                # Try to find another port if this one fails
                if auto_find_port:
                    try:
                        new_port = self.find_free_port(port + 1)
                        print(f"Retrying with port {new_port}")
                        self._current_port = new_port
                        # Recursive retry with new port
                        self.start_server(layout, new_port, auto_find_port=False)
                    except:
                        print("Failed to find alternative port")
        
        self._server_thread = threading.Thread(target=start_server_thread, daemon=True)
        self._server_thread.start()
        
        return port
    
    def get_current_port(self):
        """Get the port the server is currently running on"""
        return self._current_port

    def next_tick_update(self, update):
        self.doc.add_next_tick_callback(update)

def iframe_server(layout, PORT=5008, width='100%', height='100%'):
    server_manager.stop_existing_server()
    port = server_manager.start_server(layout, PORT)
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
    
def update_cds(cds, key, val):
    def update():
        cds.data[key] = val
        
    doc = server_manager.doc
    doc.add_next_tick_callback(update)

# Global server manager instance
server_manager = BokehServerManager()
    