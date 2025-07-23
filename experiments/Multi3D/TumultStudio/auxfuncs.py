import bokeh.palettes
from functools import wraps
import logging
import threading
import numpy as np
import os
from line_profiler import LineProfiler

module_dir = os.path.dirname(os.path.abspath(__file__)) + os.sep
logger = logging.getLogger(__name__)

def run_in_thread(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        out = {}
        def runner():
            out['result'] = func(*args, **kwargs)

        thread = threading.Thread(target=runner)
        thread.start()
        thread.out = out
        return thread  # so you can optionally join or check .is_alive()
    return wrapper

def log_callback_errors(func):
    """Minimal decorator to log exceptions in callbacks."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.exception(f"Error in callback '{func.__name__}'")
    return wrapper

    from line_profiler import LineProfiler

def profile_lines(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        profiler = LineProfiler()
        profiler.add_function(func)
        profiler.enable_by_count()
        
        try:
            result = func(*args, **kwargs)
        finally:
            profiler.disable_by_count()
            profiler.print_stats()
        
        return result
    
    return wrapper

class ColorCycle():
    def __init__(self, palette='Light'):
        pal = bokeh.palettes.__dict__[palette]
        n = max(list(pal.keys()))
        self.palette = pal[n]
        self.n = n

    def __getitem__(self, index):
        i = np.mod(index, self.n)
        return self.palette[i]

def get_crop_range(data, range_type, nmax=None):
    vmin = range_type.start
    vmax = range_type.end
    imin = np.searchsorted(data, vmin)
    imax = np.searchsorted(data, vmax)-1
    if nmax:
        step = max(1, (imax-imin)//nmax)
        return slice(imin, imax, step)
    else:
        return slice(imin, imax)

def find_smallest_missing(nums):
    return min(set(range(0, max(nums, default=0) + 2)) - set(nums))

def get_common_suffix_files(files, suffixes):
    return np.unique([[f[:-len(s)] for f in files if f.endswith(s)] for s in suffixes])

def replace_nth(string, sub, replacement, n):
    parts = string.split(sub)
    if n < 1 or n >= len(parts):
        return string  # nth occurrence doesn't exist
    return sub.join(parts[:n]) + replacement + sub.join(parts[n:])
    
def combine_selections(selected_indices, current_indices, mode='replace'):
    """
    Combine two lists of indices based on different selection modes.
    
    Parameters:
    -----------
    selected_indices : list
        New indices to be combined
    current_indices : list  
        Existing selected indices (e.g., source.selected.indices)
    mode : str
        Selection mode - 'replace', 'append', 'subtract', 'intersect', or 'toggle'
        
    Returns:
    --------
    list : Combined indices sorted in ascending order
    """
    # Convert to sets for efficient operations
    new_set = set(selected_indices)
    current_set = set(current_indices)
    
    if mode == 'replace':
        # Replace current selection with new selection
        result = new_set
        
    elif mode == 'append' or mode == 'union':
        # Add new indices to current selection (union)
        result = current_set.union(new_set)
        
    elif mode == 'subtract' or mode == 'remove':
        # Remove new indices from current selection
        result = current_set.difference(new_set)
        
    elif mode == 'intersect':
        # Keep only indices that are in both selections
        result = current_set.intersection(new_set)
        
    elif mode == 'toggle' or mode == 'xor':
        # Toggle: add if not present, remove if present (symmetric difference)
        result = current_set.symmetric_difference(new_set)
        
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'replace', 'append', 'subtract', 'intersect', or 'toggle'")
    
    # Convert back to sorted list
    return sorted(list(result))
