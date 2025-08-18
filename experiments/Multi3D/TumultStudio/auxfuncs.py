import bokeh.palettes
from functools import wraps
import logging
import threading
import numpy as np
import pandas as pd
import ipywidgets as widgets
from line_profiler import LineProfiler
import tempfile
import shutil
import atexit
import pathlib
import sys
from io import StringIO
import traceback

module_dir = pathlib.Path(__file__).parent
logger = logging.getLogger(__name__)
spinner_chars = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']

def run_in_thread(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        out = {
            'result': None,
            'error': None,
            'traceback': None,
            'stdout': None,
            'stderr': None,
            'logs': []
        }
        
        def runner():
            # Capture stdout/stderr
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            stdout_capture = StringIO()
            stderr_capture = StringIO()
            level = logging.getLogger().getEffectiveLevel()
            
            # Capture logging
            log_capture = StringIO()
            log_handler = logging.StreamHandler(log_capture)
            log_handler.setLevel(level)
            
            # Get the root logger
            root_logger = logging.getLogger()
            original_level = root_logger.level
            root_logger.addHandler(log_handler)
            root_logger.setLevel(level)
            
            try:
                sys.stdout = stdout_capture
                sys.stderr = stderr_capture
                
                out['result'] = func(*args, **kwargs)
                
            except Exception as e:
                out['error'] = e
                out['traceback'] = traceback.format_exc()
                logging.error(f"Error in thread: {e}")
                
            finally:
                # Restore everything
                sys.stdout = old_stdout
                sys.stderr = old_stderr
                root_logger.removeHandler(log_handler)
                root_logger.setLevel(original_level)
                
                # Store captured output
                out['stdout'] = stdout_capture.getvalue()
                out['stderr'] = stderr_capture.getvalue()
                out['logs'] = log_capture.getvalue()
        
        thread = threading.Thread(target=runner)
        thread.start()
        thread.out = out
        
        # Add the same convenience methods as before
        def get_result(timeout=None):
            thread.join(timeout)
            if thread.is_alive():
                raise TimeoutError(f"Thread did not complete within {timeout} seconds")
            if out['error']:
                raise out['error']
            return out['result']
        
        def print_output():
            if out['stdout']:
                print("=== Thread stdout ===")
                print(out['stdout'])
            if out['stderr']:
                print("=== Thread stderr ===", file=sys.stderr)
                print(out['stderr'], file=sys.stderr)
            if out['logs']:
                print("=== Thread logs ===")
                print(out['logs'])
            if out['error']:
                print("=== Thread error ===", file=sys.stderr)
                print(out['traceback'], file=sys.stderr)
        
        thread.get_result = get_result
        thread.print_output = print_output

        return thread
    
    return wrapper

def log_callback_errors(func):
    """Minimal decorator to log exceptions in callbacks."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        if len(args)>0 and isinstance(args[0], widgets.Button):
            logging.info(f"Button '{args[0].description}' clicked.")

        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.exception(f"Error in callback '{func.__name__}'")
    return wrapper


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


def debounce(delay=0.5):
    def decorator(func):
        timer = None
        
        @wraps(func)
        def wrapper(*args, **kwargs):
            nonlocal timer
            
            # Cancel existing timer
            if timer is not None:
                timer.cancel()
                
            # Start new timer
            timer = threading.Timer(delay, func, args=args, kwargs=kwargs)
            timer.start()
            
        return wrapper
    return decorator

class ColorCycle():
    def __init__(self, palette='Light'):
        pal = bokeh.palettes.__dict__[palette]
        n = max(list(pal.keys()))
        self.palette = np.array(pal[n])
        self.n = n

    def __getitem__(self, index):
        i = np.mod(index, self.n)
        return self.palette[i]
    
    def get_new_color(self, used_colors):
        """Get the next unused color from the palette."""
        colors = self[np.arange(self.n)]
        used_colors = np.array(used_colors)
        i = np.argmin([np.sum(used_colors == c) for c in colors])
        return self[i]
    
def slicesorted(data, vmin, vmax):
    imin = np.searchsorted(data, vmin, side='left')
    imax = np.searchsorted(data, vmax, side='right')
    return slice(imin, imax)

def get_crop_range(data, range_type, nmax=None, margin=0):
    vmin = range_type.start
    vmax = range_type.end
    dv = (vmax - vmin)*margin
    slc = slicesorted(data, vmin-dv, vmax+dv)
    
    if nmax:
        step = max(1, (slc.stop-slc.start)//nmax)
        slc = slice(slc.start, slc.stop, step)
    
    return slc

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def beyond_sorted(data, val):
    return (val < data[0]) or (val > data[-1])

def find_smallest_missing(nums):
    return min(set(range(0, max(nums, default=0) + 2)) - set(nums))

def get_common_suffix_files(files, suffixes):
    return np.unique([[f.with_suffix('') for f in files if f.suffix == s] for s in suffixes])

def replace_nth(string, sub, replacement, n):
    parts = string.split(sub)
    if n < 1 or n >= len(parts):
        return string  # nth occurrence doesn't exist
    return sub.join(parts[:n]) + replacement + sub.join(parts[n:])
    
def unique_folder_name(folder):
    folder = pathlib.Path(folder)
    
    if folder.exists():
        i = 1
        while True:
            new_folder = folder.with_name(f"{folder.name}_{i}")
            if not new_folder.exists():
                break
            i += 1
        folder = new_folder
    
    return folder

def get_unique_name(name, existing_names, fmt):
    """Generate a unique name by appending a number if necessary."""
    if name not in existing_names:
        return name
    
    i = 1
    while True:
        new_name = name + fmt.format(i)
        if new_name not in existing_names:
            return new_name
        i += 1


def create_temp_dir_for_kernel(base_path):
    """Create a temporary directory that gets cleaned up when the kernel exits"""
    temp_dir = pathlib.Path(tempfile.mkdtemp(dir=base_path))
    data_path = 'data' / temp_dir
    data_path.mkdir(parents=True, exist_ok=True)
    logging.info(f"Created temporary directory: {temp_dir}")
    logging.info(f"Created temporary directory: {data_path}")
    
    def cleanup():
        for td in [temp_dir, data_path]:
            if td.exists():
                logging.info(f"Cleaning up temporary directory: {td}")
                shutil.rmtree(td)
    
    # Register cleanup for kernel shutdown
    atexit.register(cleanup)
    return temp_dir


def link_files(source_dir, target_dir, pattern):
    """Create symbolic links for all files starting with 'pop_' from source to target"""
    source_path = pathlib.Path(source_dir)
    target_path = pathlib.Path(target_dir)
    
    # Create target directory if it doesn't exist
    target_path.mkdir(parents=True, exist_ok=True)
    
    # Find all files with pattern
    files = source_path.glob(pattern)
    
    for source_file in files:
        if source_file.is_file():  # Only link files, not directories
            target_link = target_path / source_file.name
            
            # Remove existing link/file if it exists
            if target_link.exists() or target_link.is_symlink():
                target_link.unlink()
            
            # Create symbolic link
            target_link.symlink_to(source_file.resolve())
            logging.info(f"Linked: {source_file.name} -> {target_link}")


def advance_spinner(spinner, key):
    val = getattr(spinner, key, ' ')
    s = val[0]
    if s not in spinner_chars:
        val = spinner_chars[0] + ' ' + val
    else:
        i = spinner_chars.index(s)
        i = (i + 1) % len(spinner_chars)
        s = spinner_chars[i]
        val = s + val[1:]
    setattr(spinner, key, val)


def cross_match_dataframes(df1, df2, on, by=None, direction='nearest', tolerance=0, **kwargs):
    """
    Merge columns from df2 to df1 using exact matches on 'by' columns and 
    nearest/backward/forward matches on 'on' column.
    
    Parameters:
    -----------
    df1 : pd.DataFrame
        Left dataframe (base dataframe to merge to)
    df2 : pd.DataFrame
        Right dataframe (source of additional columns)
    on : str
        Column name for nearest matching (usually numerical/datetime)
    by : str or list of str, optional
        Column(s) that must match exactly before applying nearest match
    direction : str, default 'nearest'
        Direction for merge_asof ('backward', 'forward', 'nearest')
    tolerance : pd.Timedelta or float, default 0
        Maximum distance between matches for the 'on' column (passed to merge_asof)
    **kwargs
        Additional arguments passed to pd.merge_asof
    
    Returns:
    --------
    pd.DataFrame
        Merged dataframe
    """
    
    # Convert by to list if it's a single column
    if by is not None:
        if isinstance(by, str):
            by = [by]
    
    if by:
        dfs = []
        # Group both dataframes by the exact match columns
        for vals, g1 in df1.groupby(by):
            # Handle single column case (vals will be scalar, not tuple)
            if len(by) == 1:
                filter_condition = df2[by[0]] == vals
            else:
                # Multiple columns - create filter condition
                filter_condition = pd.Series([True] * len(df2))
                for col, val in zip(by, vals):
                    filter_condition &= (df2[col] == val)
            
            g2 = df2[filter_condition].drop(columns=by)
            
            # Skip if no matching rows in df2
            if g2.empty:
                # Add rows with NaN values for df2 columns
                g1_extended = g1.copy()
                for col in df2.columns:
                    if col not in g1.columns:
                        g1_extended[col] = pd.NA
                dfs.append(g1_extended)
                continue
            
            # Sort both groups by the 'on' column
            g1 = g1.sort_values(on)
            g2 = g2.sort_values(on)
            
            # Perform merge_asof
            matched = pd.merge_asof(g1, g2, on=on, direction=direction, tolerance=tolerance, **kwargs)
            dfs.append(matched)
        
        df = pd.concat(dfs, ignore_index=True)
    else:
        # No groupby needed - direct merge_asof
        df1_sorted = df1.sort_values(on)
        df2_sorted = df2.sort_values(on)
        
        df = pd.merge_asof(df1_sorted, df2_sorted, on=on, direction=direction, tolerance=tolerance, **kwargs)
    
    return df