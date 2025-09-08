import subprocess
import threading
from functools import wraps
import time
import sys
import logging
from .auxfuncs import *
from . import ServerManager
import signal
import ipywidgets as widgets
import cProfile
import pstats
import io
import pathlib
from pyinstrument import Profiler
import traceback

last_profiler = None
profile_checkbox = widgets.Checkbox(False, description='profile')

class CancellableExecution:
    def __init__(self, cancel_button=None, spinner=None, progress=None):
        self.cancel_event = threading.Event()
        self.process = None                  
        self._process_lock = threading.Lock()
        self.progress = progress
        self.reset_spinner = None
        
        if cancel_button:
            cancel_button.on_click(self.cancel_callback)
            
        if spinner:
            cid = ServerManager.server.doc.add_periodic_callback(spinner.spin, 100)
            
            def reset_spinner():
                try:
                    ServerManager.server.doc.remove_periodic_callback(cid)
                    spinner.reset()
                except Exception as e:
                    logging.warning(f"Failed to reset spinner: {e}")

            self.reset_spinner = reset_spinner

        if self.progress:
            self.progress_initial = self.progress.layout.display
            self.progress.layout.display = 'flex'

    @property
    def cancelled(self):
        return self.cancel_event.is_set()

    def run_command(self, command, **kwargs):
        if self.cancelled:
            return None
            
        try:
            with self._process_lock:
                self.process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    **kwargs
                )
            
            if hasattr(self, 'progress') and hasattr(self, 'monitor'):
                # Wait for completion or cancellation
                while self.process.poll() is None:
                    self.monitor()
                    threading.Event().wait(0.1)
            
            # Process completed naturally
            return self._log_and_return_process()
            
        except Exception as e:
            logging.error(f"Failed to run command: {e}")
            return None

    def _terminate_current_process(self):
        with self._process_lock:
            if self.process and self.process.poll() is None:
                self.process.terminate()
                try:
                    self.process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    self.process.kill()

    def _log_and_return_process(self):
        if not self.cancelled:
            stdout, stderr = self.process.communicate()
            
        if self.process.returncode == 0:
            final_message = "✅ Command completed successfully"
            if stdout:
                final_message += f"\nOutput: {stdout.strip()[:100]}..."
        elif self.process.returncode == -signal.SIGTERM:
            final_message = "🛑 Command was cancelled by user"
        elif self.process.returncode == -signal.SIGKILL:
            final_message = "🛑 Command was forcefully killed"
        else:
            final_message = f"❌ Command failed (exit code: {self.process.returncode})"
            if stderr:
                final_message += f"\nError: {stderr.strip()[:100]}..."
        
        logging.info(final_message)
        return self.process

    @log_callback_errors
    def cancel_callback(self, btn):
        logging.info("Cancellation requested")
        self.cancel_event.set()
        self._terminate_current_process()

    def finalize(self):
        if self.reset_spinner:
            self.reset_spinner()

        if self.progress:
            self.progress.layout.display = self.progress_initial
            self.progress.value = 0
            self.progress.max = 1


def execute_in_thread(cancel=None, spinner=None, progress=None):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Create context manager
            ctx = CancellableExecution(cancel, spinner, progress)

            def async_execution():
                global last_profiler
                profile = profile_checkbox.value
                if profile:
                    pr = Profiler()
                    pr.start()

                try:
                    func(*args, execution=ctx, **kwargs)
                    ctx.finalize()

                except Exception as e:
                    logging.error(f"Error in threaded execution: {e}")
                    traceback.print_exc()

                if profile:
                    pr.stop()
                    last_profiler = pr

            threading.Thread(target=async_execution, daemon=True).start()
        
        return wrapper
    return decorator


def get_profile_results(unicode=True, show_all=False, **kwargs):
    if last_profiler is None:
        return "No profiling data available."
    else:
        print(last_profiler.output_text(unicode=unicode, show_all=show_all, **kwargs))