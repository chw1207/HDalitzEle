import json
import numpy as np
from tools.plotUtils import color

# https://stackoverflow.com/a/70214170
class trial_count:
    """
    A decorator that will count and print how many times the decorated function was called
    """

    def __init__(self, inline_func):
        self.call_count = 0
        self.inline_func = inline_func

    def __call__(self, *args, **kwargs):
        self.call_count += 1
        self._print_call_count()
        return self.inline_func(*args, **kwargs)

    def _print_call_count(self):
        print(color.CYAN + "[Trial {}] ".format(self.call_count) + color.END, flush=True)

# https://stackoverflow.com/a/57915246
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)