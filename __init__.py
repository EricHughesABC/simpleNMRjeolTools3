# src/__init__.py
"""
SimpleNMR JEOL Tools
"""
__version__ = "1.0.0"

# Import everything from subdirectories into the main namespace
try:
    from core import *
except ImportError:
    pass

try:
    from gui import *
except ImportError:
    pass
    
try:
    from parsers import *
except ImportError:
    pass
    
try:
    from utils import *
except ImportError:
    pass

# Also make submodules available
try:
    from . import core
    from . import gui
    from . import parsers  
    from . import utils
except ImportError:
    pass