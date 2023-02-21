
import os
import sys

DATA_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), './data'))

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import subsetFasta

