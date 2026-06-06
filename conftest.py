import os
import sys

# Make the copied src/ importable as `src.<module>` (matches the upstream tests).
sys.path.insert(0, os.path.dirname(__file__))
