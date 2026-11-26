import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from hypothesis import settings

# Derandomized so a run on a loaded machine or in CI tests the same cases as
# on a laptop; no deadline because step time depends on machine load.
settings.register_profile("repo", derandomize=True, deadline=None, max_examples=60)
settings.load_profile("repo")
