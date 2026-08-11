#!/usr/bin/env python
# ============================================================================
# tobias_fork.py -- run TOBIAS on macOS
# ----------------------------------------------------------------------------
# WHY THIS EXISTS. TOBIAS's multiprocessing assumes the "fork" start method,
# which is the default on Linux. Since Python 3.8 macOS defaults to "spawn",
# which pickles the process object -- and TOBIAS passes a logger through, so
# every run dies immediately with:
#     _pickle.PicklingError: logger cannot be pickled
#   (tobias/tools/atacorrect.py -> logger.start_logger_queue())
# Forcing "fork" restores the Linux behaviour.
#
# OBJC_DISABLE_INITIALIZE_FORK_SAFETY is required alongside it: macOS aborts a
# forked child that touches an already-initialised Objective-C runtime, which
# numpy/matplotlib do. Without it the workers die with
# "objc[...]: +[__NSCFConstantString initialize] may have been in progress".
#
# MPLBACKEND=Agg for the same reason as run_peakqc.py -- the interactive
# "macosx" backend BLOCKS FOREVER when run detached, it does not error.
#
# Usage: identical to the TOBIAS CLI, e.g.
#   python pipeline/tobias_fork.py ATACorrect --bam ... --genome ... --cores 8
# ============================================================================
import multiprocessing
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES")

if __name__ == "__main__":
    multiprocessing.set_start_method("fork", force=True)
    from tobias.TOBIAS import main
    sys.exit(main())
