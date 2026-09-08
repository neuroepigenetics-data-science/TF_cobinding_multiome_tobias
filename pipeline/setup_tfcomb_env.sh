#!/usr/bin/env bash
# ============================================================================
# setup_tfcomb_env.sh -- build the NATIVE arm64 environment for TF-COMB
# ----------------------------------------------------------------------------
# *** WHY A SEPARATE ENV, AND WHY NATIVE ***
# The `tobias` conda env -- and the whole miniconda3 install under it -- is
# x86_64 running under Rosetta on Apple Silicon: `platform.machine()` returns
# `x86_64` on arm64 hardware. That single fact is the entire reason TF-COMB
# would not install, and it is NOT the compiler problem the error message
# suggests. The chain is
#     tf-comb -> qnorm -> numba -> llvmlite
# and llvmlite 0.49.0 ships THIS wheel and no x86_64 macOS counterpart:
#     llvmlite-0.49.0-cp311-cp311-macosx_12_0_arm64.whl
# Under Rosetta pip finds no wheel, falls back to the sdist, tries to build
# against a matching LLVM, and dies naming llvmlite. Natively it downloads in
# seconds. llvmlite is the symptom; the architecture is the cause.
#
# The `tobias` env is NEVER repaired in place: it produced every result the
# project has, and TF-COMB does not need to share an env with it.
#
# *** WHY uv AND NOT conda ***
# `conda` itself is the Rosetta x86_64 install, so envs made with it inherit
# the wrong subdir unless forced. uv is already present as a native
# aarch64-apple-darwin binary and can fetch a managed arm64 CPython directly.
#
# NOTE the --python-preference only-managed flag: without it uv silently picks
# up the Intel Homebrew python at /usr/local/opt/python@3.11 and you get an
# x86_64 venv on arm64 hardware -- the exact trap this script exists to avoid.
# The arch assertion below is not decoration; it caught precisely that.
#
# *** THE tf-comb PATCH ***
# tf-comb 1.1.1 calls `pd.melt(..., var_name=["TF2"], ...)` in two places.
# pandas requires var_name to be a SCALAR and raises
#     ValueError: var_name=['TF2'] must be a scalar.
# on every pandas that satisfies tf-comb's own declared floor of >=2.3.0 (and
# on 2.2.x below it), so this cannot be pinned around -- downgrading far enough
# to accept a list means abandoning the declared dependency entirely. The bug
# is in tf-comb, the fix is two characters, and the result is identical to the
# obvious intent: a column named "TF2". Applied idempotently after install.
#
# pandas is held BELOW 3.0 as well: 3.0 is a major bump against a package last
# released for the 2.x line, and 2.3.3 is the newest version satisfying
# tf-comb's declared requirement.
#
# *** MULTIPROCESSING ***
# Scripts that USE this env must set the "fork" start method, the same reason
# pipeline/tobias_fork.py exists -- see the header of that file. TF-COMB's
# count_within()/market_basket() spawn workers, and under macOS "spawn" they
# re-import __main__; a script run from a heredoc or stdin then forks endlessly.
# Always run from a real .py file with an `if __name__ == "__main__":` guard.
#
# Usage:  bash pipeline/setup_tfcomb_env.sh [ENV_DIR]
#         default ENV_DIR is $HOME/venvs/tfcomb
# Verify: $ENV_DIR/bin/python -c "import platform; print(platform.machine())"
#         must print arm64.
# ============================================================================
set -euo pipefail

ENV_DIR="${1:-$HOME/venvs/tfcomb}"
LOCK="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tfcomb_env.lock.txt"

command -v uv >/dev/null || { echo "ERROR: uv not on PATH (expected ~/.local/bin/uv)"; exit 1; }

echo "==> creating native arm64 venv at $ENV_DIR"
uv venv "$ENV_DIR" --python-preference only-managed --python cpython-3.11-macos-aarch64

ARCH="$("$ENV_DIR/bin/python" -c 'import platform; print(platform.machine())')"
if [ "$ARCH" != "arm64" ]; then
    echo "ERROR: venv is $ARCH, not arm64 -- TF-COMB will not install. Aborting."
    exit 1
fi
echo "==> arch OK: $ARCH"

echo "==> installing"
if [ -f "$LOCK" ]; then
    echo "    from lockfile $LOCK"
    VIRTUAL_ENV="$ENV_DIR" uv pip install -r "$LOCK"
else
    echo "    resolving fresh (no lockfile found)"
    VIRTUAL_ENV="$ENV_DIR" uv pip install tf-comb 'pandas>=2.3,<3'
fi

echo "==> patching tf-comb melt calls (var_name must be scalar)"
OBJ="$ENV_DIR/lib/python3.11/site-packages/tfcomb/objects.py"
[ -f "$OBJ" ] || { echo "ERROR: $OBJ not found"; exit 1; }
if grep -q 'var_name=\["TF2"\]' "$OBJ"; then
    [ -f "$OBJ.orig" ] || cp "$OBJ" "$OBJ.orig"
    sed -i '' 's/var_name=\["TF2"\]/var_name="TF2"/g' "$OBJ"
    echo "    patched 2 call sites"
else
    echo "    already patched (or upstream fixed it -- check before assuming)"
fi

echo "==> verifying"
"$ENV_DIR/bin/python" - <<'PYEOF'
import platform
assert platform.machine() == "arm64", platform.machine()
import numpy, numba, llvmlite, pandas, pysam, qnorm, tobias, tfcomb
from numba import njit


@njit
def _f(x):
    return x * 2 + 1


assert _f(20.0) == 41.0            # numba must JIT, not merely import
import pandas as pd
pd.melt(pd.DataFrame({"TF1": ["a"], "x": [1]}), id_vars=["TF1"], var_name="TF2", value_name="n")
print("OK  arm64 | numpy", numpy.__version__, "| numba", numba.__version__,
      "| llvmlite", llvmlite.__version__, "| pandas", pandas.__version__,
      "| tf-comb", tfcomb.__version__)
PYEOF

echo "==> done. Activate with: source $ENV_DIR/bin/activate"
