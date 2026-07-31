#!/usr/bin/env python3
"""
probe_ccAFv2.py
===============
Make ccAFv2 2.0.6 importable under Keras 3, and print its API.

Root cause of the import error: ccAFv2 ships a Keras-2 model (legacy HDF5). Keras
3 rejects the saved loss `reduction='auto'` while rebuilding the model's TRAINING
config at import. We only do inference, so the loss/optimizer are irrelevant:
loading with compile=False skips that config. ccAFv2 calls
`keras.models.load_model(inPath)` at import without compile=False, so we patch
that function before importing ccAFv2.

If this prints "IMPORT OK", the same 4-line patch goes at the top of the
cell-cycle script and we're done. If it prints "IMPORT STILL FAILED", fall back
to Keras 2: `pip install tf-keras` then set TF_USE_LEGACY_KERAS=1 before import.

Run:  conda run -n NETWORK python probe_ccAFv2.py
"""
import os
import inspect
os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '2')

import keras
import keras.models
try:
    import keras.saving as _ksaving
except Exception:
    _ksaving = None

# ---- the fix: force compile=False on every model load ----
_orig_load = keras.models.load_model
def _load_no_compile(*args, **kwargs):
    kwargs['compile'] = False
    return _orig_load(*args, **kwargs)
keras.models.load_model = _load_no_compile
if _ksaving is not None and hasattr(_ksaving, 'load_model'):
    _ksaving.load_model = _load_no_compile

print(f"keras version: {keras.__version__}")

try:
    import ccAFv2
    print("IMPORT OK\n")
except Exception as e:
    print(f"IMPORT STILL FAILED: {type(e).__name__}: {str(e)[:300]}")
    print("\n-> fall back to Keras 2: pip install tf-keras, then set "
          "TF_USE_LEGACY_KERAS=1 before importing anything.")
    raise SystemExit(1)

# ---- reveal the API so the predict call can be wired correctly ----
print("ccAFv2 top-level:", [x for x in dir(ccAFv2) if not x.startswith('_')])
try:
    import ccAFv2.ccAFv2 as cc
    print("ccAFv2.ccAFv2 :", [x for x in dir(cc) if not x.startswith('_')])
except Exception as e:
    cc = ccAFv2
    print(f"(submodule introspection skipped: {e})")

print("\ncallables (signature: first doc line):")
seen = set()
for mod in (ccAFv2, cc):
    for name in dir(mod):
        if name.startswith('_') or name in seen:
            continue
        obj = getattr(mod, name)
        if callable(obj):
            seen.add(name)
            try:
                sig = str(inspect.signature(obj))
            except (ValueError, TypeError):
                sig = "(...)"
            doc = (getattr(obj, '__doc__', '') or '').strip().replace('\n', ' ')
            print(f"  {name}{sig}\n      {doc[:200]}")

print("\nDONE -- paste this whole output back.")
