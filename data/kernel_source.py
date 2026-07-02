"""Locate/stage radiative kernel NetCDF files for the M2 lapse-rate module.

Reuses kernel data already downloaded by climkern (``python -m climkern
download``, run once inside the `pycfram-kern` conda env -- see
docs/plan.md §3.2) rather than re-implementing a Zenodo downloader. This
keeps `core/kernels.py` import-free of climkern/xesmf (docs/plan.md §3.3
"原则") while still giving every other machine on the project a simple way
to get at the same kernel files: run this once from within `pycfram-kern`,
then the plain `pycfram` env (or hqlx210) can read the staged copies with
zero extra dependencies.

Usage:
    python -m data.kernel_source stage CloudSat GFDL
    python -m data.kernel_source manifest
"""
import argparse
import hashlib
import os
import shutil
import sys

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KERNEL_DIR = os.path.join(PROJECT_ROOT, 'data', 'kernels')


def _md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def find_climkern_data_dir():
    """Locate the installed climkern package's data/kernels directory.

    Requires climkern to be importable (i.e. must be run from inside the
    `pycfram-kern` env, or any env where `pip install climkern` succeeded --
    climkern itself does not require xesmf to be *installed on disk*, only
    to *import the top-level package*, which does eagerly `import xesmf`.
    So in practice this only works inside `pycfram-kern`).
    """
    try:
        import importlib_resources
        base = importlib_resources.files('climkern')
    except Exception:
        try:
            import climkern
            base = os.path.dirname(climkern.__file__)
        except Exception as e:
            raise RuntimeError(
                "Cannot locate climkern package (import failed: %r). "
                "Run this from within the pycfram-kern conda env after "
                "`python -m climkern download`." % e)
    d = os.path.join(str(base), 'data', 'kernels')
    if not os.path.isdir(d):
        raise RuntimeError(
            "climkern found at %s but data/kernels/ missing -- "
            "run `python -m climkern download` first." % base)
    return d


def stage(names):
    """Copy TOA_<name>_Kerns.nc for each kernel name into data/kernels/,
    writing a manifest.txt with md5 sums alongside."""
    src_root = find_climkern_data_dir()
    os.makedirs(KERNEL_DIR, exist_ok=True)
    manifest_lines = []
    for name in names:
        src = os.path.join(src_root, name, 'TOA_%s_Kerns.nc' % name)
        if not os.path.exists(src):
            raise FileNotFoundError(
                "Kernel file not found: %s (available: %s)"
                % (src, os.listdir(src_root) if os.path.isdir(src_root) else '??'))
        dst_dir = os.path.join(KERNEL_DIR, name)
        os.makedirs(dst_dir, exist_ok=True)
        dst = os.path.join(dst_dir, 'TOA_%s_Kerns.nc' % name)
        shutil.copy2(src, dst)
        md5 = _md5(dst)
        manifest_lines.append("%s  %s  (from %s)" % (md5, dst, src))
        print("staged %s -> %s (md5=%s)" % (name, dst, md5))
    with open(os.path.join(KERNEL_DIR, 'manifest.txt'), 'a') as f:
        f.write('\n'.join(manifest_lines) + '\n')


def get_kernel_path(name):
    """Return the local staged path for a kernel, staging it on first use
    if climkern is importable and the file isn't already staged."""
    dst = os.path.join(KERNEL_DIR, name, 'TOA_%s_Kerns.nc' % name)
    if os.path.exists(dst):
        return dst
    stage([name])
    return dst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='cmd', required=True)
    p_stage = sub.add_parser('stage', help='Copy kernel files from climkern install')
    p_stage.add_argument('names', nargs='+', help='Kernel names, e.g. CloudSat GFDL')
    sub.add_parser('manifest', help='Print current manifest.txt')
    args = parser.parse_args()

    if args.cmd == 'stage':
        stage(args.names)
    elif args.cmd == 'manifest':
        mpath = os.path.join(KERNEL_DIR, 'manifest.txt')
        if os.path.exists(mpath):
            with open(mpath) as f:
                print(f.read())
        else:
            print("No manifest yet -- run `stage` first.")


if __name__ == '__main__':
    sys.exit(main())
