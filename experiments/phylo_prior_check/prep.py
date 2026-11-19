"""Build genus-cluster presence matrices for one PanGBank genus, reusing pangbank_scale/worker.py (run_B's
export + MMseqs2 step). usage: prep.py GENUS ID,ID,... H5DIR OUTNPZ"""
import sys, os
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../pangbank_scale/worker.py')
genus, ids, h5dir, outnpz = sys.argv[1:5]
src = open(W).read().split("\nif a.mode == 'species':")[0]
# replace evaluation part of run_B: stop after G is built and save it
src = src.replace("    nsp = len(G)\n", "    nsp = len(G)\n    np.savez_compressed(OUTNPZ, Nn=Nn, names=np.array(names), **{f'G{i}': g for i, g in enumerate(G)})\n    shutil.rmtree(wd, ignore_errors=True)\n    return\n", 1)
sys.argv = ['worker.py', 'genus', genus, ids, h5dir, os.path.dirname(outnpz)]
g = {'__name__': 'w', 'OUTNPZ': outnpz}
exec(compile(src, W, 'exec'), g)
g['run_B'](genus, ids.split(','), h5dir, os.path.dirname(outnpz))
