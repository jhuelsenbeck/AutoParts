"""Compare the empirical analyses (gall wasps, skinks, hummingbirds) of the original runs
with the new ones, summarized the same way as the original partition_summaries.tsv files.

    python3 empirical_compare.py PARTITION_DISTANCE_PROGRAM [OUTPUT_TSV]

(run from Analyses/). For each data set, E(k), and treatment of the concentration
parameter (fixed or hyper), the samples of the two replicates are pooled (as in the
paper) and for each of the four partitions (tree length, base frequencies, substitution
rates, ASRV) we compute: the mean partition (with the partitiondistance program, as the
original summaries were), the posterior mean number of components, and the size of the
95% credible set of partition schemes. We also compute, for each pair of replicates,
the largest split R-hat over the continuous parameters, as a check of convergence.
"""
import sys, os, subprocess, tempfile, collections, concurrent.futures as cf
import numpy as np

PROGRAM = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else 'empirical_comparison.tsv'
DATA = {'cynmix': ['1.2', '2.0', '4.0', '6.0', '8.0', '10.0', '10.8'],
        'skinks': ['1.2', '2.0', '4.0', '6.0', '8.0', '10.0', '10.8'],
        'hummer': ['1.2', '2.0', '4.0', '6.0', '8.0', '8.8']}
PARMS = [('length', 'Length'), ('pi', 'Pi'), ('r', 'R'), ('alpha', 'Alpha')]

def files(ds, v, hyper, new):
    tag = f'{ds}_hyper_{v}' if hyper else f'{ds}_{v}'
    if new:
        return [f'Output/empirical/{ds}/{tag}_{r}.tsv' for r in (1, 2)]
    return [f'Empirical/{ds}_analyses/{tag}_{r}.out.tsv' for r in (1, 2)]

def load(path):
    with open(path) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        rows = [l.rstrip('\n').split('\t') for l in f if l.strip()]
    return hdr, rows

def canon(rgf):
    """restricted growth function of a partition given as a list of labels"""
    m, out = {}, []
    for x in rgf:
        m.setdefault(x, len(m) + 1)
        out.append(m[x])
    return tuple(out)

def splitRhat(chains):
    """largest-variance-safe split R-hat (Gelman et al.) for a list of 1-d arrays"""
    halves = []
    for c in chains:
        n = len(c) // 2
        halves += [c[:n], c[n:2 * n]]
    m, n = len(halves), min(len(h) for h in halves)
    x = np.array([h[:n] for h in halves])
    W = x.var(axis=1, ddof=1).mean()
    B = n * x.mean(axis=1).var(ddof=1)
    if W == 0:
        return float('nan')
    return float(np.sqrt(((n - 1) / n * W + B / n) / W))

def meanPartition(name, rgfs):
    with tempfile.NamedTemporaryFile('w', suffix='.part', delete=False) as f:
        f.write(name + '\n' + ''.join(f'[{x}]\n' for x in rgfs))
        path = f.name
    out = subprocess.run([PROGRAM, path], capture_output=True, text=True).stdout
    os.unlink(path)
    for line in out.split('\n'):
        if line.startswith('Average partition = '):
            return line.split('(')[1].split(')')[0]
    return '?'

def summarize(ds, v, hyper, new):
    paths = files(ds, v, hyper, new)
    if not all(os.path.exists(p) for p in paths):
        return None
    loaded = [load(p) for p in paths]
    hdr = loaded[0][0]
    res = {}
    for short, col in PARMS:
        i = hdr.index(f'RGF({col})')
        rgfs = [r[i] for _, rows in loaded for r in rows]
        parts = collections.Counter(canon(x.split(',')) for x in rgfs)
        n = sum(parts.values())
        ks = [max(p) for p in parts.elements()]
        cum, size = 0.0, 0
        for p, c in parts.most_common():
            cum += c / n
            size += 1
            if cum >= 0.95:
                break
        res[short] = {'rgfs': rgfs, 'meanK': float(np.mean(ks)), 'ci': size}
    # convergence of the continuous parameters between the two replicates
    worst = (0.0, '')
    for j, name in enumerate(hdr):
        if name in ('Gen', 'lnL') or name.startswith('RGF') or '(Alpha)' in name:
            continue
        try:
            chains = [np.array([float(r[j]) for r in rows]) for _, rows in loaded]
        except (ValueError, IndexError):
            continue        # not a numeric column (e.g. an empty trailing column)
        rh = splitRhat(chains)
        if rh == rh and rh > worst[0]:
            worst = (rh, name)
    return res, worst

def main():
    jobs = [(ds, v, h, new) for ds in DATA for v in DATA[ds] for h in (False, True) for new in (False, True)]
    summaries = {}
    for key in jobs:
        summaries[key] = summarize(*key)
    # mean partitions with the partitiondistance program, in parallel
    tasks = {}
    with cf.ThreadPoolExecutor(max_workers=os.cpu_count()) as ex:
        for key, s in summaries.items():
            if s is None:
                continue
            for short, _ in PARMS:
                tasks[(key, short)] = ex.submit(meanPartition, short, s[0][short]['rgfs'])
        means = {k: f.result() for k, f in tasks.items()}
    with open(OUT, 'w') as f:
        f.write('dataset\tEk\tconcentration\tparameter\told_mean_partition\tnew_mean_partition\tsame\t'
                'old_meanK\tnew_meanK\told_CI_size\tnew_CI_size\n')
        for ds in DATA:
            for h in (False, True):
                for v in DATA[ds]:
                    o, nw = summaries[(ds, v, h, False)], summaries[(ds, v, h, True)]
                    for short, _ in PARMS:
                        om = means.get(((ds, v, h, False), short), '-')
                        nm = means.get(((ds, v, h, True), short), '-')
                        same = om != '-' and nm != '-' and canon(om.split(',')) == canon(nm.split(','))
                        f.write(f"{ds}\t{v}\t{'hyper' if h else 'fixed'}\t{short}\t{om}\t{nm}\t{same}\t"
                                f"{o[0][short]['meanK'] if o else float('nan'):.2f}\t{nw[0][short]['meanK'] if nw else float('nan'):.2f}\t"
                                f"{o[0][short]['ci'] if o else '-'}\t{nw[0][short]['ci'] if nw else '-'}\n")
    # convergence
    print('largest split R-hat over the continuous parameters (two replicates), old / new:')
    for ds in DATA:
        for h in (False, True):
            cells = []
            for v in DATA[ds]:
                o, nw = summaries[(ds, v, h, False)], summaries[(ds, v, h, True)]
                cells.append(f"{v}: {o[1][0] if o else float('nan'):.3f}/{nw[1][0] if nw else float('nan'):.3f}")
            print(f"  {ds:7} {'hyper' if h else 'fixed':5}  " + '  '.join(cells))
    print('wrote', OUT)

if __name__ == '__main__':
    main()
