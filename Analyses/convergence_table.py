"""Convergence of pairs of replicate runs, comparing two sets of runs of the same analyses.

    python3 convergence_table.py SHORT_DIR LONG_DIR [OUT.tsv] [--subset N]

e.g., from Analyses/:
    python3 convergence_table.py Output/empirical/skinks Output/empirical/skinks_long --subset 7
    python3 convergence_table.py Output/empirical/cynmix Output/empirical/cynmix_long

For every analysis whose two replicates (<name>_1 and <name>_2) are both finished in
LONG_DIR (a .done file), the same analysis in SHORT_DIR is summarized too. For each pair:
the largest split R-hat over the continuous parameters (all samples, and with the first
10% and 20% of each run discarded), how many parameters have R-hat above 1.01, the
smallest and median effective sample size with the first 20% discarded (the smaller of
the two replicates, per parameter), the mean log-likelihood of each replicate, and the
posterior mean number of components of each partition in each replicate. With
--subset N (for the skinks, 7: the ND1 third codon positions), also the fraction of
samples in which subset N has a tree-length table of its own, how many samples have it
merged with other subsets (in all, and in the first 10% of the run), and its mean tree
length (the third_* columns of the tsv file).
"""
import sys, os, glob, re
import numpy as np

args = sys.argv[1:]
THIRD = None    # the subset whose tree-length table is followed (--subset)
if '--subset' in args:
    k = args.index('--subset')
    THIRD = int(args[k + 1])
    del args[k:k + 2]
SHORT, LONG = args[0], args[1]
OUT = args[2] if len(args) > 2 else None
BURNIN = 0.2    # the fraction of each run discarded for the effective sample sizes

def load(path):
    with open(path) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        rows = [l.rstrip('\n').split('\t') for l in f if l.strip()]
    return hdr, [r for r in rows if len(r) == len(hdr)]

def ess(x):
    """Geyer's initial monotone sequence estimator"""
    x = np.asarray(x, float)
    n = len(x)
    x = x - x.mean()
    if not x.any():
        return float('nan')
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf /= acf[0]
    s = []
    for k in range(0, n - 1, 2):
        g = acf[k] + acf[k + 1]
        if g <= 0:
            break
        s.append(g)
    for i in range(1, len(s)):
        s[i] = min(s[i], s[i - 1])
    return n / max(-1.0 + 2.0 * sum(s), 1e-9)

def splitRhat(chains):
    halves = []
    for c in chains:
        n = len(c) // 2
        halves += [c[:n], c[n:2 * n]]
    n = min(len(h) for h in halves)
    x = np.array([h[:n] for h in halves])
    W = x.var(axis=1, ddof=1).mean()
    B = n * x.mean(axis=1).var(ddof=1)
    return float(np.sqrt(((n - 1) / n * W + B / n) / W)) if W > 0 else float('nan')

def numComponents(s):
    return len(set(s.split(',')))

def summarize(paths):
    (hdr, a), (_, b) = load(paths[0]), load(paths[1])
    out = {}
    cols = []
    for j, name in enumerate(hdr):
        if name in ('Gen', 'lnL') or name.startswith('RGF') or '(Alpha)' in name:
            continue
        try:
            x = np.array([float(r[j]) for r in a]); y = np.array([float(r[j]) for r in b])
        except ValueError:
            continue
        if x.std() == 0 and y.std() == 0:
            continue
        cols.append((name, x, y))
    for label, frac in [('all', 0.0), ('drop10', 0.1), ('drop20', 0.2)]:
        rh = []
        for name, x, y in cols:
            rh.append((splitRhat([x[int(len(x) * frac):], y[int(len(y) * frac):]]), name))
        rh.sort(reverse=True)
        out[f'rhat_{label}'] = rh[0][0]
        out[f'rhat_{label}_param'] = rh[0][1]
        out[f'n_above_1.01_{label}'] = sum(r > 1.01 for r, _ in rh)
    es = [min(ess(x[int(len(x) * BURNIN):]), ess(y[int(len(y) * BURNIN):])) for _, x, y in cols]
    out['n_params'] = len(cols)
    out['ess_min'] = float(np.nanmin(es))
    out['ess_median'] = float(np.nanmedian(es))
    out['lnL'] = [float(np.mean([float(r[1]) for r in rows])) for rows in (a, b)]
    out['K'] = {}
    for col in ['Length', 'Pi', 'R', 'Alpha']:
        i = hdr.index(f'RGF({col})')
        out['K'][col] = [float(np.mean([numComponents(r[i]) for r in rows])) for rows in (a, b)]
    if THIRD is not None and 'RGF(Length)' in hdr and f'L({THIRD},1)' in hdr:
        i = hdr.index('RGF(Length)')
        alone, merged = [], []
        for rows in (a, b):
            flags = [(lambda g: g.count(g[THIRD - 1]) == 1)(r[i].split(',')) for r in rows]
            alone.append(float(np.mean(flags)))
            # samples in which the third positions share a table: in all, and in the first 10%
            tenth = len(flags) // 10
            merged.append((len(flags) - sum(flags), tenth - sum(flags[:tenth])))
        out['third_alone'] = alone
        out['third_merged'] = merged
        k = hdr.index(f'L({THIRD},1)')
        out['third_length'] = [float(np.mean([float(r[k]) for r in rows])) for rows in (a, b)]
    return out

def analyses(d):
    names = set()
    for p in glob.glob(os.path.join(d, '*_1.done')):
        base = os.path.basename(p)[:-len('_1.done')]
        if os.path.exists(os.path.join(d, base + '_2.done')):
            names.add(base)
    return names

def setting(name):
    m = re.match(r'(\w+?)_(hyper_)?([\d.]+)$', name)
    return ('hyper' if m.group(2) else 'fixed', float(m.group(3))) if m else ('?', 0.0)

names = sorted(analyses(LONG), key=lambda n: (setting(n)[0], setting(n)[1]))
table = []
for name in names:
    for label, d in [('short', SHORT), ('long', LONG)]:
        paths = [os.path.join(d, f'{name}_{r}.tsv') for r in (1, 2)]
        if all(os.path.exists(p) for p in paths):
            table.append((name, label, summarize(paths)))

def f2(v): return '-' if v is None else f'{v:.2f}'
lines = ['| Analysis | Run | Largest R-hat (all) | > 1.01 | Largest R-hat (first 10% dropped) | > 1.01 | Largest R-hat (first 20% dropped) | > 1.01 | Smallest ESS (first 20% dropped) | Median ESS | Mean lnL (rep 1 / 2) | Subset alone | Samples with it merged: all (first 10%) | Subset tree length | Components: length, pi, r, alpha (rep 1 / 2) |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
for name, label, s in table:
    K = ', '.join(f"{s['K'][c][0]:.1f}/{s['K'][c][1]:.1f}" for c in ['Length', 'Pi', 'R', 'Alpha'])
    third = f"{s['third_alone'][0]:.2f} / {s['third_alone'][1]:.2f}" if 'third_alone' in s else '-'
    found = ' / '.join(f'{m} ({m10})' for m, m10 in s['third_merged']) if 'third_merged' in s else '-'
    tl = f"{s['third_length'][0]:.1f} / {s['third_length'][1]:.1f}" if 'third_length' in s else '-'
    lines.append(f"| {name} | {label} | {s['rhat_all']:.3f} ({s['rhat_all_param']}) | {s['n_above_1.01_all']} of {s['n_params']} | "
                 f"{s['rhat_drop10']:.3f} ({s['rhat_drop10_param']}) | {s['n_above_1.01_drop10']} | "
                 f"{s['rhat_drop20']:.3f} ({s['rhat_drop20_param']}) | {s['n_above_1.01_drop20']} | {s['ess_min']:.0f} | {s['ess_median']:.0f} | "
                 f"{s['lnL'][0]:.1f} / {s['lnL'][1]:.1f} | {third} | {found} | {tl} | {K} |")
print('\n'.join(lines))
if OUT:
    with open(OUT, 'w') as f:
        f.write('analysis\trun\trhat_all\trhat_all_param\tn_above_1.01_all\trhat_drop10\trhat_drop10_param\tn_above_1.01_drop10\t'
                'rhat_drop20\trhat_drop20_param\tn_above_1.01_drop20\tn_params\tess_min_drop20\tess_median\tlnL_1\tlnL_2\tthird_alone_1\tthird_alone_2\tthird_merged_1\tthird_merged_2\t'
                'third_length_1\tthird_length_2\t' + '\t'.join(f'K_{c}_{r}' for c in ['Length', 'Pi', 'R', 'Alpha'] for r in (1, 2)) + '\n')
        for name, label, s in table:
            g = lambda k, i: s[k][i] if k in s else ''
            f.write('\t'.join(str(v) for v in [name, label, s['rhat_all'], s['rhat_all_param'], s['n_above_1.01_all'],
                    s['rhat_drop10'], s['rhat_drop10_param'], s['n_above_1.01_drop10'],
                    s['rhat_drop20'], s['rhat_drop20_param'], s['n_above_1.01_drop20'], s['n_params'], s['ess_min'],
                    s['ess_median'], s['lnL'][0], s['lnL'][1], g('third_alone', 0), g('third_alone', 1),
                    g('third_merged', 0)[0] if 'third_merged' in s else '', g('third_merged', 1)[0] if 'third_merged' in s else '', g('third_length', 0), g('third_length', 1)]
                    + [s['K'][c][r] for c in ['Length', 'Pi', 'R', 'Alpha'] for r in (0, 1)]) + '\n')
