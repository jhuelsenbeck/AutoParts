"""Coverage of the AutoParts simulation runs (e.g., from Analyses/:
    python3 coverage.py Simulation Output/simulation
compares the original runs in Simulation/ with the new ones in Output/simulation/), from the per-data-set "Results Summary" blocks
in the logs. Usage: coverage.py OLD_DIR NEW_DIR

For each data set: for every continuous parameter, whether its 95% CI contains the true
value (the YES/NO at the end of its line); for each of the four partitions, whether the
true partition is in the 95% credible set ("<- True Partition YES"; absent or NO = not).
Parameter coverage is averaged over all of a quantity's elements (as in the paper's
summarize_metrics.py), then over data sets. Conditions pool the _1 and _2 runs."""
import sys, re, glob, os, collections, math

QTY = [('ASRV', 'Alpha', 'Alpha'), ('Base frequencies', 'Pi', 'Pi'),
       ('Substitution rates', 'R', 'R'), ('Tree length', 'L', 'Length')]
PAPER = {  # Table 1: fixed 1, 3, 6 then gamma 1, 3, 6
    ('ASRV', 'Partition'): [0.96, 0.94, 1.00, 0.96, 0.94, 1.00],
    ('Base frequencies', 'Partition'): [1.00, 0.94, 0.94, 0.92, 0.98, 0.88],
    ('Substitution rates', 'Partition'): [0.96, 0.96, 0.98, 0.94, 0.94, 0.94],
    ('Tree length', 'Partition'): [0.88, 0.82, 0.96, 0.90, 0.88, 0.98],
    ('ASRV', 'Parameter'): [0.96, 0.96, 0.95, 0.94, 0.94, 0.97],
    ('Base frequencies', 'Parameter'): [0.97, 0.97, 0.97, 0.98, 0.98, 0.97],
    ('Substitution rates', 'Parameter'): [0.95, 0.97, 0.95, 0.97, 0.97, 0.97],
    ('Tree length', 'Parameter'): [0.92, 0.96, 0.94, 0.95, 0.97, 0.94],
}
CONDS = ['fixed_1.2', 'fixed_3.0', 'fixed_5.8', 'hyper_1.2', 'hyper_3.0', 'hyper_5.8']

cont = re.compile(r'^\s*\* (Alpha|Pi|R|L)\((\d+),(\w+)\): .*\(\S+ \S+ (YES|NO)\)\s*$')
rgfStart = re.compile(r'^\s*\* RGF\((\w+)\):')

def parse(path):
    """List of data sets, each {'part': {Length: bool,...}, 'param': {Alpha: [bools],...}}."""
    blocks = open(path, errors='replace').read().split('Results Summary:')[1:]
    out = []
    for b in blocks:
        part, param = {}, collections.defaultdict(list)
        cur = None
        for line in b.split('\n'):
            m = rgfStart.match(line)
            if m:
                cur = m.group(1)
                part[cur] = False
            m2 = cont.match(line)
            if m2:
                cur = None
                param[m2.group(1)].append(m2.group(4) == 'YES')
            if cur and 'True Partition YES' in line:
                part[cur] = True
        complete = len(part) == 4 and len(param['Alpha']) >= 6 and len(param['L']) >= 6
        if complete:
            out.append({'part': part, 'param': {k: v[:{'Alpha': 6, 'L': 6, 'Pi': 24, 'R': 36}[k]] for k, v in param.items()}})
    return out

def finalSummary(path):
    """The program's own end-of-run coverage, 'name -- mse ci', if the run finished."""
    res = {}
    for line in open(path, errors='replace'):
        m = re.match(r'^(RGF\(\w+\)|\w+\(\d+,\w+\)) -- (\S+) (\S+)\s*$', line)
        if m:
            res[m.group(1)] = float(m.group(3))
    return res

def summarize(d, pattern):
    table, ns, check = {}, {}, []
    for c in CONDS:
        sets = []
        for p in sorted(glob.glob(os.path.join(d, pattern.format(c)))):
            s = parse(p)
            sets += s
            fin = finalSummary(p)
            if fin and s:
                # check this parser against the program's own summary for each finished run
                for q, pre, rg in QTY:
                    mine = sum(x['part'][rg] for x in s) / len(s)
                    if abs(mine - fin.get(f'RGF({rg})', mine)) > 1e-9:
                        check.append((os.path.basename(p), rg, mine, fin[f'RGF({rg})']))
        ns[c] = len(sets)
        for q, pre, rg in QTY:
            if sets:
                table[(q, 'Partition', c)] = sum(x['part'][rg] for x in sets) / len(sets)
                table[(q, 'Parameter', c)] = sum(sum(x['param'][pre]) / len(x['param'][pre]) for x in sets) / len(sets)
    return table, ns, check

old, nOld, chkOld = summarize(sys.argv[1], '{}_[12].log')
new, nNew, chkNew = summarize(sys.argv[2], '{}_[12].log')
print('data sets per condition (old / new):', ', '.join(f'{c} {nOld[c]}/{nNew[c]}' for c in CONDS))
print('parser vs program summaries, mismatches: old', len(chkOld), ' new', len(chkNew), chkOld[:3], chkNew[:3])
print()
hdr = f"{'':34}" + ''.join(f'{c:>16}' for c in CONDS)
print(hdr); print(f"{'':34}" + ''.join(f"{'paper old  new':>16}" for c in CONDS))
for kind in ['Partition', 'Parameter']:
    for q, pre, rg in QTY:
        cells = []
        for i, c in enumerate(CONDS):
            o = old.get((q, kind, c)); n = new.get((q, kind, c))
            f = lambda v: '  - ' if v is None else f'{v:.2f}'
            cells.append(f'{PAPER[(q, kind)][i]:.2f} {f(o)} {f(n)}'.rjust(16))
        print(f'{q + " (" + kind + ")":34}' + ''.join(cells))
print()
# overall means and a binomial check for the partition coverages
for name, t, ns in [('old', old, nOld), ('new', new, nNew)]:
    for kind in ['Partition', 'Parameter']:
        vals = [t[(q, kind, c)] for q, _, _ in QTY for c in CONDS if (q, kind, c) in t]
        print(f'{name} {kind:9}: mean coverage {sum(vals)/len(vals):.3f} over {len(vals)} cells')
lo = [c for q, _, _ in QTY for c in CONDS if (q, 'Partition', c) in new and nNew[c] > 0
      and (new[(q, 'Partition', c)] - 0.95) / math.sqrt(0.95 * 0.05 / nNew[c]) < -2.5]
print('new partition cells more than 2.5 binomial s.e. below 0.95:',
      [(q, c, round(new[(q, 'Partition', c)], 2), nNew[c]) for q, _, _ in QTY for c in CONDS
       if (q, 'Partition', c) in new and nNew[c] > 0 and (new[(q, 'Partition', c)] - 0.95) / math.sqrt(0.95 * 0.05 / nNew[c]) < -2.5])
