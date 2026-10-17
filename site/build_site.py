"""Build the static GitHub Pages site (docs/index.html) from the public summary results.

Reads ONLY files that are already in this repository (results/, figures/). No raw or
patient-level clinical data are read, no model is fitted and no new statistic is computed,
except simple summaries of already published numbers (bootstrap percentiles, histogram bins).

Every headline number is checked against the values recorded in the step 5-8 progress
notes; the build stops if any check fails.

Usage (from the repository root):
    python site/build_site.py
Optional input: results/09c/km_curves_09c.json (Kaplan-Meier step coordinates, aggregate).
If it is absent, the page shows the published KM figures instead of interactive curves.
"""
import base64
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / 'results'
F = ROOT / 'figures'
TEMPLATE = ROOT / 'site' / 'template.html'
OUT = ROOT / 'docs' / 'index.html'
REPO_URL = 'https://github.com/gurbuzmericarda0-arch/dlbcl-prognostic-score'
SITE_VERSION = '09c-1.0'

COHORT_LABEL = {'GSE181063_test': 'GSE181063 test', 'GSE10846': 'GSE10846', 'GSE87371': 'GSE87371'}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def clean(x):
    """JSON-safe: NaN/inf -> None, numpy -> python, rounded floats."""
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return None if (math.isnan(x) or math.isinf(x)) else float(f'{x:.6g}')
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def num(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) or math.isinf(v) else v


checks = []


def check(label, value, expected, tol):
    ok = value is not None and abs(value - expected) <= tol
    checks.append({'check': label, 'value': value, 'expected': expected, 'ok': ok})


# ------------------------------------------------------------------ provenance
lock = json.loads((R / '04' / 'plan_lock_04.json').read_text()) if (R / '04' / 'plan_lock_04.json').exists() else {}
try:
    commit = subprocess.run(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], capture_output=True,
                            text=True, check=True).stdout.strip()
except Exception:
    commit = None
meta = {
    'site_version': SITE_VERSION,
    'repo_url': REPO_URL,
    'source_commit': commit,
    'plan_lock_sha256': sha256(R / '04' / 'plan_lock_04.json'),
    'model_sha256': sha256(R / '04' / 'model_04.json'),
}
check('plan lock prefix', 1.0 if meta['plan_lock_sha256'].startswith('8d87cf887e2d') else 0.0, 1.0, 0)
check('model prefix', 1.0 if meta['model_sha256'].startswith('661e2b80be77') else 0.0, 1.0, 0)

# ------------------------------------------------------------------ cohorts (Table S1)
s1 = pd.read_csv(R / '09b' / 'table_s1_characteristics_wide_09b.csv')
cohorts = {'columns': list(s1.columns), 'rows': s1.fillna('').astype(str).values.tolist()}

# ------------------------------------------------------------------ nested CV (training)
ncv = json.loads((R / '04' / 'nested_cv_summary_04b.json').read_text())
nested = {'mean_c': ncv['min']['mean_c'], 'sd_c': ncv['min']['sd_c'], 'min_c': ncv['min']['min_c'],
          'max_c': ncv['min']['max_c'], 'units': ncv['units']}
check('nested CV C', nested['mean_c'], 0.6117, 0.0005)

# ------------------------------------------------------------------ validation (steps 5-6)
s5 = pd.read_csv(R / '05' / 'sensitivity_05.csv')
s6 = pd.read_csv(R / '06' / 'sensitivity_06.csv')
validation = []
for _, r in s5.iterrows():
    primary = r['analysis'].startswith('primary')
    validation.append({
        'cohort': 'GSE181063 test', 'role': 'internal', 'primary': primary,
        'analysis': 'primary' if primary else r['analysis'],
        'n': int(r['n']), 'deaths': int(r['deaths']), 'hr': r['hr_per_sd'], 'lo': r['hr_ci_low'],
        'hi': r['hr_ci_high'], 'p': r['p_wald'], 'c': r['harrell_c'], 'c_lo': num(r['c_ci_low']),
        'c_hi': num(r['c_ci_high']), 'uno_c': r['uno_c'], 'tau': 60.0, 'logrank_p': r['logrank_p'],
        'calibration_slope': num(r['calibration_slope'])})
for _, r in s6.iterrows():
    primary = r['analysis'] == 'primary'
    validation.append({
        'cohort': r['cohort'], 'role': 'external' if r['cohort'] == 'GSE10846' else 'supporting external',
        'primary': primary, 'analysis': r['analysis'], 'n': int(r['n']), 'deaths': int(r['deaths']),
        'hr': r['hr_per_sd'], 'lo': r['hr_ci_low'], 'hi': r['hr_ci_high'], 'p': r['p_wald'],
        'c': r['harrell_c'], 'c_lo': num(r['c_ci_low']), 'c_hi': num(r['c_ci_high']),
        'uno_c': r['uno_c'], 'tau': num(r['uno_tau_months']), 'logrank_p': r['logrank_p'],
        'hr_high_vs_low': num(r['hr_high_vs_low'])})
iv = json.loads((R / '05' / 'internal_validation_05.json').read_text())['primary']
validation[0]['hr_high_vs_low'] = iv['hr_high_vs_low']['hr']
prim = {v['cohort']: v for v in validation if v['primary']}
check('H1 HR test', prim['GSE181063 test']['hr'], 1.93, 0.005)
check('H1 HR CI low', prim['GSE181063 test']['lo'], 1.49, 0.005)
check('test C', prim['GSE181063 test']['c'], 0.664, 0.0005)
check('H2 HR GSE10846', prim['GSE10846']['hr'], 1.69, 0.005)
check('C GSE10846', prim['GSE10846']['c'], 0.659, 0.0005)
check('HR GSE87371', prim['GSE87371']['hr'], 1.49, 0.005)
check('C GSE87371', prim['GSE87371']['c'], 0.630, 0.0005)
check('n GSE10846', prim['GSE10846']['n'], 229, 0)
check('mini-CHOP HR', [v for v in validation if 'mini' in v['analysis']][0]['hr'], 1.39, 0.005)

# ------------------------------------------------------------------ KM group summaries
km_groups = []
k5 = pd.read_csv(R / '05' / 'km_test_05.csv')
for _, r in k5.iterrows():
    km_groups.append({'cohort': 'GSE181063 test',
                      'grouping': 'median' if r['grouping'].startswith('training median') else 'tertile',
                      'group': r['group'], 'n': int(r['n']), 'deaths': int(r['deaths']),
                      'median_os': num(r['median_os_months']), 'os5': num(r['os_5y']),
                      'os5_lo': num(r['os_5y_ci_low']), 'os5_hi': num(r['os_5y_ci_high'])})
k6 = pd.read_csv(R / '06' / 'km_external_06.csv')
for _, r in k6.iterrows():
    ci = json.loads(str(r['os_5y_ci_95']).replace('None', 'null')) if isinstance(r['os_5y_ci_95'], str) else [None, None]
    km_groups.append({'cohort': r['cohort'], 'grouping': 'median' if 'median' in r['grouping'] else 'tertile',
                      'group': r['group'], 'n': int(r['n']), 'deaths': int(r['deaths']),
                      'median_os': num(r['median_os_months']), 'os5': num(r['os_5y']),
                      'os5_lo': num(ci[0]), 'os5_hi': num(ci[1])})
check('test high n', [g for g in km_groups if g['cohort'] == 'GSE181063 test' and g['group'] == 'high'][0]['n'], 90, 0)

km_path = R / '09c' / 'km_curves_09c.json'
km_curves = json.loads(km_path.read_text()) if km_path.exists() else None
km_figures = None
if km_curves is None:
    figs = {('GSE181063 test', 'median'): F / '05' / 'km_test_median_05.png',
            ('GSE181063 test', 'tertile'): F / '05' / 'km_test_tertiles_05.png',
            ('GSE10846', 'median'): F / '06' / 'km_GSE10846_median_06.png',
            ('GSE10846', 'tertile'): F / '06' / 'km_GSE10846_tertiles_06.png',
            ('GSE87371', 'median'): F / '06' / 'km_GSE87371_median_06.png',
            ('GSE87371', 'tertile'): F / '06' / 'km_GSE87371_tertiles_06.png'}
    km_figures = {f'{c}|{g}': 'data:image/png;base64,' + base64.b64encode(p.read_bytes()).decode()
                  for (c, g), p in figs.items()}

# ------------------------------------------------------------------ negative / leakage control
nc = json.loads((R / '05' / 'negative_control_05.json').read_text())
nc_c = pd.read_csv(R / '05' / 'negative_control_05.csv')['test_c'].tolist()
lk = pd.read_csv(R / '05' / 'leakage_check_05.csv')
negative = {'observed_c': nc['observed_c'], 'perm_c': nc_c, 'mean_c': nc['mean_c'], 'max_c': nc['max_c'],
            'p_empirical': nc['p_empirical'], 'n': nc['n_permutations'], 'share_empty': nc['share_empty_models'],
            'leakage_mean_c': float(lk['c'].mean()), 'leakage_min_c': float(lk['c'].min()),
            'leakage_max_c': float(lk['c'].max()), 'leakage_units': int(len(lk))}
check('perm mean C', negative['mean_c'], 0.493, 0.0005)
check('leakage mean C', negative['leakage_mean_c'], 0.499, 0.0005)

# ------------------------------------------------------------------ clinical comparison (step 7)
cc = json.loads((R / '07' / 'clinical_comparison_07.json').read_text())
cidx = pd.read_csv(R / '07' / 'c_index_models_07.csv')
hr7 = pd.read_csv(R / '07' / 'hr_lp_07.csv')
lrt = pd.read_csv(R / '07' / 'lrt_07.csv')
boot = pd.read_csv(R / '07' / 'delta_c_bootstrap_07.csv')
clinical = {'cohorts': []}
edges = np.round(np.arange(-0.04, 0.0801, 0.004), 4)
for key, label in COHORT_LABEL.items():
    sub = cidx[cidx['cohort'] == key].set_index('model')['harrell_c'].to_dict()
    pop = cc['populations'][key]
    b = boot[key].dropna().to_numpy()
    counts, _ = np.histogram(np.clip(b, edges[0], edges[-1] - 1e-9), bins=edges)
    h = hr7[hr7['cohort'] == key].set_index('model')
    l_ = lrt[lrt['cohort'] == key].set_index('comparison')
    clinical['cohorts'].append({
        'key': key, 'label': label, 'n_all': pop['n'], 'deaths_all': pop['deaths'],
        'missing_ipi': pop['missing_ipi'], 'n': int(cidx[cidx['cohort'] == key]['n'].iloc[0]),
        'deaths': int(cidx[cidx['cohort'] == key]['deaths'].iloc[0]), 'c': sub,
        'delta_c': sub['IPI + COO + LP'] - sub['IPI + COO'],
        'delta_c_ci': cc['delta_c_bootstrap'][key]['delta_c_ci_95'],
        'delta_c_hist': {'edges': edges.tolist(), 'counts': counts.tolist(), 'n': int(len(b))},
        'hr': {m: {'hr': h.loc[m, 'hr_per_sd'], 'lo': h.loc[m, 'ci_low'], 'hi': h.loc[m, 'ci_high']} for m in h.index},
        'lrt': {c: {'chi2': l_.loc[c, 'chi2'], 'p': l_.loc[c, 'p']} for c in l_.index},
        'spearman_lp_ipi': cc['overlap_descriptive'][key]['spearman_lp_ipi'],
    })
clinical['h3'] = cc['H3']
by = {c['key']: c for c in clinical['cohorts']}
check('H3 LRT p', by['GSE181063_test']['lrt']['IPI + COO vs IPI + COO + LP']['p'], 0.0067, 0.0001)
check('adj HR test', by['GSE181063_test']['hr']['IPI + COO + LP']['hr'], 1.63, 0.005)
check('adj HR GSE10846', by['GSE10846']['hr']['IPI + COO + LP']['hr'], 1.29, 0.005)
check('adj HR GSE87371', by['GSE87371']['hr']['IPI + COO + LP']['hr'], 1.17, 0.005)
check('dC test', by['GSE181063_test']['delta_c'], 0.020, 0.0006)
check('dC CI low test', by['GSE181063_test']['delta_c_ci'][0], -0.004, 0.0006)
check('IPI+COO C test', by['GSE181063_test']['c']['IPI + COO'], 0.727, 0.0006)

# ------------------------------------------------------------------ biology (step 8, exploratory)
rk = pd.read_csv(R / '08' / 'rank_test_08.csv')
hall = rk[rk['collection'] == 'HALLMARK'].copy()
hall['label'] = hall['gene_set'].str.replace('HALLMARK_', '', regex=False).str.replace('_', ' ').str.title()
hallmark = [{'set': r['label'], 'n': int(r['n_genes']), 'z': r['z_gene_sampling'], 'mean_rho': r['mean_rho_in_set'],
             'robust': bool(r['robust']), 'fdr_perm': r['fdr_bh'], 'fdr_gs': r['fdr_gene_sampling'],
             'leading': r['leading_genes'] if isinstance(r['leading_genes'], str) else ''}
            for _, r in hall.sort_values('z_gene_sampling', ascending=False).iterrows()]
check('Hallmark robust count', float(sum(h['robust'] for h in hallmark)), 11, 0)
check('E2F z', [h for h in hallmark if h['set'] == 'E2F Targets'][0]['z'], 12.3, 0.05)
robust_counts = rk.groupby('collection')['robust'].agg(['sum', 'size']).reset_index()
robust_counts = [{'collection': r['collection'], 'robust': int(r['sum']), 'total': int(r['size'])}
                 for _, r in robust_counts.iterrows()]

sig = pd.read_csv(R / '08' / 'signature_associations_08.csv')
signatures = [{'name': r['signature'], 'group': r['group'], 'rho': r['rho_with_lp'], 'lo': r['ci95_low'],
               'hi': r['ci95_high'], 'rho_removed': num(r['rho_model_genes_removed']), 'fdr': num(r['fdr_bh']),
               'n_genes': int(r['n_genes_found'])} for _, r in sig.iterrows()]
check('Lymph2Cx rho', [s for s in signatures if s['name'] == 'LYMPH2CX_COO_SCORE'][0]['rho'], 0.46, 0.005)

# ------------------------------------------------------------------ model genes (Table S2)
s2 = pd.read_csv(R / '09b' / 'table_s2_model_coefficients_09b.csv')
genes = [{'rank': int(r['rank']), 'entrez': str(r['entrez']), 'symbol': r['symbol_current'],
          'coef': r['coef_log_hr_per_train_sd'], 'hr': r['hr_per_sd'], 'freq': r['bootstrap_selection_freq']}
         for _, r in s2.iterrows()]
check('n genes', len(genes), 44, 0)
check('positive coefs', float(sum(g['coef'] > 0 for g in genes)), 26, 0)

# ------------------------------------------------------------------ assemble
failed = [c for c in checks if not c['ok']]
for c in checks:
    print(('PASS' if c['ok'] else 'FAIL'), c['check'], c['value'], 'expected', c['expected'])
if failed:
    raise SystemExit(f'{len(failed)} check(s) failed; site not built.')

flow_png = 'data:image/png;base64,' + base64.b64encode((F / '09b' / 'flow_diagram_09b.png').read_bytes()).decode()
data = clean({'meta': meta, 'flow_png': flow_png, 'cohorts': cohorts, 'nested_cv': nested, 'validation': validation,
              'km_groups': km_groups, 'km_curves': km_curves, 'km_figures': km_figures,
              'negative': negative, 'clinical': clinical, 'hallmark': hallmark,
              'robust_counts': robust_counts, 'signatures': signatures, 'genes': genes,
              'checks': {'n': len(checks), 'passed': len(checks) - len(failed)}})
payload = json.dumps(data, separators=(',', ':'), allow_nan=False).replace('</', '<\\/')
tpl = TEMPLATE.read_text(encoding='utf-8')
assert tpl.count('/*__SITE_DATA__*/null') == 1 and tpl.count('<!--BODY-->') == 1, 'template markers missing'
tpl = tpl.replace('/*__SITE_DATA__*/null', payload)
head, body = tpl.split('<!--BODY-->')
head = head.split('-->', 1)[1] if head.lstrip().startswith('<!--') else head   # drop the template note
html = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        + head.strip() + '\n</head>\n<body>\n' + body.strip() + '\n</body>\n</html>\n')
OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(html, encoding='utf-8')
(OUT.parent / '.nojekyll').write_text('')
print(f'Wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size / 1024:.0f} KB); '
      f'{len(checks)} checks passed; KM: {"interactive curves" if km_curves else "published figures"}')
