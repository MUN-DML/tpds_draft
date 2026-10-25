"""Fig. 5: measured absolute FCT; unavailable runs are None, never synthesized."""
from pathlib import Path
import hashlib
import json
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/tpds-fig5-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter, FuncFormatter
import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
RUNS = Path('/home/hongw/Desktop/ns-allinone-3.19/prolet/output')
SCHEMES = ['DuoGran', 'ConWeave', 'CONGA']
WORKLOADS = {
    'AliStorage': {
        'flows': 156777,
        'edges': [0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 250000, 500000, 750000, 1000000, 1500000, 2000000],
        'files': {
            'DuoGran': RUNS / 'duogran-default-subflow-20261019/alistorage/ewma_q/run_out_fct.txt',
            'ConWeave': RUNS / 'duogran-conga-conweave-20261003/conweave/run_out_fct.txt',
            'CONGA': RUNS / 'duogran-conga-conweave-20261003/conga/run_out_fct.txt',
        },
    },
    'Meta Hadoop': {
        'flows': 52469,
        'edges': [0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 250000, 1000000, 2000000, 4000000, 10000000],
        'files': {
            'DuoGran': RUNS / 'duogran-default-subflow-20261019/hadoop/ewma_q/run_out_fct.txt',
            'ConWeave': None,
            'CONGA': None,
        },
    },
}
STYLES = {
    'DuoGran': dict(color='#D55E00', marker='s', linestyle='--', linewidth=1.5),
    'ConWeave': dict(color='#0072B2', marker='o', linestyle='-', linewidth=1.2),
    'CONGA': dict(color='#009E73', marker='^', linestyle='-.', linewidth=1.2),
}
ALI_AXES = {
    'mean_ms': dict(knee=30.0, compression=150.0, lower=10.0, upper=1250.0,
                    reference_compression=4.5, reference_upper=1500.0,
                    ticks=[10, 20, 30, 300, 1000]),
    'p99_ms': dict(knee=120.0, compression=400.0, lower=40.0, upper=2500.0,
                   reference_compression=12.0, reference_upper=3500.0,
                   ticks=[40, 80, 120, 1000]),
}
plt.rcParams.update({
    'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 9.5,
    'axes.labelsize': 9.5, 'axes.titlesize': 9.5,
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'axes.linewidth': .6, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'savefig.facecolor': 'white',
})
def linear_log_scale(knee, compression, upper, reference_compression, reference_upper):
    """Keep the approved low-range positions and normalize only the upper range."""
    reference_top = knee + reference_compression * np.log1p(
        (reference_upper - knee) / reference_compression)
    gain = (reference_top - knee) / np.log1p((upper - knee) / compression)

    def forward(values):
        y = np.asarray(values, dtype=float)
        tail = knee + gain * np.log1p(np.maximum(y - knee, 0) / compression)
        return np.where(y <= knee, y, tail)

    def inverse(values):
        z = np.asarray(values, dtype=float)
        tail = knee + compression * np.expm1(np.maximum(z - knee, 0) / gain)
        return np.where(z <= knee, z, tail)

    return forward, inverse


summary = {
    'metric': 'absolute FCT in milliseconds',
    'settings': '128-host Spine-Leaf, 100 Gbps, OS2, DCQCN, PFC, 4096-byte payload, 80% network load, seed 1, 10 ms injection',
    'cohort': 'all completed flows including drain; no warmup exclusion',
    'grouping': 'fixed workload-specific byte bins; lower exclusive, upper inclusive; shared by all schemes',
    'quantile': 'numpy linear interpolation, p=0.99; computed from individual FCTs in each bin',
    'missing': 'None means no matching measurements; no curve is drawn',
    'axes': 'logarithmic x at actual bin upper bounds (KB); AliStorage piecewise linear/logarithmic y; Meta Hadoop logarithmic y',
    'alistorage_y_transform': {'formula': 'f(y)=y for y<=k; f(y)=k+A*ln(1+(y-k)/s) for y>k; y in microseconds',
                               'normalization': 'A fixes the displayed height of the low linear range to the approved reference scale.',
                               'all_measured_values_visible': False,
                               'clipping': 'Values above the displayed limits are clipped; raw statistics remain unchanged.',
                               'clipped_points': []},
    'workloads': {},
}
fig, axes = plt.subplots(2, 2, figsize=(7.16, 3.0))
fig.subplots_adjust(left=.085, right=.975, top=.89, bottom=.19, hspace=1.25, wspace=.27)
for row, (name, cfg) in enumerate(WORKLOADS.items()):
    rows = {}
    metadata = {}
    for scheme, path in cfg['files'].items():
        if path is None:
            metadata[scheme] = None
            continue
        data = np.loadtxt(path, dtype=np.int64)
        assert data.shape == (cfg['flows'], 8)
        order = np.lexsort(tuple(data[:, i] for i in (5, 3, 2, 1, 0, 4)))
        data = data[order]
        assert len(np.unique(data[:, :6], axis=0)) == len(data)
        assert np.all(data[:, 4] > 0) and np.all(data[:, 6:] > 0)
        rows[scheme] = data
        fct = data[:, 6] / 1e6
        metadata[scheme] = {
            'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'flow_count': len(data), 'mean_ms': float(fct.mean()),
            'p99_ms': float(np.quantile(fct, .99, method='linear')),
        }
    reference = rows['DuoGran']
    for data in rows.values():
        assert np.array_equal(data[:, :6], reference[:, :6]), 'mismatched flow identities'
        assert np.array_equal(data[:, 7], reference[:, 7]), 'mismatched ideal FCT'
    sizes = reference[:, 4]
    bins = []
    masks = []
    for lo, hi in zip(cfg['edges'][:-1], cfg['edges'][1:]):
        mask = (sizes > lo) & (sizes <= hi)
        assert mask.any()
        masks.append(mask)
        stats = {}
        for scheme in SCHEMES:
            if scheme not in rows:
                stats[scheme] = None
                continue
            fct = rows[scheme][mask, 6] / 1e6
            stats[scheme] = {'mean_ms': float(fct.mean()),
                             'p99_ms': float(np.quantile(fct, .99, method='linear'))}
        bins.append({'lower_exclusive_bytes': lo, 'upper_inclusive_bytes': hi,
                     'n': int(mask.sum()), 'statistics': stats})
    assert np.all(np.sum(masks, axis=0) == 1), 'bins must cover every flow exactly once'
    summary['workloads'][name] = {'inputs': metadata, 'bins': bins}
    x = np.array(cfg['edges'][1:], dtype=float) / 1000
    for col, metric in enumerate(['mean_ms', 'p99_ms']):
        ax = axes[row, col]
        panel_pos = ax.get_position()
        for scheme in ['CONGA', 'ConWeave', 'DuoGran']:
            if scheme not in rows:
                continue
            y = [1000 * b['statistics'][scheme][metric] for b in bins]
            ax.plot(x, y, **STYLES[scheme], markersize=3.6,
                    markerfacecolor='white', markeredgewidth=.85, label=scheme)
        fig.text(panel_pos.x0 + panel_pos.width/2, panel_pos.y0 - .65*panel_pos.height,
                 f'({chr(97+row*2+col)}) {name}', ha='center', va='top', fontsize=9.5)
        ax.set_ylabel('Average FCT (μs)' if col == 0 else 'P99 FCT (μs)', labelpad=4)
        ax.set_xlabel('Flow size upper bound (KB)', labelpad=3)
        ax.set_xscale('log')
        ax.set_yscale('log')
        xticks = [1, 10, 100, 1000] if row == 0 else [1, 10, 100, 1000, 10000]
        ax.set_xticks(xticks)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_xlim(x[0] / 1.20, x[-1] * 1.45)
        plotted = [1000 * b['statistics'][scheme][metric]
                   for b in bins for scheme in SCHEMES if scheme in rows]
        if row == 0:
            axis_cfg = ALI_AXES[metric]
            transform = linear_log_scale(axis_cfg['knee'], axis_cfg['compression'],
                                         axis_cfg['upper'], axis_cfg['reference_compression'],
                                         axis_cfg['reference_upper'])
            reference_top = axis_cfg['knee'] + axis_cfg['reference_compression'] * np.log1p(
                (axis_cfg['reference_upper'] - axis_cfg['knee']) / axis_cfg['reference_compression'])
            assert np.isclose(transform[0](axis_cfg['upper']), reference_top)
            low_values = np.linspace(axis_cfg['lower'], axis_cfg['knee'], 100)
            assert np.array_equal(transform[0](low_values), low_values)
            summary['alistorage_y_transform']['average' if col == 0 else 'p99'] = {
                'k': axis_cfg['knee'], 's': axis_cfg['compression'],
                'A': float((reference_top - axis_cfg['knee']) / np.log1p(
                    (axis_cfg['upper'] - axis_cfg['knee']) / axis_cfg['compression'])),
                'limits_us': [axis_cfg['lower'], axis_cfg['upper']],
                'low_range_height_fraction': float((axis_cfg['knee'] - axis_cfg['lower']) /
                                                   (reference_top - axis_cfg['lower'])),
                'reference_compression': axis_cfg['reference_compression'],
                'reference_upper_us': axis_cfg['reference_upper'],
            }
            ax.set_yscale('function', functions=transform)
            ax.set_yticks(axis_cfg['ticks'])
            ax.set_ylim(axis_cfg['lower'], axis_cfg['upper'])
            for scheme in SCHEMES:
                values = [1000 * b['statistics'][scheme][metric] for b in bins]
                if scheme in ('DuoGran', 'ConWeave'):
                    assert min(values) >= axis_cfg['lower'] and max(values) <= axis_cfg['upper']
                for upper_bytes, value in zip(cfg['edges'][1:], values):
                    if value > axis_cfg['upper']:
                        summary['alistorage_y_transform']['clipped_points'].append({
                            'scheme': scheme, 'metric': metric,
                            'bin_upper_bytes': upper_bytes, 'value_us': value,
                        })
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
            ax.yaxis.set_minor_formatter(NullFormatter())
            check = np.geomspace(min(plotted), max(plotted), 1000)
            assert np.all(np.diff(transform[0](check)) > 0)
            assert np.allclose(transform[1](transform[0](check)), check)
        else:
            ax.set_ylim(min(plotted) / 6, max(plotted) * 1.18)
            ax.set_yticks([5, 20, 200, 2000] if col == 0 else [20, 100, 1000, 10000])
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(2, 5)))
            ax.yaxis.set_minor_formatter(NullFormatter())
        ax.grid(axis='y', which='major', color='#cccccc', linewidth=.5)
        ax.set_axisbelow(True)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(length=2.5, width=.6, pad=2)
handles = [Line2D([0], [0], **STYLES[s], markersize=4, markerfacecolor='white', label=s) for s in SCHEMES]
fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.53, 1.0), ncol=3,
           frameon=False, handlelength=2.8, columnspacing=2.5, fontsize=9.5)
for ext in ['png', 'pdf']:
    fig.savefig(OUT / f'workload_fct_comparison.{ext}', dpi=250)
fig.savefig(ROOT / 'Source/figures/workload_fct_comparison.pdf')
plt.close(fig)
(OUT / 'plot_data.json').write_text(json.dumps(summary, indent=2)+'\n')
print('Verified same-flow identities and ideal FCTs, all-flow coverage, no fabricated missing curves.')
for name, data in summary['workloads'].items():
    print(name, {s: None if v is None else v['flow_count'] for s,v in data['inputs'].items()})
