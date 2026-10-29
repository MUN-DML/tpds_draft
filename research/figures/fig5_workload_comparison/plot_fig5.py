"""Complete ConWeave/CONGA baselines with the original DuoGran curves retained."""
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
ARCHIVE = Path(os.environ.get('PROLET_LOG_ARCHIVE',
    ROOT.parent / 'ns-allinone-3.19/prolet-experiment-logs'))
SCHEMES = ['DuoGran', 'ConWeave', 'CONGA']
WORKLOADS = {
    'AliStorage': {
        'flows': {'DuoGran': 156777, 'ConWeave': 1567425, 'CONGA': 1567425},
        'edges': [0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 250000, 500000, 750000, 1000000, 1500000, 2000000],
        'records': {
            'DuoGran': 'runs/AliStorage2019/DuoGran/output-duogran-owl-admission-20261020-alistorage-before-07fa4cab/record.json',
            'ConWeave': 'runs/AliStorage2019/ConWeave/reports-matrix-cdf-100ms-20261029-conweave-alistorage-ccffb81b/record.json',
            'CONGA': 'runs/AliStorage2019/CONGA/output-conga-paper-128-100ms-20261023-4a42b096/record.json',
        },
    },
    'Meta Hadoop': {
        'flows': {'DuoGran': 52469, 'ConWeave': 525310, 'CONGA': 525310},
        'edges': [0, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 250000, 1000000, 2000000, 4000000, 10000000],
        'records': {
            'DuoGran': 'runs/MetaHadoop/DuoGran/output-duogran-owl-admission-20261020-hadoop-before-96bf77b0/record.json',
            'ConWeave': 'runs/MetaHadoop/ConWeave/reports-matrix-cdf-100ms-20261029-conweave-metahadoop-70a92184/record.json',
            'CONGA': 'runs/MetaHadoop/CONGA/reports-matrix-cdf-100ms-20261029-conga-metahadoop-5caf6ceb/record.json',
        },
    },
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def resolve_record(relative, scheme, flow_count, arrival_window_ms):
    record = json.loads((ARCHIVE / relative).read_text())
    assert record['scheme'] == scheme and record['status'] == 'complete'
    assert record['completed_flows'] == record['expected_flows'] == flow_count
    assert record['hosts'] == 128 and record['network_load_percent'] == 80
    assert record['arrival_window_ms'] == arrival_window_ms and record['seed'] == 1
    c = record['configuration']
    expected = {'WORKLOAD_MODE': 'cdf', 'CC_MODE': '1', 'ENABLE_PFC': '1',
                'ENABLE_IRN': '0', 'PACKET_PAYLOAD_SIZE': '4096',
                'HAS_WIN': '0', 'LOAD': '80', 'RANDOM_SEED': '1'}
    assert all(c.get(k) == v for k, v in expected.items())
    assert record['topology_sha256'] == '0dddc4f3ae673139b895ff2f875befe82b234cc48e325bd8a164d9cddedc12e2'
    artifacts = {a.get('config_key'): a for a in record['artifacts'] if a.get('config_key')}
    for key in ('TOPOLOGY_FILE', 'FLOW_FILE', 'FCT_OUTPUT_FILE'):
        assert digest(ARCHIVE / artifacts[key]['archive_path']) == artifacts[key]['sha256']
    assert artifacts['FLOW_FILE']['sha256'] == record['input_sha256']
    return record, ARCHIVE / artifacts['FCT_OUTPUT_FILE']['archive_path']


STYLES = {
    'DuoGran': dict(color='#D55E00', marker='s', linestyle='--', linewidth=1.5),
    'ConWeave': dict(color='#0072B2', marker='o', linestyle='-', linewidth=1.2),
    'CONGA': dict(color='#009E73', marker='^', linestyle='-.', linewidth=1.2),
}
ALI_AXES = {
    'mean_ms': dict(knee=30.0, compression=3.0, lower=10.0, upper=2200.0,
                    ticks=[10, 20, 30, 100, 1000]),
    'p99_ms': dict(knee=120.0, compression=40.0, lower=40.0, upper=4000.0,
                   ticks=[40, 80, 120, 500, 3000]),
}
META_AXES = {
    'mean_ms': dict(knee=20.0, compression=1.5, lower=10.0, upper=8500.0,
                    ticks=[10, 15, 20, 100, 5000]),
    'p99_ms': dict(knee=100.0, compression=4.0, lower=40.0, upper=40000.0,
                   ticks=[40, 70, 100, 500, 20000]),
}
plt.rcParams.update({
    'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix', 'font.size': 9.5,
    'axes.labelsize': 9.5, 'axes.titlesize': 9.5,
    'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'axes.linewidth': .6, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'savefig.facecolor': 'white',
})
def linear_log_scale(knee, compression):
    """Linear lower range and a smooth, progressively compressed upper range.

    Both the value and derivative agree at the knee. The former independent
    tail normalization sharply reduced the slope immediately above the knee.
    """
    def forward(values):
        y = np.asarray(values, dtype=float)
        tail = knee + compression * np.log1p(np.maximum(y - knee, 0) / compression)
        return np.where(y <= knee, y, tail)

    def inverse(values):
        z = np.asarray(values, dtype=float)
        tail = knee + compression * np.expm1(np.maximum(z - knee, 0) / compression)
        return np.where(z <= knee, z, tail)

    return forward, inverse


summary = {
    'metric': 'absolute FCT in milliseconds',
    'settings': '128-host Spine-Leaf, 100 Gbps, OS2, DCQCN, PFC, 4096-byte payload, 80% network load, seed 1; injection windows are recorded per scheme',
    'cohort': 'all completed flows including drain; no warmup exclusion',
    'grouping': 'fixed workload-specific byte bins; lower exclusive, upper inclusive; shared by all schemes',
    'quantile': 'numpy linear interpolation, p=0.99; computed from individual FCTs in each bin',
    'cohorts': 'ConWeave and CONGA have paired 100-ms inputs; original DuoGran 10-ms curves are retained as historical results, with independent flow counts per bin',
    'axes': 'logarithmic x at actual bin upper bounds (KB); both workloads use a linear lower y range and progressively logarithmic upper-range compression',
    'alistorage_y_transform': {'formula': 'f(y)=y for y<=k; f(y)=k+s*ln(1+(y-k)/s) for y>k; y in microseconds',
                               'normalization': 'The upper branch has unit slope at the knee, then compresses progressively; ConWeave has at least 5% displayed vertical headroom.',
                               'all_measured_values_visible': False,
                               'clipping': 'ConWeave and DuoGran are entirely inside the axes. CONGA may exceed the upper boundary; raw statistics are unchanged.',
                               'clipped_points': []},
    'workloads': {},
}
summary['metahadoop_y_transform'] = dict(summary['alistorage_y_transform'], clipped_points=[])
fig, axes = plt.subplots(2, 2, figsize=(7.16, 3.0))
fig.subplots_adjust(left=.085, right=.975, top=.89, bottom=.19, hspace=1.25, wspace=.27)
for row, (name, cfg) in enumerate(WORKLOADS.items()):
    rows = {}
    metadata = {}
    records = {}
    for scheme, relative in cfg['records'].items():
        window = 10 if scheme == 'DuoGran' else 100
        record, path = resolve_record(relative, scheme, cfg['flows'][scheme], window)
        records[scheme] = record
        data = np.loadtxt(path, dtype=np.int64)
        assert data.shape == (cfg['flows'][scheme], 8)
        order = np.lexsort(tuple(data[:, i] for i in (5, 3, 2, 1, 0, 4)))
        data = data[order]
        assert len(np.unique(data[:, :6], axis=0)) == len(data)
        assert np.all(data[:, 4] > 0) and np.all(data[:, 6:] > 0)
        rows[scheme] = data
        fct = data[:, 6] / 1e6
        metadata[scheme] = {
            'path': str(path), 'sha256': digest(path),
            'archive_record': relative, 'input_sha256': record['input_sha256'],
            'arrival_window_ms': window, 'network_load_percent': 80,
            'flow_count': len(data), 'mean_ms': float(fct.mean()),
            'p99_ms': float(np.quantile(fct, .99, method='linear')),
        }
    assert len({records[s]['input_sha256'] for s in ('ConWeave', 'CONGA')}) == 1, 'mismatched baseline traffic inputs'
    reference = rows['ConWeave']
    for data in (rows[s] for s in ('ConWeave', 'CONGA')):
        assert np.array_equal(data[:, :6], reference[:, :6]), 'mismatched flow identities'
        assert np.array_equal(data[:, 7], reference[:, 7]), 'mismatched ideal FCT'
    bins = []
    covered = {s: np.zeros(len(data), dtype=np.int8) for s, data in rows.items()}
    for lo, hi in zip(cfg['edges'][:-1], cfg['edges'][1:]):
        stats = {}
        counts = {}
        for scheme in SCHEMES:
            data = rows[scheme]
            mask = (data[:, 4] > lo) & (data[:, 4] <= hi)
            assert mask.any()
            covered[scheme] += mask
            fct = data[mask, 6] / 1e6
            counts[scheme] = int(mask.sum())
            stats[scheme] = {'mean_ms': float(fct.mean()),
                             'p99_ms': float(np.quantile(fct, .99, method='linear'))}
        assert counts['ConWeave'] == counts['CONGA']
        bins.append({'lower_exclusive_bytes': lo, 'upper_inclusive_bytes': hi,
                     'counts': counts, 'statistics': stats})
    assert all(np.all(v == 1) for v in covered.values()), 'bins must cover every flow exactly once'
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
                 f'({chr(97+row*2+col)}) {"Average FCT" if col == 0 else "P99 FCT"} comparison in {name}',
                 ha='center', va='top', fontsize=9.5)
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
        axis_cfg = (ALI_AXES if row == 0 else META_AXES)[metric]
        transform_metadata = summary['alistorage_y_transform' if row == 0 else 'metahadoop_y_transform']
        transform = linear_log_scale(axis_cfg['knee'], axis_cfg['compression'])
        low_values = np.linspace(axis_cfg['lower'], axis_cfg['knee'], 100)
        assert np.array_equal(transform[0](low_values), low_values)
        epsilon = 1e-5
        right_slope = (transform[0](axis_cfg['knee'] + epsilon) - axis_cfg['knee']) / epsilon
        assert np.isclose(right_slope, 1.0, atol=1e-5)
        transform_metadata['average' if col == 0 else 'p99'] = {
            'k': axis_cfg['knee'], 's': axis_cfg['compression'],
            'A': axis_cfg['compression'],
            'limits_us': [axis_cfg['lower'], axis_cfg['upper']],
            'low_range_height_fraction': float((axis_cfg['knee'] - axis_cfg['lower']) /
                (transform[0](axis_cfg['upper']) - axis_cfg['lower'])),
            'slope_at_knee': 1.0,
        }
        conweave_max = max(b['statistics']['ConWeave'][metric] * 1000 for b in bins)
        display_span = transform[0](axis_cfg['upper']) - transform[0](axis_cfg['lower'])
        headroom = float((transform[0](axis_cfg['upper']) - transform[0](conweave_max)) / display_span)
        assert headroom >= .05, 'ConWeave endpoint/marker needs visible headroom'
        transform_metadata['average' if col == 0 else 'p99']['conweave_top_margin_fraction'] = headroom
        ax.set_yscale('function', functions=transform)
        ax.set_yticks(axis_cfg['ticks'])
        ax.set_ylim(axis_cfg['lower'], axis_cfg['upper'])
        for scheme in SCHEMES:
            values = [1000 * b['statistics'][scheme][metric] for b in bins]
            assert min(values) >= axis_cfg['lower']
            if scheme in ('ConWeave', 'DuoGran'):
                assert max(values) <= axis_cfg['upper']
            for upper_bytes, value in zip(cfg['edges'][1:], values):
                if value > axis_cfg['upper']:
                    transform_metadata['clipped_points'].append({
                        'scheme': scheme, 'metric': metric,
                        'bin_upper_bytes': upper_bytes, 'value_us': value,
                    })
        ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:g}'))
        ax.yaxis.set_minor_formatter(NullFormatter())
        check = np.geomspace(min(plotted), max(plotted), 1000)
        assert np.all(np.diff(transform[0](check)) > 0)
        assert np.allclose(transform[1](transform[0](check)), check)
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
print('Verified raw checksums, paired baseline flow identities/ideal FCTs, independent DuoGran bins, and full visible flow coverage.')
for name, data in summary['workloads'].items():
    print(name, {s: None if v is None else v['flow_count'] for s,v in data['inputs'].items()})
