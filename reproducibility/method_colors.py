"""One colour per DE method, shared by every arm's figures.

The arms label the same method differently ("DELN", "LN test", "LN's $t$-test";
"t-test", "Scanpy t-test (log1p)"; ...), so `method_color` accepts any of those
labels. An unknown label raises instead of falling back, so a new method cannot
quietly get an off-palette colour.
"""

from matplotlib.patheffects import withStroke

LN = '#082f21'            # dark green, close to black: black overlays vanish on it
LOG1P = '#8c564b'         # tab10 brown
WILCOXON = '#17becf'      # tab10 cyan
MAST = '#ff7f0e'          # tab10 orange
LOG1P_ACTUAL = '#c49c94'  # tab20 light brown: the t-test's LFC recomputed on the log1p scale

# One hue, light to dark, for an ordered variable within the log1p t-test's series.
LOG1P_RAMP = ['#c49c94', '#a8776b', LOG1P, '#5e3a32']

# Box-plot medians that read on every fill, from LN's near-black to the light brown.
MEDIANPROPS = {'color': 'white', 'linewidth': 1.5,
               'path_effects': [withStroke(linewidth=3, foreground='black')]}

_COLORS = {
    **dict.fromkeys(['ln', 'deln', 'ln_test', 'ln test', "ln's test", "ln's t-test",
                     "ln's-t-test", 'ziln', 'lognormal'], LN),
    **dict.fromkeys(['log1p', 'log1p t-test', 't-test', 't_test', 'scanpy',
                     'scanpy t-test', 'scanpy t-test (log1p)'], LOG1P),
    **dict.fromkeys(['wilcoxon', 'scanpy wilcoxon', 'scanpy wilcoxon (log1p)'], WILCOXON),
    'mast': MAST,
    'scanpy t-test (actual lfc)': LOG1P_ACTUAL,
}


def method_color(label) -> str:
    """Colour for a method label, ignoring case and LaTeX dollar signs."""
    key = str(label).replace('$', '').strip().lower()
    try:
        return _COLORS[key]
    except KeyError:
        raise KeyError(f'no standard colour for method {label!r}; add it to method_colors.py') from None
