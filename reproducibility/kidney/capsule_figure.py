"""The glomerular capsules on the kidney H&E, and one capsule's 2 um spots split into two groups.

Rebuilds RECOMB Fig 2b, which was a hand composite of a matplotlib overlay and three QuPath renders
of unseeded splits. The split is spot_split's: per capsule, n_A ~ Binomial(n_spots, p) spots chosen
at random go to group A, drawn in yellow; the rest are group B.

    python -m kidney.capsule_figure      # from reproducibility/

Needs kidney/data/ from `python -m kidney.preprocess` and the full-resolution H&E,
`python -m fetch_data --arm kidney` (kidney-he-fullres, 4.3 GB).
"""

import argparse
import json

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import shapely
import shapely.affinity
import tifffile
import zarr
from matplotlib.collections import PolyCollection
from matplotlib.image import imread

import plot_style
from paths import data_dir, fig_name, require_input, results_dir

DATA = data_dir(__file__)
SPLIT_P = (0.5, 0.3, 0.1)
OUTLINE = '#ffd400'   # capsule outlines on the overview
BOUNDS = '#f2f2f2'    # a capsule's outer bounds in the close-ups
GROUP_A = '#ffd400'   # its group-A spots
BAND = 2.9            # overview width / height, as in the published strip


def load():
    obs = ad.read_h5ad(require_input(DATA / 'podocytes_2um.h5ad', what='the per-2um-spot capsule matrix',
                                     source='python -m kidney.preprocess'), backed='r').obs
    sf = json.load(open(require_input(DATA / 'binned_outputs/square_002um/spatial/scalefactors_json.json',
                                      what='the 2um bin scale factors', source='python -m fetch_data --arm kidney')))
    return obs[['shape_id', 'pxl_col_in_fullres', 'pxl_row_in_fullres', 'array_row', 'array_col']].to_numpy(), sf


def capsule_outlines(spots, half):
    """Each capsule's outline: the union of the 16 um bins its 2 um spots fall in, in full-resolution pixels.

    A 16 um bin is 8 x 8 spots. Each bin's box spans its spots' squares, grown by half a pixel so that
    neighbouring bins overlap and merge instead of leaving hairline gaps.
    """
    outlines = {}
    for sid in np.unique(spots[:, 0]):
        mine = spots[spots[:, 0] == sid]
        bins = (mine[:, 3] // 8) * 100000 + mine[:, 4] // 8
        boxes = []
        for b in np.unique(bins):
            xy = mine[bins == b, 1:3]
            boxes.append(shapely.box(xy[:, 0].min() - half - 0.5, xy[:, 1].min() - half - 0.5,
                                     xy[:, 0].max() + half + 0.5, xy[:, 1].max() + half + 0.5))
        outlines[sid] = shapely.union_all(boxes)
    return outlines


def draw_outline(ax, geom, scale, **kw):
    for poly in getattr(geom, 'geoms', [geom]):
        x, y = poly.exterior.xy
        ax.plot(np.asarray(x) * scale, np.asarray(y) * scale, **kw)


def scale_bar(ax, microns, px_per_micron, label):
    """A white bar in the lower left, in place of pixel-coordinate ticks."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()  # images run top to bottom, so y0 is the bottom
    w, h = x1 - x0, y0 - y1
    left, bottom = x0 + 0.04 * w, y0 - 0.05 * h
    ax.plot([left, left + microns * px_per_micron], [bottom, bottom], color='white', lw=2, solid_capstyle='butt')
    ax.text(left, bottom - 0.025 * h, label, color='white', fontsize=plot_style.SMALL, va='bottom')


def read_crop(tif, x0, y0, x1, y1):
    """A region of the full-resolution image, decoding only the tiles it covers."""
    store = tifffile.imread(tif, aszarr=True)
    try:
        z = zarr.open(store, mode='r')
        full = z['0'] if isinstance(z, zarr.Group) else z  # level 0 of a pyramid is full resolution
        return np.asarray(full[y0:y1, x0:x1])
    finally:
        store.close()


def densest_band(centres, width, height):
    """Top of the horizontal band of this height holding the most capsule centres."""
    ys = np.sort(centres[:, 1])
    counts = [np.sum((ys >= y) & (ys < y + height)) for y in ys]
    return ys[int(np.argmax(counts))]


def main(capsule, seed):
    spots, sf = load()
    mpp = sf['microns_per_pixel']
    half = sf['bin_size_um'] / mpp / 2
    outlines = capsule_outlines(spots, half)

    if capsule is None:  # the capsule with the median number of spots
        sizes = {sid: np.sum(spots[:, 0] == sid) for sid in outlines}
        capsule = min(sizes, key=lambda s: abs(sizes[s] - np.median(list(sizes.values()))))
    mine = spots[spots[:, 0] == capsule, 1:3]

    fig = plt.figure(figsize=(7.2, 7.2 / BAND + 2.6), layout='constrained')
    top, bottom = fig.subfigures(2, 1, height_ratios=[7.2 / BAND, 2.6])
    ax = top.subplots()
    axes = bottom.subplots(1, len(SPLIT_P))

    # Overview: the hires image, cropped to the band of the section with the most capsules
    hires = imread(require_input(DATA / 'spatial/tissue_hires_image.png', what='the hires H&E',
                                 source='python -m fetch_data --arm kidney'))
    s = sf['tissue_hires_scalef']
    centres = np.array([[g.centroid.x, g.centroid.y] for g in outlines.values()])
    x0, x1 = centres[:, 0].min() - 600, centres[:, 0].max() + 600
    height = (x1 - x0) / BAND
    y0 = densest_band(centres, x1 - x0, height) - 0.1 * height
    r0, r1, c0, c1 = (int(v) for v in (y0 * s, (y0 + height) * s, x0 * s, x1 * s))
    ax.imshow(hires[max(r0, 0):r1, max(c0, 0):c1], extent=(max(c0, 0), c1, r1, max(r0, 0)),
              interpolation='antialiased')
    for geom in outlines.values():
        draw_outline(ax, geom, s, color=OUTLINE, lw=0.6)
    ax.set_xlim(c0, c1)
    ax.set_ylim(r1, r0)
    ax.set_axis_off()
    scale_bar(ax, 500, s / mpp, '500 µm')

    # Close-ups: one capsule on the full-resolution H&E, group A in yellow
    g = outlines[capsule]
    bx0, by0, bx1, by1 = g.bounds
    pad = 0.12 * max(bx1 - bx0, by1 - by0)
    cx0, cy0, cx1, cy1 = (int(v) for v in (bx0 - pad, by0 - pad, bx1 + pad, by1 + pad))
    tif = require_input(DATA / 'Visium_HD_Human_Kidney_FFPE_tissue_image.tif', what='the full-resolution H&E',
                        source='python -m fetch_data --arm kidney')
    crop = read_crop(tif, cx0, cy0, cx1, cy1)
    rng = np.random.default_rng(seed)
    for ax, p in zip(axes, SPLIT_P):
        n_a = rng.binomial(len(mine), p)
        a = mine[rng.choice(len(mine), n_a, replace=False)] - [cx0, cy0]
        ax.imshow(crop, extent=(0, cx1 - cx0, cy1 - cy0, 0), interpolation='antialiased')
        squares = [np.array([[x - half, y - half], [x + half, y - half], [x + half, y + half], [x - half, y + half]])
                   for x, y in a]
        ax.add_collection(PolyCollection(squares, facecolors=(1, 0.83, 0, 0.25), edgecolors=GROUP_A, linewidths=0.4))
        draw_outline(ax, shapely.affinity.translate(g, -cx0, -cy0), 1, color=BOUNDS, lw=0.8)
        ax.set_title(f'$p$ = {p:g}')
        ax.set_axis_off()
    scale_bar(axes[0], 50, 1 / mpp, '50 µm')

    out = results_dir(__file__) / fig_name('capsule_split_illustration')
    fig.savefig(out, dpi=plot_style.RASTER_DPI, bbox_inches='tight')
    plt.close(fig)
    print(f'wrote {out} (capsule {capsule:g}, {len(mine)} spots)')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--capsule', type=float, default=None, help='shape_id of the capsule to split; default the median-sized one')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()
    plot_style.use()
    main(args.capsule, args.seed)
