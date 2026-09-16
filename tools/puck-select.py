import os
import sys
import json
import shutil
import atexit
import tempfile
import argparse
import subprocess
import numpy as np
# matplotlib is imported lazily below: --non_interactive exits before the figure is built,
# so terra/reconstruction.wdl can replay a selection without a plotting stack

def get_args():
    parser = argparse.ArgumentParser(description='Interactively select a disk or polygon of beads from a Puck.csv')
    parser.add_argument("puck", type=str, help="path to Puck.csv (headerless sb,x,y), local or gs://")
    parser.add_argument("-s", "--selection", type=str, default=None, help="name this selection; required for a gs:// puck, whose sidecar is uploaded to <recon_dir>/<selection>/selection.json")
    parser.add_argument("-o", "--out", type=str, default=None, help="output csv (default: Puck-selected.csv beside the input)")
    parser.add_argument("-c", "--center", type=float, nargs=2, default=None, metavar=('X', 'Y'), help="pre-seed the circle center")
    parser.add_argument("-r", "--radius", type=float, default=None, help="pre-seed the circle radius")
    parser.add_argument("-p", "--vertices", type=float, nargs='+', default=None, metavar='V', help="polygon vertices as X1 Y1 X2 Y2 ... (at least 3 pairs); overrides --center/--radius")
    parser.add_argument("-j", "--from_json", type=str, default=None, help="replay the cut recorded in a selection.json sidecar, whatever shape it holds")
    parser.add_argument("-v", "--invert", action='store_true', help="select points outside the shape")
    parser.add_argument("-m", "--marker", type=str, default=",", help="matplotlib marker for the beads")
    parser.add_argument("-n", "--non_interactive", action='store_true', help="write the selection and exit, no window")

    args, unknown = parser.parse_known_args()
    [print(f"WARNING: unknown command-line argument {u}") for u in unknown]
    return args

# Load arguments
args = get_args()
puck = args.puck                       ; print(f"puck = {puck}")
selection = args.selection             ; print(f"selection = {selection}")
out = args.out                         ; print(f"out = {out}")
center = args.center                   ; print(f"center = {center}")
radius = args.radius                   ; print(f"radius = {radius}")
vertices = args.vertices               ; print(f"vertices = {vertices}")
from_json = args.from_json             ; print(f"from_json = {from_json}")
invert = args.invert                   ; print(f"invert = {invert}")
marker = args.marker                   ; print(f"marker = {marker}")
non_interactive = args.non_interactive  ; print(f"non_interactive = {non_interactive}")

# A cut is a dict: {"shape": "circle", "center": [x, y], "radius": r, "invert": bool}
#                  {"shape": "polygon", "vertices": [[x, y], ...],  "invert": bool}
# Sidecars written before polygons existed have no "shape" key, so it defaults to circle.
replay_source_puck = None
if from_json is not None:
    with open(from_json) as f:
        d = json.load(f)
    shape = {"shape": d.get("shape", "circle"), "invert": bool(d.get("invert", False))}
    if shape["shape"] == "circle":
        shape["center"], shape["radius"] = list(d["center"]), float(d["radius"])
    else:
        shape["vertices"] = [[float(vx), float(vy)] for vx, vy in d["vertices"]]
    # Keep the recorded provenance rather than re-deriving it from the local puck path,
    # which under terra/reconstruction.wdl is just the scratch dir the file was copied to
    replay_source_puck = d.get("source_puck")
    print(f"replaying {shape['shape']} from {from_json}")
elif vertices is not None:
    assert len(vertices) % 2 == 0, "--vertices takes X Y pairs, got an odd number of values"
    assert len(vertices) >= 6, "--vertices needs at least 3 points to close a shape"
    shape = {"shape": "polygon",
             "vertices": [[vertices[i], vertices[i + 1]] for i in range(0, len(vertices), 2)],
             "invert": invert}
else:
    shape = {"shape": "circle", "center": center, "radius": radius, "invert": invert}

# A gs:// puck is pulled to a temp dir and the sidecar is pushed back beside the selection.
# recon.py writes Puck.csv inside its own UMAP..._ne2000 folder, so the path names source_puck:
#   gs://<bucket>/recon/<bcl>/<index>/<source_puck>/Puck.csv
#   gs://<bucket>/recon/<bcl>/<index>/<selection>/selection.json
sidecar_uri = None
if puck.startswith("gs://"):
    assert selection, "a gs:// puck requires --selection"
    recon_uri, source_puck, puck_name = puck.rsplit("/", 2)
    sidecar_uri = f"{recon_uri}/{selection}/selection.json"
    tmpdir = tempfile.mkdtemp(prefix="puck-select-")
    atexit.register(shutil.rmtree, tmpdir, True)
    local = os.path.join(tmpdir, puck_name)
    print(f"downloading {puck} ...")
    subprocess.run(["gcloud", "storage", "cp", puck, local], check=True)
    puck = local
    if out is None:
        out = os.path.abspath(f"{selection}-selected.csv")
else:
    source_puck = os.path.basename(os.path.dirname(os.path.abspath(puck)))

assert os.path.isfile(puck), f"{puck} not found"
if out is None:
    out = os.path.join(os.path.dirname(os.path.abspath(puck)), 'Puck-selected.csv')
sidecar = os.path.splitext(out)[0] + '.json'
if non_interactive:
    assert from_json is not None or vertices is not None or (center is not None and radius is not None), \
        "--non_interactive requires --from_json, or --vertices, or --center and --radius"
if replay_source_puck is not None:
    source_puck = replay_source_puck
print(f"source_puck = {source_puck}")

# Load the puck - keep the raw lines so the output is written back verbatim
with open(puck) as f:
    lines = f.read().splitlines()
lines = [l for l in lines if l]
n = len(lines)
assert n > 0, f"{puck} is empty"
xy = np.fromiter((float(v) for l in lines for v in l.split(',', 1)[1].split(',')),
                 dtype=np.float64, count=2 * n).reshape(-1, 2)
x, y = xy[:, 0], xy[:, 1]
print(f"{n} beads: x [{x.min():.1f}, {x.max():.1f}], y [{y.min():.1f}, {y.max():.1f}]")

def circle_mask(center, radius):
    """Boolean mask of the beads inside the circle."""
    if center is None or radius is None or radius <= 0:
        return np.zeros(n, dtype=bool)
    dx = x - center[0]
    dy = y - center[1]
    return dx * dx + dy * dy <= radius * radius

def polygon_mask(verts):
    """Boolean mask of the beads inside the closed polygon.

    Crossing-number test, vectorized over the beads and looped over the edges, so the
    peak allocation is one bead-sized array rather than beads x edges. Written in plain
    numpy on purpose: matplotlib.path would do this too, but --non_interactive must keep
    working without a plotting stack (see the lazy import below).
    """
    verts = np.asarray(verts, dtype=np.float64)
    if verts.shape[0] < 3:
        return np.zeros(n, dtype=bool)
    # Bounding box first - on a 2M-bead puck this keeps the per-edge work off everything
    # the polygon could not possibly contain
    (xlo, ylo), (xhi, yhi) = verts.min(axis=0), verts.max(axis=0)
    cand = np.flatnonzero((x >= xlo) & (x <= xhi) & (y >= ylo) & (y <= yhi))
    px, py = x[cand], y[cand]
    inside = np.zeros(cand.size, dtype=bool)
    for (x1, y1), (x2, y2) in zip(verts, np.roll(verts, 1, axis=0)):
        # Slice to the straddling beads before dividing: a horizontal edge straddles
        # nothing, so the y2 - y1 == 0 division never happens
        i = np.flatnonzero((y1 > py) != (y2 > py))
        if i.size:
            inside[i] ^= px[i] < x1 + (py[i] - y1) * (x2 - x1) / (y2 - y1)
    m = np.zeros(n, dtype=bool)
    m[cand] = inside
    return m

def mask_of(shape):
    """Boolean mask of the beads inside (or outside, if inverted) the shape."""
    if shape["shape"] == "polygon":
        m = polygon_mask(shape["vertices"])
    else:
        m = circle_mask(shape["center"], shape["radius"])
    return ~m if shape["invert"] else m

def shape_json(shape):
    """The shape half of the sidecar. A circle emits exactly the keys it always has, so
    sidecars stay readable by Terra method snapshots that predate polygons."""
    if shape["shape"] == "polygon":
        return {"shape": "polygon",
                "vertices": [[float(vx), float(vy)] for vx, vy in shape["vertices"]]}
    return {"center": [float(shape["center"][0]), float(shape["center"][1])],
            "radius": float(shape["radius"])}

def reproduce_flags(shape):
    if shape["shape"] == "polygon":
        flags = "--vertices " + " ".join(f"{v:.6g}" for xy in shape["vertices"] for v in xy)
    else:
        flags = f"--center {shape['center'][0]:.6g} {shape['center'][1]:.6g} --radius {shape['radius']:.6g}"
    return flags + (" --invert" if shape["invert"] else "")

def write_selection(shape):
    m = mask_of(shape)
    with open(out, 'w') as f:
        f.writelines(l + '\n' for l in (lines[i] for i in np.flatnonzero(m)))
    print(f"wrote {m.sum()}/{n} beads to {out}")

    # The sidecar is the source of truth for the cut: terra/reconstruction.wdl replays it
    # through --non_interactive --from_json to rebuild this same selection in the cloud,
    # and terra/submit.py discovers submittable selections by listing these files.
    with open(sidecar, 'w') as f:
        json.dump({"source_puck": source_puck,
                   **shape_json(shape),
                   "invert": bool(shape["invert"]),
                   "n_selected": int(m.sum()),
                   "n_total": n}, f, indent=2)
        f.write('\n')
    print(f"wrote {sidecar}")
    if sidecar_uri is not None:
        subprocess.run(["gcloud", "storage", "cp", sidecar, sidecar_uri], check=True)
        print(f"uploaded {sidecar_uri}")

    print(f"  reproduce with: {reproduce_flags(shape)}")

if non_interactive:
    write_selection(shape)
    sys.exit(0)

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import PolygonSelector

# A pre-seeded cut opens in its own mode with its own geometry already in place
if shape["shape"] == "circle":
    center, radius = shape["center"], shape["radius"]
invert = shape["invert"]

# Default the circle to a disk covering the bulk of the beads
if center is None:
    center = [float(np.median(x)), float(np.median(y))]
if radius is None:
    radius = float(np.percentile(np.hypot(x - center[0], y - center[1]), 95))

# Free up the keys we bind below from matplotlib's default keymap
# (s = savefig, r/left/right = home/back/forward in the view history, p = pan)
for rc, drop in [('keymap.save', 's'), ('keymap.home', 'r'), ('keymap.back', 'left'),
                 ('keymap.forward', 'right'), ('keymap.pan', 'p')]:
    plt.rcParams[rc] = [k for k in plt.rcParams[rc] if k != drop]

# Set up the figure - shape it like the data so the equal-aspect axes fills the window
aspect = np.clip((x.max() - x.min()) / max(y.max() - y.min(), 1e-9), 0.5, 2.0)
fig, ax = plt.subplots(figsize=(9 * np.sqrt(aspect), 9 / np.sqrt(aspect)))
fig.canvas.manager.set_window_title(f"puck-select: {puck}")
ax.plot(x, y, marker, color='0.75', linestyle='none', rasterized=True)
highlight, = ax.plot([], [], marker, color='red', linestyle='none', rasterized=True)
circle = Circle(center, radius, fill=False, edgecolor='blue', linewidth=1.5, linestyle='--')
ax.add_patch(circle)
ax.set_aspect('equal')
ax.set_xlabel("xcoord")
ax.set_ylabel("ycoord")
fig.text(0.5, 0.015,
         "circle: drag = set circle (press at center)   scroll = radius (+shift = fine)\n"
         "polygon: click = drop vertex, click the first one again to close, then drag handles   esc = restart\n"
         "p = circle/polygon   arrows = nudge (+shift = fine)   i = invert   r = reset   enter/s = save   q = quit",
         ha='center', va='bottom', fontsize=9, color='0.35')
fig.subplots_adjust(bottom=0.16)

state = {'mode': shape["shape"], 'center': list(center), 'radius': radius,
         'vertices': [list(v) for v in shape.get("vertices", [])],
         'invert': invert, 'dragging': False, 'stale': True}

# Built on first use: a live PolygonSelector swallows the button events the circle drag
# needs, so it must not exist while we are in circle mode
selector = None

def polygon_selector():
    global selector
    if selector is None:
        selector = PolygonSelector(ax, on_polygon, useblit=True,
                                   props=dict(color='blue', linestyle='--', linewidth=1.5))
        if len(state['vertices']) >= 3:
            selector.verts = state['vertices']
    return selector

def on_polygon(verts):
    """PolygonSelector callback: fires once the ring closes, and on every later edit."""
    state['vertices'] = [[float(vx), float(vy)] for vx, vy in verts]
    state['stale'] = True
    refresh_selection()

def current_shape():
    if state['mode'] == 'polygon':
        return {"shape": "polygon", "vertices": state['vertices'], "invert": state['invert']}
    return {"shape": "circle", "center": state['center'], "radius": state['radius'],
            "invert": state['invert']}

def set_mode(mode):
    state['mode'] = mode
    circle.set_visible(mode == 'circle')
    sel = polygon_selector() if mode == 'polygon' else selector
    if sel is not None:
        sel.set_active(mode == 'polygon')
        sel.set_visible(mode == 'polygon')
    state['stale'] = True
    refresh_selection()
    fig.canvas.draw_idle()

def draw_circle():
    """Cheap redraw: move the outline only, mark the selection as needing a recount."""
    circle.set_center(state['center'])
    circle.set_radius(state['radius'])
    state['stale'] = True
    set_title(count=None)
    fig.canvas.draw_idle()

def set_title(count):
    head = f"selected: {count:,} / {n:,}" if count is not None else f"selected: ... / {n:,}"
    if state['mode'] == 'polygon':
        detail = f"polygon ({len(state['vertices'])} vertices)"
    else:
        cx, cy = state['center']
        detail = f"center ({cx:.1f}, {cy:.1f})  radius {state['radius']:.1f}"
    ax.set_title(f"{head}    {detail}" + ("  [INVERTED]" if state['invert'] else ""))

def refresh_selection():
    """Expensive redraw: recompute the mask over all beads and repaint the highlight."""
    if not state['stale']:
        return
    m = mask_of(current_shape())
    highlight.set_data(x[m], y[m])
    state['stale'] = False
    set_title(int(m.sum()))
    fig.canvas.draw_idle()

def toolbar_active():
    tb = fig.canvas.manager.toolbar
    return bool(tb and getattr(tb, 'mode', ''))

def on_press(event):
    if state['mode'] != 'circle' or event.inaxes is not ax or event.button != 1 or toolbar_active():
        return
    state['dragging'] = True
    state['center'] = [event.xdata, event.ydata]
    state['radius'] = 0.0
    draw_circle()

def on_motion(event):
    if not state['dragging'] or event.inaxes is not ax:
        return
    cx, cy = state['center']
    state['radius'] = float(np.hypot(event.xdata - cx, event.ydata - cy))
    draw_circle()

def on_release(event):
    if not state['dragging']:
        return
    state['dragging'] = False
    refresh_selection()

def on_scroll(event):
    if state['mode'] != 'circle':
        return
    step = 1.005 if (event.key == 'shift') else 1.05
    state['radius'] *= step if event.button == 'up' else 1 / step
    draw_circle()
    refresh_selection()

def nudge(dx, dy, frac):
    """Translate the live shape. The circle scales its step by the radius, the polygon by
    its bounding box, so a tap moves each by a comparable fraction of its own size."""
    if state['mode'] == 'polygon':
        if len(state['vertices']) < 3:
            return
        v = np.asarray(state['vertices'])
        step = frac * float(max(v.max(axis=0) - v.min(axis=0)))
        state['vertices'] = (v + [dx * step, dy * step]).tolist()
        polygon_selector().verts = state['vertices']
        state['stale'] = True
    else:
        state['center'][0] += dx * frac * state['radius']
        state['center'][1] += dy * frac * state['radius']
        draw_circle()
    refresh_selection()

def on_key(event):
    key = event.key or ''
    nudges = {'left': (-1, 0), 'right': (1, 0), 'up': (0, 1), 'down': (0, -1)}
    base = key.split('+')[-1]
    if base in nudges:
        nudge(*nudges[base], 0.001 if key.startswith('shift+') else 0.01)
    elif key == 'p':
        set_mode('circle' if state['mode'] == 'polygon' else 'polygon')
    elif key == 'i':
        state['invert'] = not state['invert']
        state['stale'] = True
        refresh_selection()
    elif key == 'r':
        state['invert'] = False
        if state['mode'] == 'polygon':
            state['vertices'] = []
            polygon_selector().clear()
            state['stale'] = True
        else:
            state['radius'] = 0.0
            draw_circle()
        refresh_selection()
    elif key in ('enter', 's'):
        refresh_selection()
        write_selection(current_shape())

for name, handler in [('button_press_event', on_press), ('motion_notify_event', on_motion),
                      ('button_release_event', on_release), ('scroll_event', on_scroll),
                      ('key_press_event', on_key)]:
    fig.canvas.mpl_connect(name, handler)

set_mode(state['mode'])
plt.show()
