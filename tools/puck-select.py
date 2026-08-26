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
    parser = argparse.ArgumentParser(description='Interactively select a disk of beads from a Puck.csv')
    parser.add_argument("puck", type=str, help="path to Puck.csv (headerless sb,x,y), local or gs://")
    parser.add_argument("-s", "--selection", type=str, default=None, help="name this selection; required for a gs:// puck, whose sidecar is uploaded to <recon_dir>/<selection>/selection.json")
    parser.add_argument("-o", "--out", type=str, default=None, help="output csv (default: Puck-selected.csv beside the input)")
    parser.add_argument("-c", "--center", type=float, nargs=2, default=None, metavar=('X', 'Y'), help="pre-seed the circle center")
    parser.add_argument("-r", "--radius", type=float, default=None, help="pre-seed the circle radius")
    parser.add_argument("-v", "--invert", action='store_true', help="select points outside the circle")
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
invert = args.invert                   ; print(f"invert = {invert}")
marker = args.marker                   ; print(f"marker = {marker}")
non_interactive = args.non_interactive  ; print(f"non_interactive = {non_interactive}")

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
    assert center is not None and radius is not None, "--non_interactive requires --center and --radius"
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

def mask_of(center, radius, invert):
    """Boolean mask of the beads inside (or outside, if inverted) the circle."""
    if center is None or radius is None or radius <= 0:
        return np.zeros(n, dtype=bool)
    dx = x - center[0]
    dy = y - center[1]
    m = dx * dx + dy * dy <= radius * radius
    return ~m if invert else m

def write_selection(center, radius, invert):
    m = mask_of(center, radius, invert)
    with open(out, 'w') as f:
        f.writelines(l + '\n' for l in (lines[i] for i in np.flatnonzero(m)))
    print(f"wrote {m.sum()}/{n} beads to {out}")

    # The sidecar is the source of truth for the cut: terra/reconstruction.wdl replays it
    # through --non_interactive to rebuild this same selection in the cloud, and
    # terra/submit.py discovers submittable selections by listing these files.
    with open(sidecar, 'w') as f:
        json.dump({"source_puck": source_puck,
                   "center": [float(center[0]), float(center[1])],
                   "radius": float(radius),
                   "invert": bool(invert),
                   "n_selected": int(m.sum()),
                   "n_total": n}, f, indent=2)
        f.write('\n')
    print(f"wrote {sidecar}")
    if sidecar_uri is not None:
        subprocess.run(["gcloud", "storage", "cp", sidecar, sidecar_uri], check=True)
        print(f"uploaded {sidecar_uri}")

    print(f"  reproduce with: --center {center[0]:.6g} {center[1]:.6g} --radius {radius:.6g}" + (" --invert" if invert else ""))

if non_interactive:
    write_selection(center, radius, invert)
    sys.exit(0)

import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Default the circle to a disk covering the bulk of the beads
if center is None:
    center = [float(np.median(x)), float(np.median(y))]
if radius is None:
    radius = float(np.percentile(np.hypot(x - center[0], y - center[1]), 95))

# Free up the keys we bind below from matplotlib's default keymap
# (s = savefig, r/left/right = home/back/forward in the view history)
for rc, drop in [('keymap.save', 's'), ('keymap.home', 'r'), ('keymap.back', 'left'), ('keymap.forward', 'right')]:
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
         "drag = set circle (press at center)   scroll = radius (+shift = fine)   arrows = nudge center (+shift = fine)\n"
         "i = invert   r = reset   enter/s = save   q = quit",
         ha='center', va='bottom', fontsize=9, color='0.35')
fig.subplots_adjust(bottom=0.13)

state = {'center': list(center), 'radius': radius, 'invert': invert, 'dragging': False, 'stale': True}

def draw_circle():
    """Cheap redraw: move the outline only, mark the selection as needing a recount."""
    circle.set_center(state['center'])
    circle.set_radius(state['radius'])
    state['stale'] = True
    set_title(count=None)
    fig.canvas.draw_idle()

def set_title(count):
    cx, cy = state['center']
    head = f"selected: {count:,} / {n:,}" if count is not None else f"selected: ... / {n:,}"
    ax.set_title(f"{head}    center ({cx:.1f}, {cy:.1f})  radius {state['radius']:.1f}"
                 + ("  [INVERTED]" if state['invert'] else ""))

def refresh_selection():
    """Expensive redraw: recompute the mask over all beads and repaint the highlight."""
    if not state['stale']:
        return
    m = mask_of(state['center'], state['radius'], state['invert'])
    highlight.set_data(x[m], y[m])
    state['stale'] = False
    set_title(int(m.sum()))
    fig.canvas.draw_idle()

def toolbar_active():
    tb = fig.canvas.manager.toolbar
    return bool(tb and getattr(tb, 'mode', ''))

def on_press(event):
    if event.inaxes is not ax or event.button != 1 or toolbar_active():
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
    step = 1.005 if (event.key == 'shift') else 1.05
    state['radius'] *= step if event.button == 'up' else 1 / step
    draw_circle()
    refresh_selection()

def on_key(event):
    key = event.key or ''
    nudges = {'left': (-1, 0), 'right': (1, 0), 'up': (0, 1), 'down': (0, -1)}
    base = key.split('+')[-1]
    if base in nudges:
        frac = 0.001 if key.startswith('shift+') else 0.01
        dx, dy = nudges[base]
        state['center'][0] += dx * frac * state['radius']
        state['center'][1] += dy * frac * state['radius']
        draw_circle()
        refresh_selection()
    elif key == 'i':
        state['invert'] = not state['invert']
        state['stale'] = True
        refresh_selection()
    elif key == 'r':
        state['radius'] = 0.0
        state['invert'] = False
        draw_circle()
        refresh_selection()
    elif key in ('enter', 's'):
        refresh_selection()
        write_selection(state['center'], state['radius'], state['invert'])

for name, handler in [('button_press_event', on_press), ('motion_notify_event', on_motion),
                      ('button_release_event', on_release), ('scroll_event', on_scroll),
                      ('key_press_event', on_key)]:
    fig.canvas.mpl_connect(name, handler)

refresh_selection()
plt.show()
