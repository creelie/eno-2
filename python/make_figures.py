"""Black-and-white figures for the paper. Run: python3 make_figures.py OUTDIR"""
import sys, itertools, collections
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from construct import build

out = sys.argv[1] if len(sys.argv) > 1 else "."
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12})
DARK, MID, LIGHT, PALE = "0.15", "0.45", "0.75", "0.93"

def cube_faces(o, s=(1, 1, 1)):
    x, y, z = o; a, b, c = s
    v = np.array([[x, y, z], [x + a, y, z], [x + a, y + b, z], [x, y + b, z],
                  [x, y, z + c], [x + a, y, z + c], [x + a, y + b, z + c], [x, y + b, z + c]])
    f = [[0, 1, 2, 3], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [1, 2, 6, 5], [0, 3, 7, 4]]
    return [v[i] for i in f]

def add_cube(ax, o, s=(1, 1, 1), fc=MID, ec="k", alpha=1.0, lw=0.6):
    pc = Poly3DCollection(cube_faces(o, s), facecolors=fc, edgecolors=ec, linewidths=lw, alpha=alpha)
    ax.add_collection3d(pc)

def box_wire(ax, dims, off=(0, 0, 0), color="0.3", lw=0.9):
    a, b, c = dims; x0, y0, z0 = off
    for e in itertools.product([0, 1], repeat=3):
        for k in range(3):
            if e[k] == 0:
                p = np.array(e, float); q = p.copy(); q[k] = 1
                P = np.array([x0, y0, z0]) + p * [a, b, c]; Q = np.array([x0, y0, z0]) + q * [a, b, c]
                ax.plot(*zip(P, Q), color=color, lw=lw)

def clean3d(ax):
    ax.set_axis_off()
    ax.set_box_aspect(None)

# 1. same-parity family in [4]x[4]x[5]
def fig_family():
    dims = (4, 4, 5)
    fig = plt.figure(figsize=(6, 5.2)); ax = fig.add_subplot(projection="3d")
    for c in itertools.product(*[range(1, k) for k in dims]):
        if all(t % 2 for t in c):
            add_cube(ax, c, fc=DARK, ec="k", alpha=0.85)
        elif all(t % 2 == 0 for t in c):
            add_cube(ax, c, fc=LIGHT, ec="k", alpha=0.85)
    pts = np.array(list(itertools.product(*[range(1, k + 1) for k in dims])))
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=6, c="k", depthshade=False)
    box_wire(ax, (3, 3, 4), off=(1, 1, 1), color="0.4", lw=0.8)
    ax.set_xlim(1, 4); ax.set_ylim(1, 4); ax.set_zlim(1, 5); ax.set_box_aspect((3, 3, 4))
    ax.view_init(elev=22, azim=-58); ax.set_axis_off()
    fig.savefig(f"{out}/fig_family.png", dpi=220, bbox_inches="tight", pad_inches=0.02); plt.close(fig)

# 2./3. sparsest configurations as unit cubes in a lattice box
def grains3d(dims, name, elev=24, azim=-55):
    g = build(list(dims))
    fig = plt.figure(figsize=(5, 5)); ax = fig.add_subplot(projection="3d")
    box_wire(ax, dims, off=(0.5, 0.5, 0.5), color="0.35")
    for k in range(3):
        for t in range(1, dims[k]):
            # faint grid lines on the three back faces
            pass
    pts = np.array(list(itertools.product(*[range(1, k + 1) for k in dims])))
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], s=4, c="0.55", depthshade=False)
    for s in sorted(g, key=lambda s: (s[0] + s[1], s[2])):
        add_cube(ax, tuple(t - 0.4 for t in s), (0.8, 0.8, 0.8), fc=MID, ec="k", lw=0.7)
    ax.set_xlim(0.5, dims[0] + 0.5); ax.set_ylim(0.5, dims[1] + 0.5); ax.set_zlim(0.5, dims[2] + 0.5)
    ax.set_box_aspect(dims); ax.view_init(elev=elev, azim=azim); ax.set_axis_off()
    fig.savefig(f"{out}/{name}", dpi=220, bbox_inches="tight", pad_inches=0.02); plt.close(fig)

# 4. sparsest configuration on 9x8 as bars
def fig_minimal():
    g = set(build([9, 8]))
    fig = plt.figure(figsize=(6, 5)); ax = fig.add_subplot(projection="3d")
    for x in range(1, 10):
        for y in range(1, 9):
            add_cube(ax, (x - 0.45, y - 0.45, 0), (0.9, 0.9, 0.08), fc=PALE, ec="0.6", lw=0.4)
    for (x, y) in sorted(g, key=lambda s: (-s[1], s[0])):
        add_cube(ax, (x - 0.35, y - 0.35, 0.08), (0.7, 0.7, 1.0), fc=MID, ec="k", lw=0.6)
    ax.set_xlim(0.5, 9.5); ax.set_ylim(0.5, 8.5); ax.set_zlim(0, 2); ax.set_box_aspect((9, 8, 2), zoom=1.1)
    ax.set_xticks(range(1, 10)); ax.set_yticks(range(1, 9)); ax.set_zticks([])
    ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor((1, 1, 1, 0)); a.pane.set_edgecolor((1, 1, 1, 0))
    ax.grid(False); ax.zaxis.line.set_color((1, 1, 1, 0))
    ax.view_init(elev=32, azim=-62)
    fig.savefig(f"{out}/fig_minimal.png", dpi=220, bbox_inches="tight", pad_inches=0.05); plt.close(fig)

# 5. last-toppling times on a unit cell
def fig_square():
    fig = plt.figure(figsize=(5.2, 5.2)); ax = fig.add_subplot(projection="3d", computed_zorder=False)
    P = {"v": (0, 0), "v'": (0, 1), "u": (1, 0), "w": (1, 1)}
    T = {"v": 2.2, "v'": 2.8, "u": 3.6, "w": 7.8}
    lab = {"v": "t_v", "v'": "t_{v'}", "u": "t_u", "w": "t_w"}
    ax.add_collection3d(Poly3DCollection([[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]],
                        facecolors=PALE, edgecolors="k", linewidths=1.2, zorder=1))
    ax.plot([1, 1], [0.95, 0.12], [0.02, 0.02], color="k", lw=2.2, zorder=4)
    ax.quiver(1, 0.3, 0.02, 0, -0.25, 0, color="k", lw=2.2, arrow_length_ratio=0.9, zorder=4)
    for k, (x, y) in P.items():
        col = "k" if k == "w" else MID
        ax.plot([x, x], [y, y], [0, T[k]], color=col, lw=2.8 if k == "w" else 2.0, zorder=5)
        ax.scatter([x], [y], [T[k]], s=42, c=col, depthshade=False, zorder=6)
        dx = -0.22 if k == "w" else 0.05
        ax.text(x + dx, y + 0.05, T[k] + 0.3, f"${lab[k]}$", color=col, zorder=7)
        if k == "w":
            ax.text(x - 0.12, y + 0.1, -0.9, "$w$", zorder=7)
        else:
            ax.text(x + (0.08 if x else -0.16), y + (0.1 if y else -0.16), -0.5, f"${k}$", zorder=7)
    ax.text(1.1, 0.45, -1.4, "grain to $u$", zorder=7)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_zlim(0, 8.5); ax.set_box_aspect((1, 1, 1.25))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([0, 2, 4, 6, 8])
    ax.set_zlabel("time of last toppling", labelpad=10); ax.tick_params(axis="z", pad=6)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor((1, 1, 1, 0)); a.pane.set_edgecolor((1, 1, 1, 0))
    ax.grid(False); ax.xaxis.line.set_color((1, 1, 1, 0)); ax.yaxis.line.set_color((1, 1, 1, 0))
    ax.view_init(elev=20, azim=-52)
    fig.savefig(f"{out}/fig_square.png", dpi=220, bbox_inches="tight", pad_inches=0.08); plt.close(fig)

# 6. induction step of the layer construction
def fig_layers():
    fig = plt.figure(figsize=(7.5, 4.6)); ax = fig.add_axes([0.0, 0.08, 1.0, 0.8], projection="3d")
    add_cube(ax, (0, 0, 0), (5, 3, 2.4), fc=LIGHT, ec="k", alpha=0.55, lw=0.9)
    add_cube(ax, (5, 0, 0), (1, 3, 2.4), fc=DARK, ec="k", alpha=0.75, lw=0.9)
    add_cube(ax, (6, 0, 0), (1, 3, 2.4), fc="white", ec="k", alpha=0.25, lw=0.9)
    for y, z in ((0.8, 0.7), (1.5, 1.4), (2.2, 2.0)):
        ax.quiver(6.15, y, z, 1.2, 0, 0, color="k", lw=1.2, arrow_length_ratio=0.3)
    ax.set_xlim(0, 7.5); ax.set_ylim(0, 3); ax.set_zlim(0, 2.4); ax.set_box_aspect((7.5, 3, 2.4), zoom=0.95)
    ax.view_init(elev=20, azim=-62); ax.set_axis_off()
    fig.text(0.03, 0.92, "layers $1,\\dots,\\ell_1$:\nsparsest by induction", fontsize=12, va="top")
    fig.text(0.97, 0.92, "layer $\\ell_1+2$:\nemptied along axis 1", fontsize=12, va="top", ha="right")
    fig.text(0.5, 0.06, "layer $\\ell_1+1$ (dark): reset to full, then made sparsest\nas a box of dimension $d-1$", fontsize=12, ha="center", va="top")
    fig.savefig(f"{out}/fig_layers.png", dpi=220, bbox_inches="tight", pad_inches=0.05); plt.close(fig)

# 7. Oslo relaxation times at L = 5
def fig_relaxation(L=5):
    def inR(s):
        for x in range(L):
            if s[x] == 0:
                j = x + 1
                while j < L and s[j] == 1: j += 1
                if j == L or s[j] != 2: return False
        return True
    cnt = collections.Counter()
    for s in itertools.product(range(3), repeat=L):
        g = sum((i + 1) * s[i] for i in range(L)); K = L * (L + 1) - g + (s[-1] == 2)
        cnt[(g, K, s[-1] == 2, inR(s))] += 1
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    xs = np.array([0, L * (L + 1)])
    ax.plot(xs, L * (L + 1) - xs, color="0.5", lw=1, label="$K=L(L+1)-g$")
    ax.plot(xs, L * (L + 1) - xs + 1, color="0.5", lw=1, ls=":", label="$K=L(L+1)-g+1$")
    ax.axhline(L * (L + 1) / 2, color="k", ls="--", lw=1, label="$L(L+1)/2=15$")
    for (g, K, top, rec), n in sorted(cnt.items(), key=lambda kv: kv[0][3]):
        m = "s" if top else "o"; size = (18 + 9 * n ** 0.8) * (0.45 if rec else 1.0)
        ax.scatter([g], [K], marker=m, s=size, facecolors=("k" if rec else "white"),
                   edgecolors=("0.45" if top else "k"), linewidths=1.2, zorder=3)
    for lab, m, fc, ec in (("$\\sigma_5\\neq2$, recurrent", "o", "k", "k"), ("$\\sigma_5\\neq2$, transient", "o", "white", "k"),
                           ("$\\sigma_5=2$, recurrent", "s", "k", "0.45"), ("$\\sigma_5=2$, transient", "s", "white", "0.45")):
        ax.scatter([], [], marker=m, facecolors=fc, edgecolors=ec, s=40, label=lab)
    ax.set_xlabel(r"number of grains $g(\sigma)$"); ax.set_ylabel(r"relaxation time $K(\sigma)$")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=9.5, loc="upper right")
    fig.savefig(f"{out}/fig_relaxation.png", dpi=220, bbox_inches="tight", pad_inches=0.05); plt.close(fig)

fig_family(); grains3d((4, 4, 4), "fig_444.png"); grains3d((5, 5, 5), "fig_box3d.png")
fig_minimal(); fig_square(); fig_layers(); fig_relaxation()
print("figures written to", out)
