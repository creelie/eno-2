import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'code'))
import construct

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
plt.rcParams.update({"font.family": "serif", "font.size": 11, "mathtext.fontset": "cm"})

def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# ---------- Figure 1: sparsest recurrent configuration on a 9 x 8 rectangle ----------
L, M = 9, 8
b = construct.run((L, M))
fig = plt.figure(figsize=(7.2, 5.2))
ax = fig.add_subplot(111, projection="3d")
for x in range(L):
    for y in range(M):
        ax.bar3d(x, y, 0, 0.92, 0.92, 0.04, color="#d9d9d9", edgecolor="#9a9a9a", linewidth=0.3, shade=False)
        if b.z[(x, y)]:
            ax.bar3d(x + 0.08, y + 0.08, 0.04, 0.76, 0.76, 0.9, color="#c0603a", edgecolor="#5a2a18", linewidth=0.4, shade=True)
ax.set_box_aspect((L, M, 2.2))
ax.set_xticks(np.arange(L) + 0.46); ax.set_xticklabels([str(i + 1) for i in range(L)])
ax.set_yticks(np.arange(M) + 0.46); ax.set_yticklabels([str(i + 1) for i in range(M)])
ax.set_zticks([])
ax.set_xlabel("$x$"); ax.set_ylabel("$y$")
ax.view_init(elev=38, azim=-58)
ax.grid(False)
for p in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
    p.set_visible(False)
ax.set_title(f"{b.grains()} grains on the $9\\times 8$ rectangle, $\\lceil 8\\cdot 7/2\\rceil = 28$", pad=2)
save(fig, "fig_minimal.png")

# ---------- Figure 2: the last-toppling corner of a unit square ----------
fig = plt.figure(figsize=(6.4, 5.4))
ax = fig.add_subplot(111, projection="3d")
corners = {(0, 0): 3.0, (1, 0): 5.0, (1, 1): 8.0, (0, 1): 1.6}
names = {(0, 0): "v", (1, 0): "u", (1, 1): "w", (0, 1): "v'"}
sq = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
ax.add_collection3d(Poly3DCollection([sq], facecolor=(0.85,0.9,0.96,0.35), edgecolor="#335", linewidth=1.2, zorder=0))
for (x, y), t in corners.items():
    col = "#c0603a" if (x, y) == (1, 1) else "#335c8a"
    ax.plot([x, x], [y, y], [0, t], color=col, linewidth=2.2, zorder=5)
    ax.scatter([x], [y], [t], color=col, s=40, depthshade=False)
    ax.text(x + 0.06, y, t + 0.3, "$t_{%s}$" % names[(x, y)], fontsize=12, color=col)
    dx = -0.16 if x == 0 else 0.07
    dy = -0.16 if y == 0 else 0.07
    ax.text(x + dx, y + dy, 0.0, "$%s$" % names[(x, y)], fontsize=14, color="black")
ax.plot([1, 1], [1, 0], [0.01, 0.01], color="#c0603a", linewidth=4.0)
ax.quiver(1, 0.75, 0.6, 0, -0.6, -0.55, color="#c0603a", linewidth=2.2, arrow_length_ratio=0.25)
ax.text(1.18, 0.55, 1.6, "grain to $u$", fontsize=11, color="#c0603a")
ax.set_xlim(-0.3, 1.3); ax.set_ylim(-0.3, 1.3); ax.set_zlim(0, 9)
ax.set_zlabel("time of last toppling", labelpad=6)
ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([0, 2, 4, 6, 8])
ax.view_init(elev=22, azim=-62)
for p in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
    p.set_visible(False)
ax.grid(False)
ax.set_box_aspect((1, 1, 1.4))
save(fig, "fig_square.png")

# ---------- Figure 3: the layer step in three dimensions ----------
def cuboid(ax, x0, x1, y0, y1, z0, z1, fc, ec, alpha):
    X = [x0, x1]; Y = [y0, y1]; Z = [z0, z1]
    faces = []
    for (a, b_, c) in [(0, 1, 2)]:
        pass
    v = np.array([[x, y, z] for x in X for y in Y for z in Z])
    idx = [[0, 1, 3, 2], [4, 5, 7, 6], [0, 1, 5, 4], [2, 3, 7, 6], [0, 2, 6, 4], [1, 3, 7, 5]]
    ax.add_collection3d(Poly3DCollection([v[i] for i in idx], facecolor=fc, edgecolor=ec, linewidth=0.8, alpha=alpha))
fig = plt.figure(figsize=(7.0, 4.6))
ax = fig.add_subplot(111, projection="3d")
cuboid(ax, 0, 5, 0, 4, 0, 3, "#9fb8d6", "#2c3e50", 0.35)
cuboid(ax, 5, 6, 0, 4, 0, 3, "#f2c14e", "#7a5a00", 0.55)
cuboid(ax, 6, 7, 0, 4, 0, 3, "#e9e9e9", "#555", 0.35)
fig.text(0.10, 0.80, "layers $1,\\dots,\\ell_1$:\nsparsest by induction", fontsize=10, color="#2c3e50")
fig.text(0.40, 0.12, "layer $\\ell_1+1$: brought to a sparsest state\nof the $(d-1)$-dimensional model", fontsize=10, color="#7a5a00")
fig.text(0.70, 0.80, "layer $\\ell_1+2$:\nemptied along axis 1", fontsize=10, color="#333333")
for y in (1, 3):
    for z in (1, 2):
        ax.quiver(6.5, y, z, -0.9, 0, 0, color="#c0603a", linewidth=1.4, arrow_length_ratio=0.3)
        ax.quiver(6.5, y, z, 0.9, 0, 0, color="#c0603a", linewidth=1.4, arrow_length_ratio=0.3)
ax.set_xlim(0, 8); ax.set_ylim(-1, 5); ax.set_zlim(-1, 4)
ax.set_box_aspect((8, 6, 5))
ax.set_axis_off()
ax.view_init(elev=24, azim=-60)
save(fig, "fig_layers.png")

# ---------- Figure 4: sparsest recurrent configuration in a 5x5x5 box ----------
D = (5, 5, 5)
bb = construct.run(D)
fig = plt.figure(figsize=(6.0, 5.6))
ax = fig.add_subplot(111, projection="3d")
filled = np.zeros(D, dtype=bool)
for s, v in bb.z.items():
    if v:
        filled[s] = True
ax.voxels(filled, facecolors="#c0603a", edgecolors="#5a2a18", linewidth=0.4)
grid = np.ones(D, dtype=bool)
ax.voxels(grid, facecolors=(0.6, 0.6, 0.7, 0.04), edgecolors=(0.3, 0.3, 0.4, 0.12), linewidth=0.3)
ax.set_box_aspect((1, 1, 1))
ax.set_axis_off()
ax.view_init(elev=24, azim=-52)
ax.set_title(f"{bb.grains()} grains in the $5\\times5\\times5$ box, $4^3/4 = 16$", pad=0)
save(fig, "fig_box3d.png")
print("grains 2D", b.grains(), "3D", bb.grains())
