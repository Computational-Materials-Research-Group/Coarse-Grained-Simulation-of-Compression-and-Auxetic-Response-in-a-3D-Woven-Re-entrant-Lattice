"""
build_reentrant_woven.py - 3D woven RE-ENTRANT (auxetic) lattice as a beam
network for LAMMPS (BPM package, bond_style bpm/rotational).

Same graph-based weaving recipe as the woven simple-cubic lattice, but the
parent beam lattice is a re-entrant (bow-tie) honeycomb:

  1. parent lattice : re-entrant (hourglass) honeycomb in the x-z plane, Nx x Nz cells,
                      horizontal ribs h, near-vertical struts l (cell height = 2 l cos(theta)),
                      copied into Ny+1 layers along y and joined by y-beams (length Ly)
  2. topology graph : nodes = rib ends, edges = ribs, struts, y-beams
  3. node graph     : in-plane, every angular sector between neighbouring beams
                      carries one fibre turn (2-D analogue of the octahedron edges
                      of the cubic node); out of plane, each in-plane beam's +y / -y
                      fibre turns into the y-beam
  4. beam bundles   : in-plane beams carry 4 helical fibres, y-beams carry one fibre
                      per in-plane beam at that node; radius R_eff, twist n_rev
  5. woven lattice  : fibres run beam -> node turn -> beam ... continuously;
                      at the lattice boundary a fibre with no partner beam ends

Default geometry follows the drawing: ribs 17 mm, struts 14.7 mm, cell 17 x 28.7 mm,
beam thickness 1 mm (lengths in mm). Compression is along z, perpendicular to the ribs.

Outputs: woven_lattice.data (atom_style bpm/sphere), woven_lattice.bonds
(bond coefficients + v_dfib), woven_lattice_ovito.data (atom_style bond, for OVITO).
Units: lj (non-dimensional).
"""
import numpy as np
from scipy.spatial import cKDTree

# ======================= INPUTS ===========================================
Nx, Nz     = 4, 4          # hourglass cells in x (width) and z (height = loading direction)
Ny         = 4             # cells through the thickness (Ny+1 layers joined by y-beams)
h          = 17.0          # horizontal rib length (cell width)          [mm]
l          = 14.7          # inclined strut length                       [mm]
cell_h     = 28.7          # cell height = 2 l cos(theta)                [mm]
Ly         = 17.0          # y-beam length (layer spacing)               [mm]
t_beam     = 3.0           # bundle thickness 2 R_eff + d [mm]; 3-5 = weave clearly visible, 1 = thin solid-like beams
n_quarter  = 2             # in-plane beams: n_rev = n_quarter / 4 revolutions (0 = straight fibres)
d_ratio    = 1.2           # fibre diameter / R_eff (auto-reduced if fibres overlap)
wrap       = 0.7           # radius of fibre turns round wide (> 150 deg) node sectors, x node size
ds_ratio   = 0.8           # bead spacing / fibre diameter
fuse_nodes = True          # True: fibres bonded where they touch inside the nodes (needed for auxetic response)
fuse_gap   = 0.5           #   fuse bead pairs of different fibres closer than (1 + fuse_gap) x d
# beam material (per fibre)
E, nu, rho = 1.0, 0.30, 1.0
sigma_f    = 0.25          # failure stress/E -> ~25 % fibre failure strain (0.05 = brittle: waists snap at ~2 %)
zeta       = 0.3           # fraction of critical damping
# ==========================================================================

assert cell_h < 2 * l, "cell_h must be < 2 l"
th = -np.arccos(cell_h / 2 / l)                    # re-entrant angle (negative)
theta_deg = np.degrees(th)
lc, ls = l * np.cos(th), -l * np.sin(th)          # strut extent across / back along the rib (> 0)
H = 2 * h - 2 * ls                                 # period along the ribs
assert H > 0 and h > 2 * ls, "struts too long / angle too steep: ribs would cross"
R = t_beam / (2 + d_ratio)
R_ratio = R / l
n_rev = n_quarter / 4.0
d_target = d_ratio * R
ds = ds_ratio * d_target
Y = np.array([0., 1., 0.])

# ---- parent lattice (one layer), then copied in y --------------------------
# rib (i, k): horizontal, at z = i*lc, from x = k*H + (i%2)*H/2 to x + h
# strut: right end of rib (i,k) -> left end of rib (i+-1, k + i%2)
# full hourglass columns k = 0..Nx-1; the waist ribs (odd i) join neighbouring columns,
# so the last column has none on its right
ribs = [(i, k) for i in range(2 * Nz + 1) for k in range(Nx) if not (i % 2 and k == Nx - 1)]
def xleft(i, k): return k * H + (i % 2) * H / 2
nodes2d = {}                                        # key -> (x, z)
for (i, k) in ribs:
    nodes2d[(i, k, 0)] = (xleft(i, k), i * lc)      # 0 = rib left end
    nodes2d[(i, k, 1)] = (xleft(i, k) + h, i * lc)  # 1 = rib right end
# waist vertices of the outer hourglass sides (left and right edges), without rib stubs
for i in range(1, 2 * Nz, 2):
    nodes2d[(i, -1, 1)] = (xleft(i, -1) + h, i * lc)
    nodes2d[(i, Nx - 1, 0)] = (xleft(i, Nx - 1), i * lc)
edges2d = [((i, k, 0), (i, k, 1)) for (i, k) in ribs]
for (i, k, e) in list(nodes2d):
    if e == 1:
        for j in (i - 1, i + 1):
            if (j, k + i % 2, 0) in nodes2d:
                edges2d.append(((i, k, 1), (j, k + i % 2, 0)))
used2d = {n for ed in edges2d for n in ed}
nodes2d = {kk: v for kk, v in nodes2d.items() if kk in used2d}

nodes, nid = [], {}
for jy in range(Ny + 1):
    for key2, (x, z) in nodes2d.items():
        nid[(jy,) + key2] = len(nodes); nodes.append(np.array([x, jy * Ly, z]))
nodes = np.array(nodes)
beams = []                                          # (node A, node B, kind)
for jy in range(Ny + 1):
    for a, b in edges2d:
        beams.append((nid[(jy,) + a], nid[(jy,) + b], "p"))
for jy in range(Ny):
    for key2 in nodes2d:
        beams.append((nid[(jy,) + key2], nid[(jy + 1,) + key2], "y"))

nb = len(beams)
bdir = np.array([(nodes[b] - nodes[a]) / np.linalg.norm(nodes[b] - nodes[a]) for a, b, _ in beams])
blen = np.array([np.linalg.norm(nodes[b] - nodes[a]) for a, b, _ in beams])
ends = [[] for _ in nodes]                          # node -> [(beam, end, outward dir)]
for ib, (a, b, _) in enumerate(beams):
    ends[a].append((ib, 0, bdir[ib])); ends[b].append((ib, 1, -bdir[ib]))

def inplane_out(n):                                 # in-plane beams at node, ccw order in x-z
    lst = [(ib, e, d) for ib, e, d in ends[n] if beams[ib][2] == "p"]
    return sorted(lst, key=lambda t: np.arctan2(t[2][2], t[2][0]))
def dkey(d): return (round(d[0], 4), round(d[2], 4))

# ---- beam frames and fibre slots ----------------------------------------------
frames, slot_ang, slot_lab = [], [], []             # slot_lab[ib] = label -> slot index
for ib, (a, b, kind) in enumerate(beams):
    d = bdir[ib]
    if kind == "p":
        e1, e2 = Y, np.cross(d, Y)                  # e2 = ccw in-plane normal of d
        ang = np.array([0, .5, 1, 1.5]) * np.pi     # slots +y, +n, -y, -n
        lab = {}
    else:
        e1, e2 = np.array([0., 0, 1]), np.array([1., 0, 0])
        dirs = [t[2] for t in inplane_out(a)]       # in-plane beams at this (x,z) site
        act = np.array([np.arctan2(v @ e2, v @ e1) for v in dirs])
        order = np.argsort(act); K = len(dirs)
        act = act[order]
        off = np.angle(np.mean(np.exp(1j * (act - 2 * np.pi * np.arange(K) / K))))
        ang = off + 2 * np.pi * np.arange(K) / K    # evenly spread, same cyclic order
        lab = {dkey(dirs[o]): s for s, o in enumerate(order)}
    frames.append((e1, e2)); slot_ang.append(ang); slot_lab.append(lab)

def svec(ib, s):
    e1, e2 = frames[ib]; a = slot_ang[ib][s]
    return np.cos(a) * e1 + np.sin(a) * e2

# node region size: fibre bundles of neighbouring beams must not overlap
delta = np.zeros(len(nodes))
for n in range(len(nodes)):
    ds_ = [t[2] for t in ends[n]]
    amin = np.pi / 2
    for p in range(len(ds_)):
        for q in range(p + 1, len(ds_)):
            amin = min(amin, np.arccos(np.clip(ds_[p] @ ds_[q], -1, 1)))
    delta[n] = R / np.tan(amin / 2) + d_target

# ---- helices -------------------------------------------------------------------
slot_pt, slot_tan, curves = {}, {}, []
for ib, (a, b, kind) in enumerate(beams):
    d, (e1, e2), ang = bdir[ib], frames[ib], slot_ang[ib]
    K = len(ang)
    dA, dB = delta[a], delta[b]
    Lh = blen[ib] - dA - dB
    assert Lh > 2 * ds, f"beam {ib} too short for R_eff: increase h/l/Ly or reduce t_beam"
    nr = n_rev if kind == "p" else round(n_rev * K) / K   # closest twist the bundle allows
    best = None
    for m in range(K):
        for N in range(4):
            D = np.array([(ang[(k + m) % K] - ang[k]) % (2 * np.pi) + 2 * np.pi * N for k in range(K)])
            err = abs(D.mean() - 2 * np.pi * nr)
            if best is None or err < best[0] - 1e-9:
                best = (err, m, D)
    _, m, D = best
    s = np.linspace(dA, blen[ib] - dB, 400)
    for k in range(K):
        t = ang[k] + D[k] * (s - dA) / Lh
        P = nodes[a][None, :] + s[:, None] * d + R * (np.cos(t)[:, None] * e1 + np.sin(t)[:, None] * e2)
        T = d[None, :] + R * (D[k] / Lh) * (-np.sin(t)[:, None] * e1 + np.cos(t)[:, None] * e2)
        sA, sB = (ib, 0, k), (ib, 1, (k + m) % K)
        slot_pt[sA], slot_tan[sA] = P[0], T[0] / np.linalg.norm(T[0])
        slot_pt[sB], slot_tan[sB] = P[-1], -T[-1] / np.linalg.norm(T[-1])
        curves.append((sA, sB, P))

# ---- node turns ---------------------------------------------------------------------
def nearest_slot(ib, e, v):
    return (ib, e, int(np.argmax([svec(ib, s) @ v for s in range(len(slot_ang[ib]))])))
turns = []
for n in range(len(nodes)):
    ip = inplane_out(n)
    if len(ip) >= 2:                                # one fibre turn per angular sector
        for q in range(len(ip)):
            (ib, e, di), (jb, f, dj) = ip[q], ip[(q + 1) % len(ip)]
            ccw = np.array([-di[2], 0, di[0]]); cw = np.array([dj[2], 0, -dj[0]])
            ai, aj = np.arctan2(di[2], di[0]), np.arctan2(dj[2], dj[0])
            phi = (aj - ai) % (2 * np.pi) if len(ip) > 1 else 2 * np.pi
            way = None
            if phi > np.radians(150):               # wide sector: wrap round the node
                b = ai + phi / 2
                bv = np.array([np.cos(b), 0, np.sin(b)])
                way = (nodes[n] + wrap * np.hypot(delta[n], R) * bv, np.array([-bv[2], 0, bv[0]]))
            turns.append((nearest_slot(ib, e, ccw), nearest_slot(jb, f, cw), way))
    yup = [(ib, e) for ib, e, d in ends[n] if beams[ib][2] == "y" and e == 0]
    ydn = [(ib, e) for ib, e, d in ends[n] if beams[ib][2] == "y" and e == 1]
    for ib, e, di in ip:                            # +y / -y fibres turn into the y-beams
        for (yb, ye), sgn in [(t, 1) for t in yup] + [(t, -1) for t in ydn]:
            turns.append((nearest_slot(ib, e, sgn * Y), (yb, ye, slot_lab[yb][dkey(di)]), None))
    # front/back face: no y-beam on that side, so the face fibres of the in-plane
    # beams turn into each other instead (widest-angle pairs first; one may stay free)
    for sgn, ylist in [(1, yup), (-1, ydn)]:
        if ylist or len(ip) < 2:
            continue
        pairs = sorted(((ip[p][2] @ ip[q][2], p, q) for p in range(len(ip)) for q in range(p + 1, len(ip))))
        free = set(range(len(ip)))
        for _, p, q in pairs:
            if p in free and q in free:
                free -= {p, q}
                turns.append((nearest_slot(ip[p][0], ip[p][1], sgn * Y),
                              nearest_slot(ip[q][0], ip[q][1], sgn * Y), None))
used = [s for t in turns for s in t[:2]]
assert len(used) == len(set(used)), "a fibre slot was used twice"
n_free = len(slot_pt) - len(used)                   # fibre ends left free (lattice boundary)

def hermite(P0, T0, P1, T1, n=300):
    c = 1.2 * np.linalg.norm(P1 - P0)
    t = np.linspace(0, 1, n)[:, None]
    h00, h10, h01, h11 = 2*t**3-3*t**2+1, t**3-2*t**2+t, -2*t**3+3*t**2, t**3-t**2
    return h00*P0 + h10*c*T0 + h01*P1 + h11*c*T1
n_helix = len(curves)                               # curves before this point are beam helices
for s0, s1, way in turns:
    P0, P1 = slot_pt[s0], slot_pt[s1]
    T0, T1 = -slot_tan[s0], slot_tan[s1]            # into the node, then out
    if way is None:
        curves.append((s0, s1, hermite(P0, T0, P1, T1)))
    else:                                           # two arcs through a waypoint round the node
        W, TW = way
        curves.append((s0, s1, np.vstack([hermite(P0, T0, W, TW), hermite(W, TW, P1, T1)[1:]])))

# ---- beads: shared slot beads + resampled interior beads ------------------------------
pos, bonds, slot_id, at_node = [], [], {}, []
for sl, p in slot_pt.items():
    slot_id[sl] = len(pos); pos.append(p); at_node.append(True)
for kc, (sA, sB, P) in enumerate(curves):
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    sc = np.concatenate([[0], np.cumsum(seg)])
    n = max(1, int(round(sc[-1] / ds)))
    ids = [slot_id[sA]]
    for v in np.linspace(0, sc[-1], n + 1)[1:-1]:
        ids.append(len(pos)); at_node.append(kc >= n_helix)
        pos.append(np.array([np.interp(v, sc, P[:, k]) for k in range(3)]))
    ids.append(slot_id[sB])
    bonds += list(zip(ids[:-1], ids[1:]))
pos = np.array(pos); bonds = np.array(bonds); at_node = np.array(at_node)
pos[:, 2] -= pos[:, 2].min(); pos[:, 0] -= pos[:, 0].min(); pos[:, 1] -= pos[:, 1].min()

# ---- fibres = connected components (molecule id) -------------------------------------
parent = list(range(len(pos)))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]; i = parent[i]
    return i
for i, j in bonds:
    parent[find(i)] = find(j)
roots = {r: k + 1 for k, r in enumerate(sorted({find(i) for i in range(len(pos))}))}
mol = np.array([roots[find(i)] for i in range(len(pos))])

# ---- fibre diameter: largest value <= target that keeps different fibres apart --------
tree = cKDTree(pos)
pairs = tree.query_pairs(1.5 * d_target, output_type="ndarray")
pairs = pairs[mol[pairs[:, 0]] != mol[pairs[:, 1]]]
dmin = np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1).min() if len(pairs) else np.inf
d_fibre = min(d_target, 0.97 * dmin)
bl = np.linalg.norm(pos[bonds[:, 0]] - pos[bonds[:, 1]], axis=1)
btype = np.ones(len(bonds), int)
if fuse_nodes:                                      # joints: bonds between touching fibres in the nodes
    fp = tree.query_pairs((1 + fuse_gap) * d_fibre, output_type="ndarray")
    fp = fp[(mol[fp[:, 0]] != mol[fp[:, 1]]) & at_node[fp[:, 0]] & at_node[fp[:, 1]]]
    bonds = np.vstack([bonds, fp]); btype = np.concatenate([btype, 2 * np.ones(len(fp), int)])
nbt = btype.max()
maxb = np.bincount(bonds.ravel(), minlength=len(pos)).max()
print(f"re-entrant cells {Nx}x{Ny}x{Nz} (x y z), {Ny + 1} layers in y; rib {h}, strut {l}, "
      f"cell {h} x {2*lc:.2f}, theta={theta_deg:.2f} deg")
print(f"beam thickness {2*R + d_fibre:.3f}, R_eff={R:.3f}, bead spacing {ds:.3f}, n_rev={n_quarter}/4, fibre diameter {d_fibre:.3f} "
      f"(target {d_target:.3f}; closest different-fibre beads {dmin:.3f})")
if d_fibre < 0.6 * R:
    print("  WARNING: fibres very thin - fibres cross somewhere; change t_beam, d_ratio or lengths")
print(f"beads {len(pos)}, bonds {len(bonds)}, fibres {len(roots)}, free fibre ends {n_free}, beams {nb}, "
      f"bond length {bl.min():.3f}-{bl.max():.3f}")
if fuse_nodes:
    print(f"fused node joints: {int((btype == 2).sum())} bonds between touching fibres (bond type 2)")
print(f"lattice size {np.ptp(pos[:,0]):.1f} x {np.ptp(pos[:,1]):.1f} x {np.ptp(pos[:,2]):.1f} (x y z)")

# ---- beam coefficients (one bond type) --------------------------------------------------
a_r = d_fibre / 2
A, I, J = np.pi * a_r**2, np.pi * a_r**4 / 4, np.pi * a_r**4 / 2
G = E / (2 * (1 + nu))
m = rho * A * ds
Irot = 0.4 * m * a_r**2
kr, ks, kt, kb = E * A, 0.9 * G * A, G * J / ds, E * I / ds
Fcr, Fcs, Tct, Tcb = sigma_f * A, 0.6 * sigma_f * A, 0.6 * sigma_f * J / a_r, sigma_f * I / a_r
gn = zeta * 2 * np.sqrt(kr / ds * m / 2)
gs = zeta * 2 * np.sqrt(ks / ds * m / 2)
gr = zeta * 2 * np.sqrt(kb * Irot / 2)
gt = zeta * 2 * np.sqrt(kt * Irot / 2)
with open("woven_lattice.bonds", "w") as f:
    f.write("# bond_coeff type Kr Ks Kt Kb Fcr Fcs Tct Tcb gnorm gslide groll gtwist\n")
    f.write(f"bond_coeff 1 {kr:.6g} {ks:.6g} {kt:.6g} {kb:.6g} {Fcr:.6g} {Fcs:.6g} "
            f"{Tct:.6g} {Tcb:.6g} {gn:.6g} {gs:.6g} {gr:.6g} {gt:.6g}\n")
    if nbt == 2:                                    # fused joints: same material as the fibres
        f.write(f"bond_coeff 2 {kr:.6g} {ks:.6g} {kt:.6g} {kb:.6g} {Fcr:.6g} {Fcs:.6g} "
                f"{Tct:.6g} {Tcb:.6g} {gn:.6g} {gs:.6g} {gr:.6g} {gt:.6g}\n")
    f.write(f"variable dfib equal {d_fibre:.6g}\n")
    # contact stiffness scaled with the fibre, and a stable timestep for beads of this size
    kn0 = 2.2 * E * d_fibre
    w = max(np.sqrt(kn0 / (m / 2)), np.sqrt(kr / ds / (m / 2)), np.sqrt(kb / (Irot / 2)), np.sqrt(kt / (Irot / 2)))
    f.write(f"variable kn0 equal {kn0:.6g}\n")
    f.write(f"variable dtauto equal {0.07 / w:.4g}\n")

# ---- data files ---------------------------------------------------------------------------
lo, hi = pos.min(0) - 3.0, pos.max(0) + 3.0
dens = m / (np.pi / 6 * d_fibre**3)                  # bead mass = rho*A*ds
def write(fname, style):
    with open(fname, "w") as f:
        f.write(f"woven re-entrant lattice {Nx}x{Ny}x{Nz} theta={theta_deg:.2f} t={t_beam} n_rev={n_quarter}/4\n\n")
        f.write(f"{len(pos)} atoms\n{len(bonds)} bonds\n1 atom types\n{nbt} bond types\n")
        if style == "bpm":
            f.write(f"{max(6, maxb + 2)} extra bond per atom\n{36 if nbt == 1 else 200} extra special per atom\n")
        f.write("\n")
        for k, ax in enumerate("xyz"):
            f.write(f"{lo[k]:.6f} {hi[k]:.6f} {ax}lo {ax}hi\n")
        if style == "bond":
            f.write(f"\nMasses\n\n1 {m:.6g}\n")
        f.write(f"\nAtoms # {'bpm/sphere' if style == 'bpm' else 'bond'}\n\n")
        for i, (p, mo) in enumerate(zip(pos, mol)):
            if style == "bpm":
                f.write(f"{i+1} {mo} 1 {d_fibre:.6g} {dens:.6g} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            else:
                f.write(f"{i+1} {mo} 1 {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        f.write("\nBonds\n\n")
        f.writelines(f"{k+1} {t} {i+1} {j+1}\n" for k, ((i, j), t) in enumerate(zip(bonds, btype)))
write("woven_lattice.data", "bpm")
write("woven_lattice_ovito.data", "bond")
print("wrote woven_lattice.data, woven_lattice.bonds, woven_lattice_ovito.data")
