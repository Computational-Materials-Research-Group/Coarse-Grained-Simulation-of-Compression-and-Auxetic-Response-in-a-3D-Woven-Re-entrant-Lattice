"""
make_lattice_session.py - ready-made OVITO sessions for the woven re-entrant
lattice compression run.

    python make_lattice_session.py            # writes the sessions
    python make_lattice_session.py render     # also renders first/last frame PNGs (front view)

Sessions (open in OVITO with File > Load Session State, then press Play):
  lattice_view.ovito    fibres as solid blue tubes, broken bonds removed
  lattice_fibres.ovito  every fibre in its own colour (shows the weave)
  lattice_stress.ovito  fibres coloured by axial stress sigma_zz (blue = compression)

Needs in the same folder: woven_lattice_ovito.data, woven_lattice.bonds,
lattice.dump, lattice_bonds.dump.   pip install ovito
"""
import sys, re
import numpy as np
import ovito
from ovito.io import import_file
from ovito.modifiers import (LoadTrajectoryModifier, AssignColorModifier, ColorCodingModifier,
                             ComputePropertyModifier)
from ovito.vis import Viewport, TachyonRenderer

BLUE = (0.25, 0.42, 0.85)
dfib = float(re.search(r"variable dfib equal ([0-9.eE+-]+)", open("woven_lattice.bonds").read()).group(1))

def base():
    p = import_file("woven_lattice_ovito.data", atom_style="bond")
    t = LoadTrajectoryModifier(); p.modifiers.append(t)
    t.source.load("lattice.dump")                                   # positions, c_nb, szz
    b = LoadTrajectoryModifier(); p.modifiers.append(b)
    b.source.load("lattice_bonds.dump",                             # bonds still intact
                  columns=["Bond Index", "Particle Identifiers.1", "Particle Identifiers.2", "Bond Type"])
    d = p.source.data
    d.particles.vis.radius = dfib / 2
    d.particles.bonds.vis.width = dfib                               # beads + bonds = smooth tubes
    d.particles.bonds.vis.use_particle_colors = True
    d.cell.vis.enabled = False
    return p, t

def reset():
    for q in list(ovito.scene.pipelines):
        q.remove_from_scene()

# 1) uniform blue tubes
p, traj = base()
p.modifiers.append(AssignColorModifier(color=BLUE))
p.add_to_scene()
ovito.scene.save("lattice_view.ovito")
nfr = traj.source.num_frames

# 2) one colour per fibre
reset()
p, _ = base()
p.modifiers.append(ColorCodingModifier(property="Molecule Identifier",
                                       gradient=ColorCodingModifier.Rainbow()))
p.add_to_scene()
ovito.scene.save("lattice_fibres.ovito")

# 3) axial stress per bead (stress/atom is stress x volume; divide by bead volume)
reset()
p, _ = base()
vb = np.pi / 6 * dfib**3
p.modifiers.append(ComputePropertyModifier(output_property="sigma_zz",
                                           expressions=[f"szz / {vb:.6g}"]))
last = p.compute(nfr - 1).particles["sigma_zz"]
lim = float(np.percentile(np.abs(last), 98)) or 1e-6
p.modifiers.append(ColorCodingModifier(property="sigma_zz", start_value=-lim, end_value=lim,
                                       gradient=ColorCodingModifier.BlueWhiteRed()))
p.add_to_scene()
ovito.scene.save("lattice_stress.ovito")
print(f"wrote lattice_view.ovito, lattice_fibres.ovito, lattice_stress.ovito ({nfr} frames, "
      f"stress colour range +-{lim:.3g})")

if len(sys.argv) > 1 and sys.argv[1] == "render":
    reset()
    p, _ = base(); p.modifiers.append(AssignColorModifier(color=BLUE)); p.add_to_scene()
    vp = Viewport(type=Viewport.Type.Ortho, camera_dir=(0, 1, 0))   # front view (x-z plane)
    vp.zoom_all((1400, 1600))
    for fr, name in [(0, "lattice_first.png"), (nfr - 1, "lattice_last.png")]:
        vp.render_image(filename=name, size=(1400, 1600), frame=fr, background=(1, 1, 1),
                        renderer=TachyonRenderer(ambient_occlusion=True, antialiasing_samples=8))
        print("rendered", name)
