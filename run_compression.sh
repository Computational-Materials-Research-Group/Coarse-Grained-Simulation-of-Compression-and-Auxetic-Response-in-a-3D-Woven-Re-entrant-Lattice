#!/usr/bin/env bash
# run_compression.sh - build the woven re-entrant 4x4x4 lattice and run the
# 50 % compression test with LAMMPS on 14 MPI processes (Ubuntu / WSL).
#
#   cd "/mnt/d/Atomisitc modeling research/woven_reentrant_4x4x4"
#   chmod +x run_compression.sh
#   ./run_compression.sh                                  # run in this terminal
#   nohup ./run_compression.sh > run.out 2>&1 &           # or: keep running after closing it
#
# Options (environment variables):
#   NP=14         number of MPI processes
#   LMP=lmp       LAMMPS executable name (lmp, lmp_mpi, or a full path)
#   EXTRA="-var erate 5e-5"   extra LAMMPS command-line settings

set -euo pipefail

NP=${NP:-14}
LMP=${LMP:-lmp}
EXTRA=${EXTRA:-}
cd "$(dirname "$0")"
echo "== folder: $(pwd)"

# ---- 1. checks ---------------------------------------------------------------
if ! command -v "$LMP" >/dev/null 2>&1; then
    for c in lmp_mpi lmp_openmpi lmp_mpich lmp; do
        command -v "$c" >/dev/null 2>&1 && LMP=$c && break
    done
fi
command -v "$LMP" >/dev/null 2>&1 || { echo "ERROR: LAMMPS executable not found (set LMP=/path/to/lmp)"; exit 1; }
command -v mpirun >/dev/null 2>&1 || { echo "ERROR: mpirun not found: sudo apt install openmpi-bin libopenmpi-dev"; exit 1; }
echo "== LAMMPS: $(command -v "$LMP")"

PKGS=$("$LMP" -h 2>/dev/null | awk '/Installed packages:/{f=1;next} f&&/^List of/{exit} f' || true)
for p in BPM GRANULAR MOLECULE; do
    if ! grep -qw "$p" <<< "$PKGS"; then
        echo "ERROR: LAMMPS package $p is not installed in $LMP (see build notes at the end of this file)"
        exit 1
    fi
done
echo "== packages BPM, GRANULAR, MOLECULE: ok"

CORES=$(nproc)
OVERSUB=""
if [ "$NP" -gt "$CORES" ]; then
    echo "WARNING: $NP processes but only $CORES cores visible; it will run, but slower"
    mpirun --version 2>/dev/null | grep -qi "open mpi" && OVERSUB="--oversubscribe"
fi

# Python: use a local virtual environment (Ubuntu blocks pip install into the system Python)
PY=python3
if ! python3 -c "import numpy, scipy, matplotlib" 2>/dev/null; then
    if [ ! -x .venv/bin/python ]; then
        echo "== creating .venv with numpy scipy matplotlib ovito"
        python3 -m venv .venv || { echo "ERROR: sudo apt install python3-venv"; exit 1; }
        .venv/bin/pip install -q numpy scipy matplotlib ovito
    fi
    PY=.venv/bin/python
fi

# ---- 2. build the lattice ------------------------------------------------------
echo "== building lattice"
"$PY" build_reentrant_woven.py | tee build.log

# ---- 3. run LAMMPS -----------------------------------------------------------------
export OMP_NUM_THREADS=1
echo "== running: mpirun -np $NP $LMP -in in.compression_lattice $EXTRA"
echo "   started $(date); progress: tail -f log.lammps   (strain = 2nd column)"
mpirun $OVERSUB -np "$NP" "$LMP" -in in.compression_lattice -log log.lammps $EXTRA
echo "== LAMMPS finished $(date)"

# ---- 4. post-processing ----------------------------------------------------------
"$PY" plot_compression.py || echo "WARNING: plotting failed"
if "$PY" -c "import ovito" 2>/dev/null; then
    QT_QPA_PLATFORM=offscreen "$PY" make_lattice_session.py || echo "WARNING: OVITO session script failed"
else
    echo "NOTE: ovito not installed for $PY - run make_lattice_session.py on Windows instead"
fi
echo "== done: stress_strain.txt, compression_results.png, lattice.dump, lattice_bonds.dump, lattice_*.ovito"

# ---------------------------------------------------------------------------------------
# Build notes - LAMMPS with MPI + BPM + GRANULAR + MOLECULE on Ubuntu (once):
#
#   sudo apt update
#   sudo apt install -y build-essential cmake git openmpi-bin libopenmpi-dev python3-venv
#   git clone -b stable --depth 1 https://github.com/lammps/lammps.git ~/lammps
#   cd ~/lammps && mkdir build && cd build
#   cmake ../cmake -D BUILD_MPI=on -D PKG_BPM=on -D PKG_GRANULAR=on -D PKG_MOLECULE=on \
#                  -D PKG_EXTRA-FIX=on -D CMAKE_BUILD_TYPE=Release
#   make -j 14
#   sudo cp lmp /usr/local/bin/lmp
# ---------------------------------------------------------------------------------------
