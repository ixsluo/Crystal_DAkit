import argparse
import logging
import os
import warnings
import shutil
import subprocess
from itertools import product
from pathlib import Path

import numpy as np
from pymatgen.io.vasp import Poscar, Potcar
from joblib import Parallel, delayed

from cdakit.log import logit


logger = logging.getLogger(__name__)


def potcar2distmat(potcar: Potcar):
    rc = [p.RCORE for p in potcar]
    distmat = np.asarray([(i + j) * 0.529177 for i, j in product(rc, repeat=2)])
    distmat = distmat.reshape(len(rc), len(rc))
    return distmat


def vasp2inputdat(poscar, potcar, dist_ratio, popsize):
    distmat = potcar2distmat(potcar) * dist_ratio
    ds = ""
    for line in distmat:
        ds += " ".join(map(str, line)) + "\n"

    inputdat = (
        f"SystemName = {''.join(poscar.site_symbols)}\n"
        f"NumberOfSpecies = {len(poscar.site_symbols)}\n"
        f"NameOfAtoms = {' '.join(poscar.site_symbols)}\n"
        f"NumberOfAtoms = {' '.join(map(str, poscar.natoms))}\n"
        "NumberOfFormula = 1 1\n"
        f"Volume = {poscar.structure.volume}\n"
        "@DistanceOfIon\n"
        f"{ds}"
        "@End\n"
        "Ialgo = 2\n"
        "PsoRatio = 0.6\n"
        f"PopSize = {popsize}\n"
        "ICode = 15\n"
        "NumberOfLbest = 4\n"
        "NumberOfLocalOptim = 3\n"
        "Command = sh submit.sh\n"
        "MaxStep = 5\n"
        "PickUp = F\n"
        "PickStep = 5\n"
        "Parallel = F\n"
        "Split = T\n"
    )
    return inputdat


@logit()
def prepare_calypso(njobs, indir, dist_ratio, popsize, calypsocmd, calypsotimeout, **kwargs):
    indir = Path(indir)
    Parallel(njobs, backend="multiprocessing")(
        delayed(prepare_calypso_one)(indir, fposcar, dist_ratio, popsize, calypsocmd, calypsotimeout)
        for fposcar in indir.rglob("POSCAR")
    )


def prepare_calypso_one(indir, fposcar, dist_ratio, popsize, calypsocmd, calypsotimeout):
    outdir = indir.with_name(f"{indir.name}.calypso")
    logger.info(f"Processing {fposcar.parent}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        poscar = Poscar.from_file(fposcar)
        potcar = Potcar.from_file(fposcar.with_name("POTCAR"))
    calypsodir = outdir.joinpath(fposcar.parent.relative_to(indir))
    logger.debug(f"calypsodir: {str(calypsodir)}")
    calypsodir.mkdir(parents=True, exist_ok=True)
    with open(calypsodir.joinpath("input.dat"), 'w') as f:
        f.write(vasp2inputdat(poscar, potcar, dist_ratio, popsize))
    for f in ["INCAR", "POTCAR", "KPOINTS"]:
        try:
            shutil.copy(fposcar.with_name(f), calypsodir)
        except FileNotFoundError:
            logger.warning(f"{str(fposcar.with_name(f))} not found")
    # run calypso
    try:
        os.remove(calypsodir.joinpath("step"))
    except FileNotFoundError:
        pass
    with open(calypsodir.joinpath("caly.log"), "w") as calylog:
        try:
            proc = subprocess.run(
                calypsocmd, stdout=calylog, stderr=subprocess.STDOUT, cwd=calypsodir,
                timeout=calypsotimeout,
            )
        except subprocess.TimeoutExpired:
            logger.error(f"CALYPSO failed in {str(calypsodir)}")
        else:
            # split POSCAR_* to subdir
            if proc.returncode != 0:
                logger.error(f"Calling {calypsocmd} failed in {calypsodir}")
            else:
                for popi in range(1, popsize + 1):
                    calcdir = calypsodir.joinpath(f"calc/{popi}")
                    calcdir.mkdir(parents=True, exist_ok=True)
                    shutil.move(calypsodir / f"POSCAR_{popi}", calcdir / "POSCAR")
                    for f in ["INCAR", "POTCAR", "KPOINTS"]:
                        try:
                            shutil.copy(calypsodir / f, calcdir)
                        except FileNotFoundError:
                            pass
        finally:
            # clean dir
            for pyfile in calypsodir.glob("*.py"):
                os.remove(pyfile)


def add_subparser(subparsers):
    subparser = subparsers.add_parser(
        __name__.split(".")[-1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="generate input.dat follow each POSCAR in subdir",
    )
    subparser.set_defaults(func=prepare_calypso)
    subparser.add_argument("indir", help="directory with each POSCAR/POTCAR/KPOINTS in a subdir")
    subparser.add_argument("-r", "--dist_ratio", type=float, default=0.7, help="distance ratio multipied on RCORE to generate DistanceOfIon")
    subparser.add_argument("-p", "--popsize", type=int, default=10, help="PopSize")
    subparser.add_argument("-c", "--calypsocmd", default="calypso.x", help="CALYPSO executable file")
    subparser.add_argument("-t", "--calypsotimeout", type=float, default=180, help="maxtime for each calypso subprocess")
