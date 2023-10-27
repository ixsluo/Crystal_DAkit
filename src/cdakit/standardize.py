import argparse
import logging
from pathlib import Path

from ase import Atoms
from ase.io import read, write
from spglib import standardize_cell

from cdakit.log import logit


logger = logging.getLogger(__name__)


@logit()
def standardize(vaspfile, symprec, *args, **kwargs):
    vaspfile = Path(vaspfile)
    atoms = read(vaspfile, format="vasp")
    bulk = (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    for isymprec in symprec:
        stddir = Path(vaspfile).parent.joinpath(f"{vaspfile}.std/{isymprec}")
        stddir.mkdir(parents=True, exist_ok=True)
        ucell = standardize_cell(bulk, False, symprec=isymprec)
        pcell = standardize_cell(bulk, True, symprec=isymprec)
        for celltag, stdcell in zip(["ucell", "pcell"], [ucell, pcell], strict=True):
            if stdcell is None:
                logger.warning(f"{vaspfile} cannot find standard {celltag} under symprec={isymprec}, using itself to replace")
                stdcell = atoms
            else:
                lattice, scaled_positions, numbers = stdcell
                stdcell = Atoms(numbers, cell=lattice, scaled_positions=scaled_positions)
            write(stddir.joinpath(f"{vaspfile.stem}.{celltag}.vasp"), stdcell)


def add_subparser(subparsers):
    subparser = subparsers.add_parser(
        __name__.split(".")[-1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="match structures to target by pymatgen, write to match.<target>.table",
    )
    subparser.set_defaults(func=standardize)
    subparser.add_argument("vaspfile", help="vaspfile to analysis, recommanded to name as *.vasp")
    subparser.add_argument("-s", "--symprec", type=float, nargs="+", default=[0.5, 0.1, 0.01], help="symprec tolerence")
