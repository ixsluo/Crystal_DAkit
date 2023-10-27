# find the symmetry and standardlized cell of a given dir

import argparse
import io
import subprocess
import sys
import shutil
from contextlib import redirect_stdout
from pathlib import Path

import pandas as pd
import spglib
from ase import Atoms
from ase.io import read, write
from joblib import Parallel, delayed
from tqdm import tqdm

from cdakit.iotools import to_format_table
from cdakit.log import logit


def atoms2cif(atoms):
    with io.BytesIO() as buffer, redirect_stdout(buffer):
        write('-', atoms, format='cif')
        cif = buffer.getvalue().decode()  # byte to string
    return cif


def atoms2vasp(atoms):
    with io.StringIO() as buffer, redirect_stdout(buffer):
        write('-', atoms, format='vasp')
        vasp = buffer.getvalue()  # byte to string
    return vasp


def get_spg_one(name: Path, atoms, symprec_list, angle_tolerance=10):
    spg_dict = {
        "name": name.name,
        "formula": atoms.get_chemical_formula("metal"),
    }
    cell = (atoms.cell[:], atoms.get_scaled_positions(), atoms.get_atomic_numbers())
    for symprec in symprec_list:
        symds = spglib.get_symmetry_dataset(cell, symprec, angle_tolerance)
        # ---- record
        if symds is not None:
            std_atoms = Atoms(
                symds["std_types"],
                cell=symds["std_lattice"],
                scaled_positions=symds["std_positions"],
            )
            spg_dict["{:.0e}".format(symprec)] = symds['number']
            spg_dict["{:.0e}".format(symprec) + "_symbol"] = symds['international']
            spg_dict["{:.0e}".format(symprec) + "_std_natoms"] = len(symds['std_types'])
            spg_dict["{:.0e}".format(symprec) + "_std_cif"] = atoms2cif(std_atoms)
            spg_dict["{:.0e}".format(symprec) + "_std_vasp"] = atoms2vasp(std_atoms)
        else:
            print(name, symprec, "Cannot find symmetry", file=sys.stderr)
            spg_dict["{:.0e}".format(symprec)] = 0
            spg_dict["{:.0e}".format(symprec) + "_symbol"] = "-"
            spg_dict["{:.0e}".format(symprec) + "_std_natoms"] = 0
            spg_dict["{:.0e}".format(symprec) + "_std_cif"] = atoms2cif(atoms)
            spg_dict["{:.0e}".format(symprec) + "_std_vasp"] = atoms2vasp(atoms)
    # ---- filter P1
    if any(spg_dict["{:.0e}".format(symprec)] > 1 for symprec in symprec_list):
        sympart = name.parent.parent.joinpath("sympart/gen")
        sympart.mkdir(exist_ok=True, parents=True)
        shutil.copy(name, sympart / name.name)
    return pd.Series(spg_dict)


def get_spg_df(fdir, symprec_list=(0.5, 0.1, 0.01)):
    flist = list(Path(fdir).glob("*.vasp"))
    symprec_list = sorted(symprec_list, reverse=True)
    ser_list = Parallel(-1, backend="multiprocessing")(
        delayed(get_spg_one)(f, read(f), symprec_list)
        for f in tqdm(flist, ncols=180, desc=f"{fdir}")
    )
    df = pd.DataFrame(ser_list)
    df = df.sort_values(by=list(map("{:.0e}".format, symprec_list)), ascending=False)
    return df


def write_std_vasp(df: pd.DataFrame, indir):
    cols = [col for col in df if col.endswith("_std_vasp")]
    prec_list = [col[:-9] for col in cols]
    for prec in prec_list:
        prec_dir = Path(indir).with_name(f"std_{prec}")
        prec_dir.mkdir(exist_ok=True)
        for _, ser in df.iterrows():
            with open(prec_dir / Path(ser['name']).name, "w") as fvasp:
                fvasp.write(ser[prec + "_std_vasp"])


@logit()
def find_spg(indirs, symprec, **kwargs):
    spgdfdict = {}
    for indir in indirs:
        df = get_spg_df(indir, symprec)
        table = to_format_table(
            df[[col for col in df if not col.endswith(("_std_cif", "_std_vasp"))]]
        )
        with open(Path(indir).with_name("spg.txt"), 'w') as f:
            f.write(table)
        write_std_vasp(df, indir)
        spgdfdict[Path(indir).parent.name] = df


def add_subparser(subparsers):
    subparser = subparsers.add_parser(
        __name__.split(".")[-1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.set_defaults(func=find_spg)
    subparser.add_argument("indirs", nargs="*", help="directiries containing *.vasp")
    subparser.add_argument("-s", "--symprec", type=float, default=[0.5, 0.1, 0.01], nargs=-1, help="symprec, only one significant digits is kept")

