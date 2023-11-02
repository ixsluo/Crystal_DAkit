import argparse
import logging
import pickle
import warnings
from itertools import chain
from pathlib import Path

import pandas as pd
from ase.formula import Formula
from joblib import Parallel, delayed
from tqdm import tqdm

from cdakit.iotools import to_format_table
from cdakit.log import logit

logger = logging.getLogger(__name__)


def parse_one_outcar(foutcar: Path) -> pd.DataFrame:
    energylist = []  # eV
    Vlist = []
    PVlist = []  # eV
    extpres = []  # kbar
    converge = False
    cputime = pd.NA
    natoms = pd.NA
    formula = pd.NA
    with open(foutcar, "r") as f:
        for line in f:
            if "NIONS" in line:
                natoms = int(line.strip().split()[-1])
            elif "POSCAR =" in line:
                formula = Formula(line.strip()[8:].replace(" ", "")).format("metal")
            elif "energy  without" in line:
                energylist.append(float(line.strip().split()[-1]))
            elif "P V=" in line:
                PVlist.append(float(line.strip().split()[-1]))
            elif "volume of cell" in line:
                Vlist.append(line.strip().split()[-1])
            elif "external pressure" in line:
                extpres.append(line.strip().split()[3])
            elif "reached required" in line:
                converge = True
            elif "CPU" in line:
                cputime = float(line.strip().split()[-1])
    if len(energylist) == 0:
        energylist = [pd.NA]
        Vlist = [pd.NA]
        extpres = [pd.NA]
    else:
        Vlist = Vlist[1:]  # drop the duplicated first one
        Vlist = Vlist[:len(energylist)]
    if len(PVlist) == 0:
        PVlist = [0] * len(energylist)
    else:
        PVlist = PVlist[:len(energylist)]
    extpres = extpres[:len(energylist)]
    convergelist = [False] * len(energylist)
    convergelist[-1] = converge
    cputime = [cputime] * len(energylist)
    natoms = [natoms] * len(energylist)
    formula = [formula] * len(energylist)
    parsed_df = pd.DataFrame(
        {
            "formula": formula,
            "energy": energylist,
            "volume": Vlist,
            "PV": PVlist,
            "extpressure": extpres,
            "converge": convergelist,
            "cputime": cputime,
            "natoms": natoms,
            "nsites": natoms,
        }
    )
    parsed_df["enthalpy"] = parsed_df["energy"] + parsed_df["PV"]
    return parsed_df


def stat_outcar_dfdict(dfdict: dict[str, pd.DataFrame]) -> pd.DataFrame:
    serlist = []
    for fname, df in dfdict.items():
        ser = pd.Series(
            {
                "formula": df.at[0, "formula"],
                "converge": df.converge.iloc[-1],
                "decreased_enth": df.at[0, "enthalpy"] - df.at[len(df) - 1, "enthalpy"],
                "ion_steps": len(df),
                "natoms": df.at[0, "natoms"],
                "nsites": df.at[0, "nsites"],
            },
            name=fname
        )
        serlist.append(ser)
    stat_df = pd.DataFrame(serlist)
    stat_df["decreased_enth_per_atom"] = stat_df["decreased_enth"] / stat_df["natoms"]
    return stat_df


@logit()
def parse_outcar(indir, njobs, *args, **kwargs):
    indir = Path(indir)
    outcars = list(chain(indir.rglob("OUTCAR"), indir.rglob("*.OUTCAR")))
    if len(outcars) == 0:
        raise ValueError("No OUTCAR or *.OUTCAR found")
    parsed_dflist = Parallel(njobs, backend="multiprocessing")(
        delayed(parse_one_outcar)(foutcar)
        for foutcar in tqdm(outcars, ncols=120)
    )
    parsed_dfdict = {
        str(foutcar.relative_to(indir)): parsed_df
        for foutcar, parsed_df in zip(outcars, parsed_dflist)
    }
    stat_df = stat_outcar_dfdict(parsed_dfdict)
    print(stat_df)
    with open(indir.joinpath("parsed_outcar.pkl"), "wb") as f:
        pickle.dump(parsed_dfdict, f)
    with open(indir.joinpath("parsed_outcar.table"), "w") as f:
        f.write(to_format_table(stat_df))


def add_subparser(subparsers):
    subparser = subparsers.add_parser(
        __name__.split(".")[-1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Parse all OUTCAR and *.OUTCAR, return list of DataFrame
to parsed_outcar.pkl and summary to parsed_outcar.table"""
    )
    subparser.set_defaults(func=parse_outcar)
    subparser.add_argument("indir", help="directory containing OUTCAR or *.OUTCAR")
