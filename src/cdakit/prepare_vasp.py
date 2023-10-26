import argparse
import logging
import warnings
from pathlib import Path

from joblib import Parallel, delayed
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import VaspInput
from pymatgen.io.vasp.sets import MPRelaxSet
from tqdm import tqdm

from cdakit.iotools import read_format_table


logger = logging.getLogger(__name__)


def prepare_task(structure, relax_path, vaspargs):
    user_incar_settings = {
        'LREAL': False,
        'ISMEAR': 0,
        'NCORE': 4,
        'NSW': vaspargs["nsw"],
        'PSTRESS': vaspargs["pstress"],
        'ISYM': vaspargs["sym"],
    }
    if vaspargs["ediff"] is not None:
        user_incar_settings["EDIFF"] = vaspargs["ediff"]
    if vaspargs["ediffg"] is not None:
        user_incar_settings["EDIFFG"] = vaspargs["ediffg"]
    if vaspargs["kspacing"] is not None:
        user_incar_settings["KSPACING"] = vaspargs["kspacing"]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mp_set = MPRelaxSet(
            structure,
            user_incar_settings=user_incar_settings,
            user_potcar_settings={"W": "W_sv"},
            user_potcar_functional="PBE_54",
        )
        vasp = VaspInput(
            incar=mp_set.incar,
            kpoints=mp_set.kpoints,
            poscar=mp_set.poscar,
            potcar=mp_set.potcar,
        )
        vasp.write_input(relax_path)


def wrapped_prepare_task(indir, uniq, uniqlevel, sf, vaspargs):
    runtype = ".scf" if vaspargs["nsw"] <= 1 else ".opt"
    if uniq is not None:
        runtype = f".uniq.{uniqlevel}" + runtype
    relax_path = indir.with_suffix(runtype).joinpath(sf.stem)
    relax_path.mkdir(exist_ok=True, parents=True)

    structure = Structure.from_file(sf)
    prepare_task(structure, relax_path, vaspargs)


def prepare_vasp_batch(
    indir, uniqfile, uniqlevel,njobs, ediff, ediffg, nsw, pstress, kspacing, sym, **kwargs
):
    vaspargs = {"ediff": ediff, "ediffg": ediffg, "nsw": nsw,
               "pstress": pstress, "kspacing": kspacing, "sym": sym}
    logger.info("You are using " + " ".join(f"{k}={v}" for k, v in vaspargs.items()))
    logger.warning("W POTCAR is replaced by W_sv")
    indir = Path(indir)
    flist = list(indir.glob("*.vasp"))
    if uniqfile is not None:
        lv = f"matcher_{uniqlevel}"
        uniqdf = read_format_table(uniqfile)
        if lv not in uniqdf.columns:
            raise KeyError(f"key '{uniqlevel}' not in {uniqfile}")
        click.echo(f"using unique key '{lv}' in {uniqfile}")
        uniqlist = list(uniqdf[uniqdf[lv]].index)
        flist = [fi for fi in flist if int(fi.stem) in uniqlist]
    Parallel(njobs, backend="multiprocessing")(
        delayed(wrapped_prepare_task)(indir, uniqfile, uniqlevel, sf, vaspargs)
        for sf in tqdm(flist, ncols=120)
    )


def add_subparser(subparsers):
    subparser = subparsers.add_parser(
        __name__.split(".")[-1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparser.set_defaults(func=prepare_vasp_batch)
    subparser.add_argument("indir", help="directory containing *.vasp")
    subparser.add_argument("-u", "--uniqfile",                                            help="unique file to read")
    subparser.add_argument("-l", "--uniqlevel", choices=["lo", "md", "st"], default="lo", help="unique level of matcher used in uniqfile")
    subparser.add_argument("-e", "--ediff", type=float,                                   help="EDIFF, autogenerate by pyamtgen if None")
    subparser.add_argument("-eg", "--ediffg", type=float,                                 help="EDIFFG")
    subparser.add_argument("-n", "--nsw", type=float, default=0,                          help="NSW")
    subparser.add_argument("-p", "--pstress", type=float, default=0,                      help="PSTRESS(kbar)")
    subparser.add_argument("-ks", "--kspacing",                                           help="KSPACING")
    subparser.add_argument("-s", "--sym", type=int, default=0,                            help="ISYM, suggest 0/2")

