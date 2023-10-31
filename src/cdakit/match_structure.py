import argparse
import logging
from pathlib import Path

import pandas as pd
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.structure import Structure

from cdakit.iotools import to_format_table
from cdakit.log import logit


logger = logging.getLogger(__name__)


def get_matchers():
    matcher_lo = StructureMatcher(ltol=0.3, stol=0.5, angle_tol=10)  # loose
    matcher_md = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5)  # midium
    matcher_st = StructureMatcher(ltol=0.1, stol=0.2, angle_tol=5)  # strict
    matchers = {
        "matcher_lo": matcher_lo,
        "matcher_md": matcher_md,
        "matcher_st": matcher_st,
    }
    return matchers


# match *.vasp with ground-truth structure(gtst) with each matcher in matchers
# calculate average rms distance if matcher
# return
#   matcher_lo matcher_lo_avgd matcher_md matcher_md_avgd matcher_st matcher_st_avgd
# 0        T/F        <float>       T/F          <float>       T/F          <float>
# 1        ...           ...        ...             ...        ...             ...
def match_structure(
    indir: Path,
    gtst: Structure,
    matchers: dict[str, StructureMatcher],
    label,
):
    f_target = indir.with_name(f"{label}.vasp")
    f_matchtable = indir.with_name(f"match.{label}.table")
    try:
        idxlist = sorted([int(f.stem) for f in indir.glob("*.vasp")])
    except Exception:
        logger.warning("cannot sorted by digit, skip sorting")
        idxlist = [f.stem for f in indir.glob("*.vasp")]

    if f_matchtable.exists():
        df = pd.read_table(f_matchtable, sep=r"\s+", index_col="index")
        if len(df) >= len(idxlist):
            return df

    st_dict = {f"{i}.vasp": Structure.from_file(indir / f"{i}.vasp") for i in idxlist}

    data = {}
    for mat_name, matcher in matchers.items():
        fitdict = {i: matcher.fit(gtst, st) for i, st in st_dict.items()}
        distdict = {
            i: (matcher.get_rms_dist(gtst, st_dict[i]) if fit else (pd.NA, pd.NA))
            for i, fit in fitdict.items()
        }

        data[mat_name] = pd.Series(fitdict)
        data[f"{mat_name}_normrms"] = pd.Series({k: v[0] for k, v in distdict.items()})
        data[f"{mat_name}_maxrms"] = pd.Series({k: v[1] for k, v in distdict.items()})

    df = pd.DataFrame(data)

    gtst.to(str(f_target), fmt="poscar")
    table_str = to_format_table(df)
    with open(f_matchtable, "w") as f:
        f.write(table_str)

    return df


@logit()
def matchtarget(indir, target, **kwargs):
    indir = Path(indir).resolve()
    target = Path(target).resolve()
    targetst = Structure.from_file(target)
    matchers = get_matchers()

    matchdf = match_structure(indir, targetst, matchers, target.stem)
    return matchdf


def add_subparser(subparsers):
    subparser = subparsers.add_parser(
        __name__.split(".")[-1],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="match structures to target by pymatgen, write to match.<target>.table",
    )
    subparser.set_defaults(func=matchtarget)
    subparser.add_argument("indir", help="directory containing *.vasp")
    subparser.add_argument("-t", "--target", required=True, help="target structure in vasp format")

