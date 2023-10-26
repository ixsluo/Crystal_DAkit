import argparse
import logging

import cdakit
import cdakit.log
import cdakit.prepare_vasp
import cdakit.find_spg
import cdakit.match_structure


def main(verbose: int, **kwargs):
    cdakit.log.init(20 - 10 * verbose)
    func = kwargs.get("func", None)
    logging.debug(f"{kwargs}")
    if func is not None:
        func(**kwargs)


def cli():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-j", "--njobs", type=int, default=1, help="n core to parallel")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    # create subparsers
    subparsers = parser.add_subparsers()
    # add subparser to subparsers
    cdakit.prepare_vasp.add_subparser(subparsers)
    cdakit.find_spg.add_subparser(subparsers)
    cdakit.match_structure.add_subparser(subparsers)
    # parse
    args, unknown_args = parser.parse_known_args()
    # calling main
    main(**vars(args))


if __name__ == "__main__":
    cli()

