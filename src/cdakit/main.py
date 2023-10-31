import argparse
import multiprocessing
import logging

import cdakit
import cdakit.find_spg
import cdakit.match_structure
import cdakit.prepare_calypso
import cdakit.prepare_vasp
import cdakit.standardize
from cdakit.log import loglistener


def main(verbose: int, **kwargs):
    queue = multiprocessing.Queue(-1)
    kwargs["logqueue"] = queue
    kwargs["logconf"] = {"level": 20 - 10 * verbose}
    # log listener
    logproc = multiprocessing.Process(target=loglistener, args=(queue,))
    logproc.start()
    # worker
    func = kwargs.get("func", None)
    if func is not None:
        worker = multiprocessing.Process(target=func, kwargs=kwargs)
        worker.start()
        worker.join()
    queue.put_nowait(None)
    logproc.join()


def cli():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-j", "--njobs", type=int, default=1, help="n core to parallel")
    parser.add_argument("-v", "--verbose", action="count", default=0)
    # create subparsers
    subparsers = parser.add_subparsers()
    # add subparser to subparsers
    cdakit.find_spg.add_subparser(subparsers)
    cdakit.match_structure.add_subparser(subparsers)
    cdakit.prepare_calypso.add_subparser(subparsers)
    cdakit.prepare_vasp.add_subparser(subparsers)
    cdakit.standardize.add_subparser(subparsers)
    # parse
    args, unknown_args = parser.parse_known_args()
    # calling main
    main(**vars(args))


if __name__ == "__main__":
    cli()

