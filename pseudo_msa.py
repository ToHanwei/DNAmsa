#!coding: utf-8

import os
import sys

from bin.commandparse import Argparse
from bin.CDSalignment import makedict
from bin.CDSalignment import buildcode
from bin.CDSalignment import mapseq
from bin.CDSalignment import pseudomap
from bin.CDSalignment import write2file


def commands():
    """
    command lin eparser
    """
    if len(sys.argv) < 2:
        # no command parma
        os.system("python pseudo_msa.py -h")
        sys.exit()
    ArgParse = Argparse()
    ArgParse.pseudo_parse()
    args = ArgParse.pseudo_args
    msafile = args.msafile
    cdsfile = args.cdsfile
    codefile = args.transcode
    outfile = args.output
    return (
        msafile, cdsfile, 
        codefile, outfile
        )


def main():
    msafile, cdsfile, codefile, outfile = commands()
    prot_dict = makedict(msafile)
    uncl_dict = makedict(cdsfile)
    code_dict = buildcode(codefile)
    resu_dict = pseudomap(prot_dict, uncl_dict, code_dict)
    write2file(resu_dict, outfile)


if __name__ == "__main__":
    main()


