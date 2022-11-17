import os
import re
import sys

import pandas as pd
import numpy as np

from Bio import Align


def makedict(
    infile: str) -> dict:
    """
    build dict from sequences
    """
    seq_dict = {}
    with open(infile) as seqf:
        seqs = seqf.read().split(">")[1:]
    for seq in seqs:
        lines = seq.split("\n")
        name = lines[0]
        seq = "".join(lines[1:])
        seq_dict[name] = seq

    return seq_dict


def buildcode(
    infile: str) -> dict:
    """
    build uncleotide code table
    """
    with open(infile) as codef:
        codes = codef.readlines()
    codes = [
        line.strip().split(",") 
        for line in codes
        ]

    codes = {
        lines[0]:lines[1:] 
        for lines in codes
        }
    return codes


def splitseq(
    indict: dict, 
    num: int =3) -> dict:
    """
    split sequence by step
    """
    names = tuple(indict.keys())
    outdict = {}
    for name in names:
        seq = indict[name]
        Length = len(seq)
        assert Length % 3 == 0
        codes = re.findall("...", seq)
        if codes[-1] in ["TAG", "TGA", "TAA"]:
            codes.remove(codes[-1])
        outdict[name] = codes
    return outdict


def mapseq(
    pdict: dict, 
    udict: dict, 
    cdict: dict,
    ) -> dict:
    """
    map MSA protein sequence to uncleotide
    """
    names = tuple(pdict.keys())
    mapdict = {}
    for name in names:
        i, mseq = 0, ""
        pseq, useq = pdict[name], udict[name]
        for aad in pseq:
            if aad == "-":
                uncl = "-" * 3
            else:
                uncl = useq[i]
                assert uncl in cdict[aad]
                i += 1
            mseq += uncl
        mapdict[name] = mseq
    return mapdict


def translate(
    frame: str, 
    code_dict: dict):
    """
    translate dna to protein
    """
    codes = re.findall('...', frame)
    pros = []
    for c in codes:
        pro = code_dict.get(c, None)
        if pro:
            pros.append(pro)
        elif 'N' in c:
            pros.append('X')
        else:
            # illegal codon
            raise KeyError
    pros = "".join(pros)
    return pros


def reverse_seq(dnaseq):
    """
    reverse the DNA sequence
    """
    basemap = {
            "A": "T",
            "G": "C",
            "C": "G",
            "T": "A",
            "N": "N",
        }
    revseq = [basemap[b] for b in dnaseq[::-1]]
    return "".join(revseq)


def generate_frame(
        dnaseq: str, 
        start_len: int,
        code_dict: dict,
        ) -> tuple:
    """
    generate frame from DNA sequence
    """
    ## think about three offsets
    #start_uncl = dnaseq[:start_len*3+2]
    ## first frame
    #frame1 = start_uncl[0:-2]
    #pros1 = translate(frame1, code_dict)
    ## second frame 
    #frame2 = start_uncl[1:-1]
    #pros2 = translate(frame2, code_dict)
    ## thrid frame 
    #frame3 = start_uncl[2:]
    #pros3 = translate(frame3, code_dict)
    #return (pros1, pros2, pros3)
    # think about three offsets
    for i in range(len(dnaseq)):
        frame = dnaseq[i:start_len*3+i]
        p = translate(frame, code_dict)
        yield i, p


def find_start_position(
    prot_dict: dict,
    uncl_dict: dict,
    code_dict: dict,
    ) -> dict:

    """
    find the location to strat translating
    """

    start_blocks = {}

    # reformat dict
    code_dict = {
        v: k for k, vs in code_dict.items()
        for v in vs
    }

    count = 0
    for name, seq in prot_dict.items():
        # split conception translation sequence
        seq = seq.replace('-', '')
        blocks = re.split("[\\\/]", seq)
        start_seq = blocks[0]
        start_len = len(start_seq)

        uncl = uncl_dict.get(name, None)
        if uncl:
            for pos, transeq in generate_frame(uncl, start_len, code_dict):
                if transeq != start_seq:
                    continue
                start_blocks[name] = pos
                break
            if start_blocks.get(name, None) is None:
                #print(uncl)
                #print(transeq)
                #print(start_seq)
                #print('******')
                # think about minus strand
                uncl = reverse_seq(uncl)
                for pos, transeq in generate_frame(uncl, start_len, code_dict):
                    if transeq != start_seq:
                        continue
                    start_blocks[name] = pos
                    uncl_dict[name] = uncl
                #print(name)
                #print(pros)
                #print(start_seq)
                #print('-'*25)
                #count += 1
    return start_blocks, uncl_dict


def getmap(
    name: str,
    msaseq: str,
    dnaseq: str,
    start: int,
    code_dict: dict,
    base2aad: dict,
    ) -> list:

    """
    mapping protein to dna
    """
    j = 0
    step = 3
    error_count = 0
    uncl_map = ""
    proseq = msaseq.replace('-', '')
    pro_len = len(proseq)
    assert "//" not in proseq
    assert "/\\" not in proseq
    assert "\\/" not in proseq
    assert "\\\\" not in proseq


    for i, aad in enumerate(msaseq):
        if aad == "-":
            code = "-" * 3
        # conceptual translation incolves
        # in the process
        elif aad == "/":
            code = "-" * 3
            step = 2
            j += 1
        # conceptual translation removes
        # a base in the process
        elif aad == "\\":
            code = "-" * 3
            step = 4
            j += 1
        # first residue is specially treated
        elif j == 0:
            code = dnaseq[:start+3]
            if start <= len(uncl_map):
                if start:
                    uncl_map = uncl_map[:-1*start]
            else:
                code = code[start-len(uncl_map):]
                uncl_map = ''
            start = start + 3
            step = 3
            j += 1
        else:
            code = dnaseq[start:start+step]
            start = start + step
            # determine where the base is inserted
            if step == 2:
                aligner = Align.PairwiseAligner()
                for c in code_dict[aad]:
                    alignments = aligner.align(code, c)
                    new_code = str(alignments[0]).split("\n")[0]
                    if new_code.count('-') == 1:
                        break
                if len(new_code) < 3:
                    raise AssertionError
                while len(new_code) > 3:
                    new_code = list(new_code)
                    new_code.remove('-')
                    error_count += 1

                code = new_code
            elif step == 3:
                if base2aad.get(code, None) != aad:
                    if aad != 'X' or 'N' not in code:
                        error_count += 1
            # determine where the base is delected
            elif step == 4:
                dele_base = code[0]
                code = code[1:]
                assert uncl_map[-3:] == "---"
                # make up the deleted bases
                uncl_map = uncl_map[:-1] + dele_base
                if base2aad.get(code, None) != aad:
                    if aad != 'X' or 'N' not in code:
                        error_count += 1
            step = 3
            j += 1

        uncl_map += code
    error_radio = error_count/pro_len
    if error_radio > 0.1:
        uncl_map = None
    return uncl_map



def pseudomap(
    prot_dict: dict,
    uncl_dict: dict,
    code_dict: dict,
    ) -> dict:

    base2aad = {
        v: k for k, vs in code_dict.items()
        for v in vs
    }

    start_indexs, uncl_dict = find_start_position(
        prot_dict,
        uncl_dict,
        code_dict,
    )
    
    result_dist = {}

    for name, start in start_indexs.items():
        msaseq = prot_dict.get(name, None)
        dnaseq = uncl_dict.get(name, None)
        if msaseq and dnaseq:
            uncl_map = getmap(name, msaseq, dnaseq, start, code_dict, base2aad)
            if uncl_map:
                result_dist[name] = uncl_map
    return result_dist


def write2file(
    indict: dict,
    outfile: str) -> None:
    """
    write map sequence to file
    """
    with open(outfile, "w") as out:
        for name in indict.keys():
            seq = ">" + name + "\n" + indict[name]+ "\n"
            out.write(seq)

