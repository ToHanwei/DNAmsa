#!coding: utf-8

import argparse


class Argparse():
    def __init__(self):
        self.single_args = ''
        self.batch_args = ''

    def pseudo_parse(self):
        """
        pseudo file parse
        run a MSA only one time
        """
        parser = argparse.ArgumentParser(description="CDS MSA tool (pseudo)")
        parser.add_argument('--msafile', '-m', type=str, metavar='', help="inputfile path, it is a MSA output file (FASTA format)")
        parser.add_argument('--cdsfile', '-c', type=str, metavar='', help="uncleotide sequence file (FASTA formfat)")
        parser.add_argument('--transcode', '-t', type=str, metavar='', help="translate code file path, default='code.csv'", default="./code.csv")
        parser.add_argument('--output', '-o', type=str, metavar='', help="output file path, CDS alignment file", default="output.fas")
        self.pseudo_args = parser.parse_args()

    def single_parse(self):
        """
        single file parse
        run a MSA only one time
        """
        parser = argparse.ArgumentParser(description="CDS MSA tool (sigle)")
        parser.add_argument('--msafile', '-m', type=str, metavar='', help="inputfile path, it is a MSA output file (FASTA format)")
        parser.add_argument('--cdsfile', '-c', type=str, metavar='', help="uncleotide sequence file (FASTA formfat)")
        parser.add_argument('--transcode', '-t', type=str, metavar='', help="translate code file path, default='code.csv'", default="./code.csv")
        parser.add_argument('--output', '-o', type=str, metavar='', help="output file path, CDS alignment file", default="output.fas")
        self.single_args = parser.parse_args()

    def batch_parse(self):
        """
        batch runing
        input folder and batch run ever file
        """
        parser = argparse.ArgumentParser(description="CDS MSA tool (batch)")
        parser.add_argument('--indir', '-i', type=str, metavar='', help="input folder path, there are MSA output files (FASTA format)")
        parser.add_argument('--cdsdir', '-c', type=str, metavar='', help=" save uncleotide sequence files (FASTA formfat), have the same name with input files")
        parser.add_argument('--transcode', '-t', type=str, metavar='', help="translate code file path", default="code.csv")
        parser.add_argument('--outdir', '-o', type=str, metavar='', help="output dir path, save CDS alignment files")
        self.batch_args = parser.parse_args()
