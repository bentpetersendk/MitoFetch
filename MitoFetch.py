#!/usr/bin/env python
import string, re
import os, sys

sys.path.insert(0, '.')
import glob, gzip
import pandas as pd
# import numpy as np
# from collections import namedtuple
# from math import log
import shutil
import subprocess
from argparse import ArgumentParser
from Bio.SeqIO.FastaIO import SimpleFastaParser
from operator import itemgetter


def exists(outfile):
    exists_already = os.path.exists(outfile)
    if exists_already:
        print(f"Re-using existing outfile: {outfile}", file=sys.stderr)
        return True
    return False


class MitochondriaFetcher:
    def __init__(self, bbnorm, basedir='.', verbose=True, threads=20, memory=196, samplename='mysample'):
        self.verbose = verbose
        self.threads = threads
        self.memory = memory
        self.bbnorm = bbnorm
        self.basedir = basedir
        self.samplename = samplename
        self.set_base_dir()

    def set_base_dir(self):
        self.print("Working in directory", self.basedir)
        if not os.path.exists(self.basedir):
            os.mkdir(self.basedir)
        os.chdir(self.basedir)

    def mkdirs(self):
        for d in 'reads assemblies ref assemblies/01.spades assemblies/02.initial_best_mito assemblies/03.mtdna'.split():
            if not os.path.exists(d): os.mkdir(d)

    def print(self, *args, sep='', end='\n'):
        print("----------------------> Progress : ", *args, file=sys.stderr, flush=True)

    def run_cmd(self, cmd):
        self.print("Running", cmd)
        status = subprocess.call(cmd, shell=True)
        if status != 0:
            print("Error", cmd, status, file=sys.stderr)
            sys.exit(status)

    def get_cmd_output(self, cmd):
        return subprocess.getoutput(cmd)

    ###method for sample processing below
    def run_on_file_SE(self, read_file):
        self.print("Running on file", read_file)
        self.mkdirs()
        self.read_file = read_file
        self.original_data_file = read_file

    # need to edit below for PE option
    # def run_on_file_PE(self, read_file):
    #    self.print("Running on file", read_file)
    #    self.mkdirs()
    #    self.read_file = read_file
    #    self.original_data_file = read_file

    def trim_SE(self):
        self.print("Trimming")
        outfile = f'reads/{samplename}.trim.fq.gz'

        if not exists(outfile):
            cmd = f"AdapterRemoval --file1 {self.read_file}  --minquality 20  --minlength 30 --trimqualities --basename reads/trim --trimns --threads {self.threads}"
            self.run_cmd(cmd)
            os.rename('reads/trim.truncated', f'reads/{samplename}.trim.fq')
            cmd = f'pigz -p {self.threads} reads/*.fq'
            self.run_cmd(cmd)
        self.read_file = 'reads/{samplename}.trim.fq.gz'

    # need to edit below for PE option
    # def trim_PE(self):
    #    self.print("Trimming")
    #    outfile = f'reads/{samplename}.trim.fq.gz'
    #    if not exists(outfile):
    #        cmd = f"AdapterRemoval --file1 {self.read_file}  --minquality 20  --minlength 30 --trimqualities --basename reads/trim --trimns --threads {self.threads}"
    #        self.run_cmd(cmd)
    #        os.rename('reads/trim.truncated', f'reads/{samplename}.trim.fq')
    #        cmd = f'pigz -p {self.threads} reads/*.fq'
    #        self.run_cmd(cmd)
    #    self.read_file = 'reads/{samplename}.trim.fq.gz'

    def initial_assembly_SE(self):
        self.print("Initial assembly")
        print(sys._getframe().f_code.co_name)
        outfile = 'assemblies/01.spades/scaffolds.fasta'
        if not exists(outfile):
            cmd = f"spades.py -s reads/{samplename}.trim.fq.gz -t {self.threads} -o assemblies/01.spades -m {self.memory}"
            self.run_cmd(cmd)

    # need to edit below for PE option
    # def initial_assembly_PE(self):
    #     self.print("Initial assembly")
    #     print(sys._getframe().f_code.co_name)
    #     outfile = 'assemblies/01.spades/scaffolds.fasta'
    #     if not exists(outfile):
    #         cmd = f"spades.py -s reads/{samplename}.trim.fq.gz -t {self.threads} -o assemblies/01.spades -m {self.memory}"
    #         self.run_cmd(cmd)

    def get_best_initial_mito_hit(self, file='assemblies/01.spades/scaffolds.fasta',
                                  outfile="assemblies/02.initial_best_mito/genome.fa"):
        self.print("Getting best mito hits from the initial assembly", file, outfile)
        print(sys._getframe().f_code.co_name)
        hits = self.blastn_initial_assembly_vs_mitodb(infile=file)

        if len(hits) < 1:
            self.print("Could not find any mitochondria hits in", file)
            sys.exit(0)
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = hits[0]

        self.print("Found %s with alignment length=%s and bitscore=%s" % (qseqid, length, bitscore))
        handle = gzip.open(file, 'rt') if file.endswith('.gz') else open(file)
        entries = list(SimpleFastaParser(handle))

        entry = [x for x in entries if x[0] == qseqid][0]  # what happens here
        if not os.path.exists('assemblies/02.initial_best_mito'): os.mkdir('assemblies/02.initial_best_mito')
        self.print("Writing to", outfile)
        fid = open(outfile, 'w+')
        print('>%s\n%s' % entry, file=fid)
        fid.close()

    def blastn_initial_assembly_vs_mitodb(self,
                                          db="/home/projects/ku-cbd/people/bpet/MitoFetch/databases/all_mitochondria.fasta",
                                          infile='assemblies/01.spades/scaffolds.fasta', min_score=150,
                                          outfile='assemblies/02.initial_best_mito/scaffold.blastn.vs.mitodb'):
        self.print("blastn", db, infile)
        if not exists(outfile):
            cmd = 'blastn -db %s -query %s -out %s -outfmt 6' % (db, infile, outfile)  ## added out
            os.system(cmd)
        res = open(outfile).readlines()
        hits = [x.split('\t') for x in res]
        keep = [x for x in hits if float(x[-1]) >= min_score]
        keep = sorted(keep, key=lambda x: float(x[11]), reverse=True)  # sorting after the last column (bitscore)
        return keep

    def create_bwa_reference(self):
        self.print("Creating BWA ref")
        if not exists('ref/genome.fa'): shutil.copy("./assemblies/02.initial_best_mito/genome.fa", './ref/genome.fa')
        cmd = "bwa index ref/genome.fa"
        self.run_cmd(cmd)

    def map_SE_to_reference(self):
        if exists('ref/genomeDNA.bam'): return
        cmd = f"bwa mem -t {self.threads} -v 1 ref/genome.fa reads/{samplename}.trim.fq.gz | samtools view -Sb - > ref/genomeDNA.bam"
        self.run_cmd(cmd)

    # need to edit below for PE option
    # def map_PE_to_reference(self):
    #    cmd = f"bwa mem -t {self.threads} -v 1 ref/genome.fa reads/{samplename}.trim.fq.gz | samtools view -Sb - > ref/genomeDNA.bam"
    #    self.run_cmd(cmd)
    #    cmd = "samtools sort ref/genomeDNA.bam -o ref/genomeDNA.sort.bam"
    #    self.run_cmd(cmd)

    def sort_bamfile(self):
        if exists('ref/genomeDNA.sort.bam'): return
        cmd = "samtools sort ref/genomeDNA.bam -o ref/genomeDNA.sort.bam"
        self.run_cmd(cmd)

    def coverage_plot(self):
        outfile = "ref/genome.cov"
        if exists(outfile): return
        cmd = "bedtools genomecov -d -ibam ref/genomeDNA.sort.bam > %s" % (outfile)
        self.run_cmd(cmd)

    def get_depth(self, file='ref/genome.cov',
                  mult_factor=5):  # why times 5??????? ... ok cording to the paper it is a balance between nuclear depth and mito depth
        df = pd.read_csv(file, header=None, sep='\t')
        depth = int(df[df.columns[2]].mean() * mult_factor)
        return depth

    def normalise_SE_read_depth_by_kmer(self):
        self.print("Normalising by kmer")
        outfile = f'reads/{samplename}.mtreads.fq'
        if exists(outfile): return
        depth = self.get_depth()
        cmd = f"{self.bbnorm} in=reads/{samplename}.trim.fq.gz passes=1 keepall lowbindepth={depth - 1} highbindepth={depth} outhigh={outfile} -Xmx{self.memory}g"
        self.run_cmd(cmd)

    # need to edit below for PE option
    # def normalise_PE_read_depth_by_kmer(self):
    #    self.print("Normalising by kmer")
    #    outfile = 'reads/mtreads.fq'
    #    if exists(outfile): return
    #    depth = self.get_depth()
    #    cmd = f"{self.bbnorm} in=reads/trim.fq.gz passes=1 keepall lowbindepth={depth-1} highbindepth={depth} outhigh={outfile} -Xmx{self.memory}g"
    #    self.run_cmd(cmd)

    def assemble_mt(self):
        self.print("asembling mito")
        outfile = 'assemblies/03.mtdna/mtdna.fasta'
        if exists(outfile): return
        if not os.path.exists('assemblies/03.mtdna'): os.mkdir('assemblies/03.mtdna')
        cmd = f'pigz -p {self.threads} reads/{samplename}.mtreads.fq'
        self.run_cmd(cmd)
        cmd = f"spades.py -s reads/{samplename}.mtreads.fq.gz -t {self.threads} -o assemblies/03.mtdna -m {self.memory}"
        # cmd = f"megahit -r reads/mtreads.fq --min-count 3 --k-min 21 --k-max 101 --k-step 20 --no-mercy -t {self.threads} --mem-flag 0 -o {outfile}"
        self.run_cmd(cmd)
        # To be deleted shutil.move("assemblies/mtdna/final.contigs.fa","assemblies/mtdna/scaffolds.fa")

    def circularise(self, file='candidates/longest.fna'):
        # not sure this is necessary ... for now ignoring
        return file


if __name__ == '__main__':
    usage = "%prog [options] file (or - for stdin)"
    parser = ArgumentParser(usage)
    parser.add_argument("-v", "--verbose", action="store_true", default=0)
    parser.add_argument("-t", "--threads", action="store", type=int, default=40, help="Threads (default: %(default)s)")
    parser.add_argument("-M", "--memory", action="store", type=int, default=190,
                        help="Memory in GB (default: %(default)s)")
    parser.add_argument("-r", "--read_file", action="store", type=str, required=True)
    parser.add_argument("-rt", "--read_type", action="store", type=str, required=True, choices=['SE', 'PE'],
                        help='Single-End (SE) or Paired-End (PE) sequencing data')
    parser.add_argument("-s", "--samplename", action="store", type=str, default='mysample')
    parser.add_argument("--bbnorm", action="store", type=str,
                        default='/home/projects/ku-cbd/people/bpet/MitoFetch/scripts/bbmap/bbnorm.sh',
                        help="Path to bbnorm.sh (default: %(default)s)")
    parser.add_argument("--trim", action="store_true", default=1, required=True,
                        help="Trim using AdapterRemoval (default: Yes)")
    parser.add_argument("--basedir", action="store", type=str, default='mitfetch_run',
                        help="Name of output directory (default: %(default)s)")
    args = parser.parse_args()
    verbose = args.verbose

    args = parser.parse_args()
    verbose = args.verbose
    read_file = os.path.realpath(args.read_file)
    read_type = args.read.type
    threads = args.threads
    memory = args.memory
    samplename = args.samplename
    bbnorm = args.bbnorm
    basedir = args.basedir

    MF = MitochondriaFetcher(basedir=basedir, verbose=verbose, bbnorm=bbnorm, threads=threads, memory=memory)

    MF.run_on_file_SE(read_file)

    # Add option for PE data
    # MF.run_on_file_PE(read_file)
    MF.trim_SE()

    # Add option for PE data
    # MF.trim_PE()

    MF.initial_assembly_SE()

    # Add option for PE data
    # MF.initial_assembly_PE()

    MF.get_best_initial_mito_hit()
    MF.create_bwa_reference()
    MF.map_SE_to_reference()

    # Add option for PE data
    # MF.map_PE_to_reference()
    # input("Method map_SE_to_reference finished - enter to move on")

    MF.sort_bamfile()
    MF.coverage_plot()
    MF.normalise_SE_read_depth_by_kmer()

    # Add option for PE data
    # MF.normalise_PE_read_depth_by_kmer()
    # input("Method normalise_PE_read_depth_by_kmer finished - enter to move on")
    MF.assemble_mt()

    MF.circularise()
    # we ignore circularising for single end

    ##MF.blastn_mitodb(outfile='assemblies/mtdna/test.out')

