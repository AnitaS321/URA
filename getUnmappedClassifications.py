#!/usr/bin/env python3
#Author: Anita S. Bowman (ASB)

import argparse
import sys
import os
import subprocess
import re
import pandas as pd

def main(argv):

    os.chdir(args.outputdir)
    
    bfile = "%s" % (args.file)
    blast_db = "%s" % (args.blastDB)

    if args.count:
        countHHV4Reads(bfile)
    else:
        getUnmappedReads(bfile)
        getBlastHits(blast_db)
        filterBlast_uniqReads()
        runKrona()
        annotate()

def countHHV4Reads(inputbam):
    cmd = "samtools view -c %s NC_007605" % (inputbam)

    if args.verbose:
        print ("Counting EBV reads in %s...\n\n" % (args.sample))
        print(cmd)

    count = subprocess.check_output(cmd, shell=True,universal_newlines=True)
    output="%s_HHVcounts.txt" % (args.sample)
    s = "%s\t%s" % (args.sample, str(count))
    output=open(output,'w')
    output.write(s)

def filterBlast_uniqReads():
    rawBlast = "%s_blast" % (args.sample)
    outBlast = "%s_blast_filt" % (args.sample)
    cmd = "sed -e '/Query/d' -e '/BLASTN/d' -e '/Database/d' -e '/0 hits found/d' -e 's#/[0-9]##g' -e '/#/d' %s > %s" % (rawBlast, outBlast)
    if args.verbose:
         print(cmd)
    subprocess.run(cmd, shell=True)



def annotate():
    cmd = "%s/bin/ktGetTaxInfo -f 2 -tax %s/taxonomy/ < %s_blast.tax > %s_tmp" % (args.krona, args.krona, args.sample, args.sample)
    tmp = "%s_tmp" % (args.sample)
    if args.verbose:
        print("Adding taxonomy name from taxonomy ID mapping file...")
        print(cmd)

    subprocess.run(cmd, shell=True)

    classifications = pd.read_table(tmp,delimiter='\t',index_col=False,header=0,skip_blank_lines=True,skiprows=[1],usecols=["rank","name"])
    original = pd.read_table("%s_blast.tax" % (args.sample),delimiter='\t',index_col=False,header=0,usecols=["#count","taxID"])
    joined = pd.concat([original,classifications],axis=1)
    print(joined)
    joined.to_csv(r"%s_HITS.txt" % (args.sample),index=False,sep='\t')


def runKrona():
    cmd = "%s/bin/ktClassifyBLAST -s %s_blast_filt -o %s_blast.tax" % (args.krona, args.sample, args.sample)

    if args.verbose:
        print("Classifying BLAST hits...\n\n")
        print(cmd)

    subprocess.run(cmd, shell=True)

def getBlastHits(blastdatabase):
    #virus-only DB from NCBI
    cmd = 'blastn -evalue 1e-10 -word_size 28 -db %s -query %s_unmapped.fastq -outfmt 7 -perc_identity 90 -num_threads 4 > %s_blast' % (blastdatabase, args.sample, args.sample)


    if args.verbose:
        print("Running BLAST...\n\n")
        print(cmd)

    subprocess.run(cmd, shell=True)

def getUnmappedReads(inputbam):
    cmd = "samtools fasta -f 0x04 -N %s > %s_unmapped.fastq" % (inputbam, args.sample)

    if args.verbose:
        print ("Gathering unmapped reads from %s...\n\n" % (args.sample))
        print(cmd)

    subprocess.run(cmd, shell=True)

    getnumunmapped = "grep '>' %s_unmapped.fastq|wc -l" % (args.sample)
    out = subprocess.run(getnumunmapped,shell=True,check=True, stdout=subprocess.PIPE, universal_newlines=True)
    bfile=out.stdout
    outfile = "%s_totalUnmapped.txt" % (args.sample)
    fH=open(outfile,'w')
    fH.write("Sample\tRun\tTotal_Unmapped_Reads\n%s\t%s" % (args.sample, bfile))
    fH.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-k", "--krona", type=str,help="Path to KronaTools directory",required=True)
    parser.add_argument("-b", "--blastDB", type=str,help="Path to blast DB",required=True)
    parser.add_argument("-f", "--file", type=str,help="Input bam file (full path)",required=True)
    parser.add_argument("-s", "--sample", type=str, help="Sample Name",required=True)
    parser.add_argument("-o", "--outputdir", type=str, help="output directory",required=True)
    parser.add_argument("-v", "--verbose", action="store_true",help="set verbose")
    parser.add_argument("-c", "--count", action="store_true", help="Count EBV reads")
    args = parser.parse_args()
    main(sys.argv)
