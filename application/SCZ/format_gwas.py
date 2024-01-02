#!/~/scratch-midway2/tools/miniconda3/bin/python
"""
format_gwas.py functions to take raw summary statistics data 
in diverse formats and convert them into the same format with
options of QC for SNPs.

"""
import sys, os
import argparse
import numpy as np
import pandas as pd
import gzip
import math


def get_delim(header):
    """
    identify the delimiter in the input txt data
    """
    
    header=header.strip()
    if "," in header:
        delim=","
    elif " " in header:
        delim=" "
    elif "\t" in header:
        delim="\t"
    else:
        sys.stderr.write("Delimiter is something else...")
    return(delim)

def format_header(header, delim):
    """
    output a list of formated column names for the GWAS header
    """
    headerDICT={ "raw_snp": ["snp", "snpid", "marker", "rsID", "rs_id", "rs_dbSNP147"],
                 "chr": ["chr", "chrom", "chromosome", "hg19chrc"],
                 "pos": ["bp", "pos", "position", "SNP_hg19", "base_pair_location", "Position(hg19)"],
                 "A1": ["a1", "A1", "effect_allele", "allele1", "Tested_Allele"],
                 "A2": ["a2", "A2", "non_effect_allele", "allele2",
                       "other_allele"],
                 "beta": ["b", "beta", "effect", "effect_size",
                         "beta_hat"],
                 "OR": ["OR", "odds_ratio", "OR(A1)"],
                 "se": ['se', "StdErr", "standard_error"],
                 "zscore": ["z", "zscore", "z_hat"],
                 "pval": ["p", "pval", "pvalue", "p.value", "P.value", "p-value", "p_value", "P-val"],
                 "freqA1": ["Freq.A1", "maf", "Freq_HapMap", "freqA1", "freq1"],
                 "freqA2": ["Freq.A2", "Freq_HapMap", "freqA2", "freq2"],
                 "N": ["N", "EffN", "totalN", "ngt"],
                 "markerName":["MarkerName", "marker"],
               }
        
    if args.ignore:
        ## ignore columns specified by the user
        ignore_cols=[i.lower() for i in args.ignore.split(',')]
    else:
        ignore_cols=[]    
    
    header=header.rstrip()    
    h_list=[i.lower() for i in header.split(delim)]
    #print(f"The old list is {h_list}")
    if len(ignore_cols)>0:
        for x in ignore_cols:
            h_list=["masked" if i==x else i for i in h_list]
    #print(f"The new list is {h_list}")
    new_headers=h_list

    for key in headerDICT:
        names=[i.lower() for i in headerDICT[key]]
        checknames=[k in names for k in h_list]
        if sum(checknames) > 0:
            #check if any column name exists in the value list that match the key
            n_col=np.where(np.array(checknames)==1)[0][0]
            new_headers[n_col]=key
   
    return(new_headers)


def main():
    print(args)
    out=gzip.open(args.out, 'wt')
    if ".gz" in args.raw_file:
        f=gzip.open(args.raw_file, 'rt')
    else:
        f=open(args.raw_file, 'r')     
    raw_header=f.readline()
    delim=get_delim(raw_header)
    header=format_header(raw_header, delim)
    min_header=["chr", "pos", "beta", "se", "A1", "A2", "raw_snp", "pval"]
    out_header=min_header + ["snp"]
    out.write("%s\n" % " ".join(out_header))
    cols=[]
    #identify what statistics are missing in the header
    diff=set(min_header).difference(set(header))
    
    if "pval" in diff and "se" in diff:
        sys.exit("Missing key statistics")
    if "chr" in diff:
        sys.stderr.write("check if chromosome for SNPs included in chromosome or snp column")
        cols=cols+[np.where(np.array(header)=="pos")[0][0]]
    if "se" not in diff:
        if "beta" in diff:
            min_header[min_header.index("beta")]="OR"
    # add other conditions that have missing se but containing pvals        
    cols=cols+[np.where(np.array(header)==i)[0][0] for i in min_header if i in header]
    #sys.stderr.write("%s\n" % ",".join(min_header))       
        
    for line in f:
        row=line.rstrip().split(delim)
        c_col=cols[min_header.index("chr")]
        bp_col=cols[min_header.index("pos")]
        snp_col=cols[min_header.index("raw_snp")]
        if "beta" not in min_header:
            a=cols[min_header.index("OR")] 
            b=cols[min_header.index("se")]
            c=cols[min_header.index("OR")] #column number for OR
            #check if log(se) has a base e or sth else  
            row[c]=format(math.log(float(row[a])), ".6f")
        if "chr" in row[c_col]:
            row[c_col]=row[c_col].split("chr")[1]
        if c_col==bp_col:
            chrom=row[c_col].split(":")[0]
            pos=row[c_col].split(":")[1]
            out.write("%s\n" % " ".join([chrom, pos]+[row[i] for i in cols[2:]]))
            continue
    
        row=row + ["%s" % ":".join(row[i] for i in [c_col, bp_col, snp_col])]

        #change snp format to rsID/chr:    
        out.write("%s\n" % " ".join([row[i] for i in cols] + [row[-1]]))
    
    f.close()
    out.close()  


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Format GWAS summary \
                                                     statistics data.')
    parser.add_argument('raw_file', type=str, help='a text file for GWAS summary \
                                        association statistics data')
    parser.add_argument('out', help='output file with full path')
    #parser.add_argument('--qc', help='perform QC on SNPs') #need to work on the QC part
    parser.add_argument('--ignore', default=None, type=str, help='comma-separated list of column names to ignore')
    args = parser.parse_args()
    main() 

