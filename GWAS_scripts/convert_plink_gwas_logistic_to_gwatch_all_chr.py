#!/usr/bin/python
import os
import sys
import argparse

def writeTests(test_list, outbase, module_name):
    tests_out_f = open(outbase + "/T1.csv", "w")
    tests_out_f.write('"col","test"\n')
    tests_dict = {}
    test_prefix = module_name + "/"
    
    for i,test in enumerate(test_list):
        col = str(i + 1)
        tests_dict[test] = col
        tests_out_f.write(col + "," + test_prefix + test + "\n")
    
    tests_out_f.close()
    
    return tests_dict

def writeSnps(bim_fname, outbase, freq_fname, snppos2rs):
    print "Writing SNPs"
    
    print "Getting MAFs"
    maf_dict = {}
    with open(freq_fname) as maf_f:
        maf_dict = dict([tuple(l.strip().split()[1:5:3]) for l in maf_f])
    print "Finished getting MAFs!"
    
    snps_out_f = open(outbase + "/T2.csv", "w")
    snps_out_f.write('"chr","nrow","snp","pos","all","maf"\n')

    cur_chr = ""
    snp_dict = {}
    with open(bim_fname) as bim_f:
        for l in bim_f:
            spl = l.strip().split()
            chr = spl[0]
            
            if cur_chr != chr:
                print "\tchr", chr
                cur_chr = chr
                snp_idx = 1
            snpid = snppos2rs.get(spl[0] + ":" + spl[3],spl[1])
            snps_out_f.write(cur_chr + "," + str(snp_idx) + "," + snpid + "," + spl[3] + "," + spl[5] + "/" + spl[4] + "," + maf_dict[spl[1]] + "\n")
            snp_dict[snpid] = str(snp_idx)
            snp_idx += 1
    snps_out_f.close()
    return snp_dict

def makeRsDict(report_fname):
    chrpos2rs = {}
    with open(report_fname) as f:
        f.readline()
        for l in f:
            spl = l.strip().split()
            chrpos2rs[spl[0] + ":" + spl[2]] = spl[1]
    return chrpos2rs

def allFilesExist(outbase):
    for chr in map(str,range(1,23)):
        t3 = outbase + "/chr" + chr + "_T3.csv"
        t2 = outbase + "/chr" + chr + "_T2.csv"
        if os.path.exists(t3) or os.path.exists(t2):
            return True
    return False

def writeRes(res_fname, test, snp_dict, test_dict, outbase):
    print "Writing tests for", test, "model"
    test_name_translation = {'D' : 'DOM', 'R' : 'REC', 'A' : 'ADD', 'CD' : 'ADD'}
    test_plink = test_name_translation[test]
    fname = outbase + "/T3.csv"
    if os.path.exists(fname):
        res_out_f = open(fname, "a")
    else:
        res_out_f = open(fname, "w")
        res_out_f.write('"chr","nrow","col","pv","qas"\n')
    with open(res_fname) as res_f:
        for l in res_f:
            spl = l.strip().split()
            if spl[4] != test_plink:
                continue
            res_out_f.write(spl[0] + "," + snp_dict[spl[1]] + "," + test_dict[test] + "," + spl[8] + "," + spl[6] + "\n")
    res_out_f.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Prepare tables for GWATCH.')
    parser.add_argument('-d','--dominant', help='Plink GWAS report for dominant model', required=False)
    parser.add_argument('-r','--recessive', help='Plink GWAS report for recessive model', required=False)
    parser.add_argument('-c','--codominant', help='Plink GWAS report for codominant model', required=False)
    parser.add_argument('-a','--allelic', help='Plink GWAS report for allelic model', required=False)
    parser.add_argument('-b','--bim', help='Plink genotype bim file', required=True)
    parser.add_argument('-f','--freq', help='Plink MAF file (--freq --nonfounders)', required=False)
    parser.add_argument('-o','--out', help='Output file name prefix', required=False)
    parser.add_argument('-i','--rsids', action = "store_true", default = False, help='Get rs ids from report files instead of bim/frq file', required=False)
    parser.add_argument('-m','--mod', help='Module name', default = "BOT", required=False)
    args = vars(parser.parse_args())
        
    bim_fname = args["bim"]
    dom_fname = args["dominant"]
    rec_fname = args["recessive"]
    codom_fname = args["codominant"]
    allelic_fname = args["allelic"]
    freq_fname = args["freq"]
    outbase = args["out"]
    rs = args['rsids']
    module_name = args['mod']
    if not outbase:
	outbase = os.getcwd()
    out_f = open(outbase + ".log", "w")
    for arg, val in args.items():
        out_f.write(arg + "\t" + str(val) + "\n")
    out_f.close()
    
    print "NB!!! USE RESULT FILES WITH RS IDS!!!"
    if allFilesExist(outbase):
        print "Tables exist! Please remove them and start again!"
        sys.exit(1)
    if rs:
    	chrpos2rs = makeRsDict(dom_fname)

    snp_dict = writeSnps(bim_fname, outbase, freq_fname, chrpos2rs)
    test_dict = writeTests(['D', 'R', 'A', 'CD'], outbase, module_name)

    writeRes(dom_fname, 'D', snp_dict, test_dict, outbase)
    writeRes(rec_fname, 'R', snp_dict, test_dict, outbase)
    writeRes(allelic_fname, 'A', snp_dict, test_dict, outbase)
    writeRes(codom_fname, 'CD', snp_dict, test_dict, outbase)