import sys

empty_out_spl=["1","x","1","1","x","1","1","x","A/T","A","0","x","0","0","0","0","x","0","0","0","0","0"]
with open(sys.argv[1]) as f:
    f.readline()
    print "PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR"
    for l in f:
        spl = l.strip().split()
        empty_spl_cp = empty_out_spl[:]
        empty_spl_cp[0] = spl[8]
        empty_spl_cp[1] = spl[1]
        empty_spl_cp[2] = spl[0]
        empty_spl_cp[3] = spl[2]
        print "\t".join(empty_spl_cp)
