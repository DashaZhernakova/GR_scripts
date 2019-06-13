#/usr/bin/python
import sys
from collections import defaultdict
import gzip

def getCallRate(fname, bim_f):
	#cr = defaultdict(list)
	mis_f = open(fname)
	bim = open(bim_f)
	cur_mis = 1
	cur_id = ""
	cur_spl = ["","","","","",""]
	cnt = 0
	mis_f.readline()
	for l in bim.readlines():
		bim_spl = l.strip().split()
		mis_spl= mis_f.readline().strip().split()
		#if (spl[1] != bim_spl[1]):
		#	print spl, bim_spl
		#	sys.exit(1)
		if cur_spl == bim_spl:
			#print "exact match", cur_spl, bim_spl
			cur_mis = 1
                        cur_id = ""
                        cur_spl = ["","","","","",""]
			continue
		if bim_spl[0] + ":" + bim_spl[3] != cur_spl[0] + ":" + cur_spl[3]:
			if cur_id != "": print cur_id
			#print "start"
			#print cur_id, cur_mis, cur_spl
			cur_mis = float(mis_spl[4])
			cur_id = bim_spl[1]
			cur_spl = bim_spl
			#print cur_id, cur_mis, cur_spl
		elif float(mis_spl[4]) <= cur_mis:
			#print "less"
			cur_mis = float(mis_spl[4])
			cur_id = bim_spl[1]
			cur_spl = bim_spl
			#print cur_id, cur_mis, cur_spl
			
getCallRate(sys.argv[1], sys.argv[2])