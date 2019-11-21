#!/usr/bin/python
import sys
import itertools as it

pop_list = []

#with open("/Users/dashazhernakova/Documents/Doby/GenomeRussia/ancientDNA/GR+Lazaridis.ind") as f:
with open(sys.argv[1]) as f:
	[pop_list.append(l.strip().split("\t")[2]) for l in f if l.strip().split("\t")[2] not in pop_list]
gr=["Adygei","Arkhangelsk","Bashkirs","Komis","Novgorod","Pskov","Rostov","Vladimir","Voronezh","Yakuts","Yaroslavl","Khantys"]
triplets = it.combinations(pop_list, 3)
for a,b,c in triplets:
	print "tmp", a, b, c
	if (c == "Mbuti") and (a in gr):
		print a + "\t" + b + "\t" + c
