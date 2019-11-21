#all_pop = ["AfontovaGora3","French","Anatolia_Neolithic_Boncuklu.SG","Chaudardes1"]
#gr = ["Adygei","Arkhangelsk","Bashkirs"]
#all_pop = ["AfontovaGora3","French","Anatolia_Neolithic_Boncuklu.SG","Chaudardes1","Komis","Arkhangelsk","Pskov","Adygei","Novgorod","Voronezh","Bashkirs","Yakuts","Rostov","Khantys","Yaroslavl","Vladimir","Udmurds","Mordvin","Tatars","Forest-Nenets","Selkups","Kets","Lezgin","North_Ossetian","Circassian","Azerbaijanis","Armenian","Iranian","Kazakhs","Croats","Ukrainians_north","Kyrgyz","Buryats","Mishar-Tatars","Chuvashes","Karelian","Belarusian","Roma","Maris","Evenks","Tabasaran","Kumyks","Chukchi","Avars","Abkhasian","Kabardin","Poles","Kryashen-Tatars","Vepsas","Nganasan","Tundra-Nenets","Evens_Sakha","Koryaks","Balkars","Altaian","Tuvinian","Mansi","Mongolian","German","Finnish","Ingrian","Hungarian","Latvian","Lithuanian","Estonian","Vietnamese_south","Vietnamese_north","Eskimo","Evens_Magadan","Shor","Albanian","Saami","Swedes","Vietnamese_central","Turkmen","Ukrainians_east","Georgian","Ukrainians_west","Motala_HG","LBK_EN","Unetice_EBA","Corded_Ware_Germany","BenzigerodeHeimburg_LN","Bell_Beaker_Germany","Karelia_HG","Halberstadt_LBA","Alberstedt_LN","Samara_Eneolithic","Samara_HG","Poltavka","Esperstedt_MN","Yamnaya_Samara","Srubnaya","Potapovka","Scythian_IA","Srubnaya_Outlier","Karsdorf_LN","Salzmuende_MN","Baalberge_MN","Pazyryk_IA","Early_Sarmatian_IA","Aldy_Bel_IA","Anatolia_Neolithic","Ukraine_Neolithic","Ukraine_Mesolithic","Yamnaya_Ukraine","Ukraine_Eneolithic","Latvia_HG","Latvia_MN","Latvia_LN","Zevakino_Chilikta_IA","Kotias","Kostenki14","Anatolia_Neolithic_Kumtepe.SG","Loschbour_published.DG","MA1_HG.SG","Mbuti","Corded_Ware_Estonia.SG","Unetice_EBA.SG","Nordic_BA.SG","Nordic_LN.SG","Yamnaya_Kalmykia.SG","Nordic_LBA.SG","Sintashta_MBA_RISE.SG","Corded_Ware_Proto_Unetice_Poland.SG","Corded_Ware_Germany.SG","Germany_Bronze_Age.SG","Andronovo.SG","Afanasievo_published.SG","Afanasievo.SG","Russia_EBA.SG","Bell_Beaker_Germany.SG","Bell_Beaker_Czech.SG","Nordic_MN_B.SG","BattleAxe_Sweden.SG","Aleut","English","Eskimo_Chaplin","Eskimo_Naukan","Eskimo_Sireniki","Even","Icelandic","Japanese","Norwegian","Spanish","Turkish","Satsurblia","Anatolia_Neolithic_Tepecik_Ciftlik.SG","Ust_Ishim_HG_published.DG"]
#gr = ["Adygei","Arkhangelsk","Bashkirs","Komis","Novgorod","Pskov","Rostov","Vladimir","Voronezh","Yakuts","Yaroslavl","Khantys"]
#gr = ["Pskov"]
#anc = ["AfontovaGora3","Anatolia_Neolithic_Boncuklu.SG","Chaudardes1","Motala_HG","LBK_EN","Unetice_EBA","Corded_Ware_Germany","BenzigerodeHeimburg_LN","Bell_Beaker_Germany","Karelia_HG","Halberstadt_LBA","Alberstedt_LN","Samara_Eneolithic","Samara_HG","Poltavka","Esperstedt_MN","Yamnaya_Samara","Srubnaya","Potapovka","Scythian_IA","Srubnaya_Outlier","Karsdorf_LN","Salzmuende_MN","Baalberge_MN","Pazyryk_IA","Early_Sarmatian_IA","Aldy_Bel_IA","Anatolia_Neolithic","Ukraine_Neolithic","Ukraine_Mesolithic","Yamnaya_Ukraine","Ukraine_Eneolithic","Latvia_HG","Latvia_MN","Latvia_LN","Zevakino_Chilikta_IA","Kotias","Kostenki14","Anatolia_Neolithic_Kumtepe.SG","Loschbour_published.DG","MA1_HG.SG","Mbuti","Corded_Ware_Estonia.SG","Unetice_EBA.SG","Nordic_BA.SG","Nordic_LN.SG","Yamnaya_Kalmykia.SG","Nordic_LBA.SG","Sintashta_MBA_RISE.SG","Corded_Ware_Proto_Unetice_Poland.SG","Corded_Ware_Germany.SG","Germany_Bronze_Age.SG","Andronovo.SG","Afanasievo_published.SG","Afanasievo.SG","Russia_EBA.SG","Bell_Beaker_Germany.SG","Bell_Beaker_Czech.SG","Nordic_MN_B.SG","BattleAxe_Sweden.SG","Satsurblia","Anatolia_Neolithic_Tepecik_Ciftlik.SG","Ust_Ishim_HG_published.DG"]
all_pop = ["Estonians","Finnish","Ingrians","Karelians","Khantys","Komis","Mansis","Maris","Saami","Seto","Udmurts","Vepsas","Hungarians","Forest-Nenets","Nganasans","Selkups","Tundra-Nenets","Arkhangelsk","Novgorod","Pskov","Rostov","Rus_Udmurtia","Vladimir","Voronezh","Yaroslavl"]
gr = ["Estonians","Finnish","Ingrians","Karelians","Khantys","Komis","Mansis","Maris","Saami","Seto","Udmurts","Vepsas","Hungarians","Forest-Nenets","Nganasans","Selkups","Tundra-Nenets"]

out_path = "//Users/dashazhernakova/Documents/Doby/GenomeRussia/freeze2019/population_genetics/Dstats/quadros/"
out = open(out_path + "all.quadros_gr_dif.txt", "w")
written = set()
cnt = 0
file_n = 1
for pop in all_pop: 
    written = set()
    for pop2 in gr:
        for pop3 in gr:
            if (not pop2 == pop3) and (not pop2 == pop) and (not pop3 == pop):
                    if not (pop2, pop3) in written:
                        out.write("Mbuti\t" + pop + "\t" + pop2 + "\t" + pop3 + "\n")
                        #print "Mbuti\t" + pop + "\t" + pop2 + "\t" + pop3
                        written.add((pop2, pop3))
                        #cnt += 1
                        #if cnt > 75000:
                        #    
                        #    out.close()
                        #    file_n += 1
                        #    cnt = 0
                        #    print "next", out_path + "all.quadros_" + str(file_n) + ".txt"
                        #    out = open(out_path + "all.quadros_" + str(file_n) + ".txt", "w")
out.close()