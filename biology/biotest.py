from Biopython import biopython
from Bio.PDB.vectors import calc_angle,calc_dihedral
import numpy


_4twm = biopython('./4twm.pdb')

amino_acid = '---ASKWTYFGPDGENSWS---KKYPSCG-GLLQSPIDLHSDILQYDASLTPLEFQGYNLSANKQFLLTNNGHSVKLNLP-SDMHIQGLQSRYSATQLHLHWGNPNDPHGSEHTVSGQHFAAELHIVHYNSDLYPDASTASNKSEGLAVLAVLIEMGSFNPSYDKIFSHLQHVKYKGQEAFVPGFNIEELLPERTAEYYRYRGSLTTPPCNPTVLWTVFRNPVQISQEQLLALETALYCTHMDDPSPREMINNFRQVQKFDERLVYTSFSQ'
amino_4twm = 'VEDEFSYIDGNPNGPENWGNLKPEWETCGKGMEQSPIQLRDNRVIFDQTLGKLRRNYRAVDAR----LRNSGHDVLVDFKGNAGSLSINRVEYQLKRIHFHSP-------SEHEMNGERFDLEAQLVHESQDQK------------RAVVSILFRFGRADPFLSDLEDFIKQFSNSQKNEINAGVVDPNQLQIDDSAYYRYMGSFTAPPCTEGISWTVMRKVATVSPRQVLLLKQAVNENAINNARPLQPTN-FRSVFYFEQLKSKLGVI-'


count, count_6qnl = 0, 0
error_dict = {}
error_list ,error_list2 = [],[]
aa_codes = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
    'PHE':'F','GLY':'G','HIS':'H','LYS':'K',
    'ILE':'I','LEU':'L','MET':'M','ASN':'N',
    'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TYR':'Y','TRP':'W'}  
keys = aa_codes.keys()
values = aa_codes.values()
aa_codes_c = dict(zip(values,keys))
may_corn = []
for x in range(len(amino_acid)):
    if amino_acid[x] != '-':
        count_6qnl += 1
    if amino_4twm[x] != '-':
        count += 1
    if amino_acid[x] == '-' or amino_4twm[x] == '-':
        continue
    if amino_acid[x] == amino_4twm[x]:
        may_corn.append((aa_codes_c[amino_4twm[x]],count))
        print(count_6qnl)
        continue
def conver(Vp1):
    Vp1 = str(Vp1).split('<')[1].split('>')[0].split('Vector')[1].split(',')
    Vp1 = [float(Vp1[0].split(' ')[1]),float(Vp1[1].split(' ')[1]),float(Vp1[2].split(' ')[1])]
    return Vp1
    
def cal_angle(p1,p2,p3):

    p1 = numpy.array(conver(p1))
    p2 = numpy.array(conver(p2))
    p3 = numpy.array(conver(p3))

    v1 = p1 - p2
    v2 = p3 - p2
    v1_l = numpy.sqrt(v1.dot(v1))
    v2_l = numpy.sqrt(v2.dot(v2))
    return numpy.degrees(numpy.arccos(v1.dot(v2)/(v1_l*v2_l)))

may_corn_list ,not_may_corn = [],[]
for x in range(len(_4twm.chain)-2):
    if x<2:
        continue
    p1 = _4twm.chain[x-2]['N'].get_vector()
    p2 = _4twm.chain[x-1]['N'].get_vector()
    p3 = _4twm.chain[x]['CA'].get_vector()
    p4 = _4twm.chain[x+1]['O'].get_vector()
    p5 = _4twm.chain[x+2]['O'].get_vector()
    degree1 = round(cal_angle(p2,p3,p4),0)
    degree2 = round(cal_angle(p1,p3,p5),0)
    # print(degree1,degree2,str(_4twm.chain[x]).split('=')[2].split(' ')[0])         
    if degree1 < 110 and degree2 < 110:
        # print(str(_4twm.chain[x]).split('=')[2].split(' ')[0])
        for y in may_corn:
            if str(_4twm.chain[x]).split('=')[2].split(' ')[0] == str(y[1]):
                may_corn_list.append(int(str(_4twm.chain[x]).split('=')[2].split(' ')[0]))
                break
        else:not_may_corn.append(str(_4twm.chain[x]).split('=')[2].split(' ')[0])


not_may_corn = []
for x in may_corn:
    if x[1] not in may_corn_list:
        not_may_corn.append(x)
    else:
        pass
print(may_corn_list)
print(not_may_corn)
        

#     if amino_acid[x] in ['D','E','Q','N'] and amino_4twm[x] in ['D','E','Q','N']:
#         count_true += 1
#     elif amino_acid[x] in ['S','T','Y'] and amino_4twm[x] in ['S','T','Y']:
#         count_true += 1
#     elif amino_acid[x] in ['R','K','H'] and amino_4twm[x] in ['R','K','H']:
#         count_true += 1
#     elif amino_acid[x] in ['A','V','I','L','G'] and amino_4twm[x] in ['A','V','I','L','G']:
#         count_true += 1
#     elif amino_acid[x] in ['F','W'] and amino_4twm[x] in ['F','W']:
#         count_true += 1
#     else:
#         error_list2.append([amino_4twm[x],amino_acid[x]])
#         if [amino_4twm[x],amino_acid[x]] not in error_list:           
#             error_list.append([amino_4twm[x],amino_acid[x]])
            
#         if amino_acid[x] not in error_dict:
#             # print('aaa')
#             error_dict.update({amino_acid[x]:[amino_4twm[x]]})
#         else:
#             error_dict[amino_acid[x]].append(amino_4twm[x])
#     # else:
#     #     if [amino_acid[x],amino_4twm[x]] not in error_list:
#     #         error_list.append([amino_acid[x],amino_4twm[x]])
# print(count,count_true)
# print(error_dict)    
# for x in error_list:
#     print(x,error_list2.count(x))
# a_list = ['D','E','Q','N','S','T','Y','R','K','H','A','V','I','L','G','F','W','M','C','P']
# for x in a_list:
#     try:
#         print(x,len(error_dict[x]),amino_acid.count(x),amino_4twm.count(x))
#     except:
#         print(x,0,amino_acid.count(x),amino_4twm.count(x))
# # for x in a_list:
# #     print(x)
# #     if x not in error_dict:
# #         print(x)
# #     if x in error_dict:
# #         print('aa')

# Hydrophobic = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
# Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR','ASP','GLU','HIS','LYS','ARG']
# # Acid = ['ASP','GLU']
# # Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR']
# # Base = ['HIS','LYS','ARG']  
# aa_codes = {
#     'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
#     'PHE':'F','GLY':'G','HIS':'H','LYS':'K',
#     'ILE':'I','LEU':'L','MET':'M','ASN':'N',
#     'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
#     'THR':'T','VAL':'V','TYR':'Y','TRP':'W'}  
# for x in range(len(Hydrophilic)):Hydrophilic[x] = aa_codes[Hydrophilic[x]]
# for x in range(len(Hydrophobic)):Hydrophobic[x] = aa_codes[Hydrophobic[x]]
# # for x in range(len(Base)):Base[x] = aa_codes[Base[x]]
# # for x in range(len(Acid)):Acid[x] = aa_codes[Acid[x]]



    
# # print(error_list)
# count, count_true,count_bictolic,count_lictobic = 0, 0, 0, 0
# for x in range(len(amino_acid)):
#     if amino_acid[x] == '-' or amino_4twm[x] == '-':
#         continue
#     count+=1
#     if amino_acid[x] in Hydrophilic and amino_4twm[x] in Hydrophilic:
#         count_true += 1
#     elif amino_acid[x] in Hydrophobic and amino_4twm[x] in Hydrophobic:
#         count_true += 1
#     elif amino_acid[x] in Hydrophobic and amino_4twm[x] in Hydrophilic:
#         count_bictolic += 1
#     elif amino_acid[x] in Hydrophilic and amino_4twm[x] in Hydrophobic:
#         count_lictobic += 1
#     # elif amino_acid[x] in Base and amino_4twm[x] in Base:
#     #     count_true += 1
#     # elif amino_acid[x] in Acid and amino_4twm[x] in Acid:
#     #     count_true += 1
# print(count,count_true,count_bictolic,count_lictobic)
