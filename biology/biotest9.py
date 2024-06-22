from Biopython import biopython
import Biopython

amino_acid = '---ASKWTYFGPDGENSWS---KKYPSCG-GLLQSPIDLHSDILQYDASLTPLEFQGYNLSANKQFLLTNNGHSVKLNLP-SDMHIQGLQSRYSATQLHLHWGNPNDPHGSEHTVSGQHFAAELHIVHYNSDLYPDASTASNKSEGLAVLAVLIEMGSFNPSYDKIFSHLQHVKYKGQEAFVPGFNIEELLPERTAEYYRYRGSLTTPPCNPTVLWTVFRNPVQISQEQLLALETALYCTHMDDPSPREMINNFRQVQKFDERLVYTSFSQ'
amino_4twm = 'VEDEFSYIDGNPNGPENWGNLKPEWETCGKGMEQSPIQLRDNRVIFDQTLGKLRRNYRAVDAR----LRNSGHDVLVDFKGNAGSLSINRVEYQLKRIHFHSP-------SEHEMNGERFDLEAQLVHESQDQK------------RAVVSILFRFGRADPFLSDLEDFIKQFSNSQKNEINAGVVDPNQLQIDDSAYYRYMGSFTAPPCTEGISWTVMRKVATVSPRQVLLLKQAVNENAINNARPLQPTN-FRSVFYFEQLKSKLGVI-'


_4twm = biopython('D:/programming_language/python/Biopython/data/1si4_del.pdb')
# surface_r4twm = _4twm.surface_residue(10)
# for x in range(len(surface_r4twm)):
#     surface_r4twm[x] = (str(surface_r4twm[x]).split(' ')[1],int(str(surface_r4twm[x]).split('=')[2].split(' ')[0]))
# # print(surface_r4twm)


# liccount,biccount = 0,0

# for x in surface_r4twm:
#     if x[0] in Biopython.Hydrophilic:
#         liccount += 1
#     elif x[0] in Biopython.Hydrophobic:
#         biccount += 1
# print(liccount,biccount,biccount/(liccount+biccount))

# liccount,biccount = 0,0
# _6qnl = biopython('D:/programming_language/python/Biopython/data/6qnl_del.pdb')
# surface_r6qnl = _6qnl.surface_residue(10)
# for x in range(len(surface_r6qnl)):
#     surface_r6qnl[x] = (str(surface_r6qnl[x]).split(' ')[1],int(str(surface_r6qnl[x]).split('=')[2].split(' ')[0]))
# for x in surface_r6qnl:
#     if x[0] in Biopython.Hydrophilic:
#         liccount += 1
#     elif x[0] in Biopython.Hydrophobic:
#         biccount += 1
# print(liccount,biccount,biccount/(liccount+biccount))

# liccount,biccount = 0,0
# _6eki = biopython('D:/programming_language/python/Biopython/data/4hba_del.pdb')
# surface_r6eki = _6eki.surface_residue(10)
# for x in range(len(surface_r6eki)):
#     surface_r6eki[x] = (str(surface_r6eki[x]).split(' ')[1],int(str(surface_r6eki[x]).split('=')[2].split(' ')[0]))
# for x in surface_r6eki:
#     if x[0] in Biopython.Hydrophilic:
#         liccount += 1
#     elif x[0] in Biopython.Hydrophobic:
#         biccount += 1
# print(liccount,biccount,biccount/(liccount+biccount))
x = 0
while x <= 10:
    _4twm = biopython('D:/programming_language/python/Biopython/data/1si4_del.pdb')
    print(_4twm.surface_area_amount(10))
    x += 1
# print(_6qnl.surface_area_amount(10))
# print(_6eki.surface_area_amount(10))