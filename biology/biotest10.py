s1a3na = 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR'
s3rgk = '-GLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDRFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG'
s1a3nb = 'VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR','ASP','GLU','HIS','LYS','ARG']
Hydrophobic = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
aa_codes = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
    'PHE':'F','GLY':'G','HIS':'H','LYS':'K',
    'ILE':'I','LEU':'L','MET':'M','ASN':'N',
    'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TYR':'Y','TRP':'W'}  
keys = aa_codes.keys()
values = aa_codes.values()
aa_codes_c = dict(zip(values,keys))
hydrophilic,hydrophobic = 0,0
for x in s1a3nb:
    if x == '-':
        continue
    elif aa_codes_c[x] in Hydrophilic:
        hydrophilic += 1
    elif aa_codes_c[x] in Hydrophobic:
        hydrophobic += 1
  
  
print(hydrophilic,hydrophobic,hydrophilic/(hydrophobic+hydrophilic),hydrophobic/(hydrophobic+hydrophilic))  