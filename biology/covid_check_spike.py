from lxml import etree
import requests
import os,shutil
from Bio.Align.Applications import ClustalwCommandline


list_merdir = os.listdir(r'E:\covid\covid')
with open('E:/covid19/Wuhan-Hu-1_Spike (1).fasta','r') as fp:
    a = fp.read()
for pdb in list_merdir:
    # if '.ent' in pdb:
    #     pdb = pdb.split('.')[0]
    #     pdb = pdb.split('pdb')[1]
    #     url='https://www.rcsb.org/fasta/entry/%s/display'%(pdb)       
    #     user_agent = {'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.105 Safari/537.36'}
    #     page_text = requests.get(url=url,headers=user_agent).text.split('\n')
    #     for x in range(len(page_text)):
    #         if 'Severe acute respiratory syndrome coronavirus 2' in page_text[x]:
    #             with open('E:/covid/covid/%s.fasta'%(pdb),'w') as fp:
    #                 fp.write(a)
    #                 fp.write('\n\n')
    #                 fp.write(page_text[x])
    #                 fp.write('\n')
    #                 fp.write(page_text[x+1])
    #             print(page_text[x+1])
    #             in_file = 'E:/covid/covid/%s.fasta'%(pdb)
    #             clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file)
    #             os.system(str(clustalw_cline))
        
    #             os.remove('E:/covid/covid/%s.fasta'%(pdb))
    #             os.remove('E:/covid/covid/%s.dnd'%(pdb))
    #     # break
    
    
    # if '.aln' in pdb:
    #     print(pdb)
    #     with open(r'E:\covid\covid\%s'%(pdb),'r') as fp:
    #         fasta = fp.read().split('\n')
    #         # print(fasta[-3],len(fasta))
    #     check_star = ''
    #     for check in range(3,len(fasta)-2):
    #         if check % 4 == 1:
    #             if check == 5:
    #                 continue
                
    #             if fasta[check].split('                           ')[1] == '':
    #                 pdb = pdb.split('.')[0]
    #                 shutil.move(r'E:\covid\covid\%s.aln'%(pdb),r'E:\covid\covid\check')
    #                 shutil.move(r'E:\covid\covid\%s.txt'%(pdb),r'E:\covid\covid\check')
    #                 shutil.move(r'E:\covid\covid\pdb%s.ent'%(pdb),r'E:\covid\covid\check')
    #                 break
    #             check_star+=fasta[check].split('                           ')[1]
    #     # print(check_star)
    #     count = 1
    #     pdb = pdb.split('.')[0]

    #     for x in check_star:
    #         if x != '*':
    #             with open(r'E:\covid\covid\%s.txt'%(pdb),'a') as fp:
    #                 fp.write(str(count))
    #                 fp.write('\n')
    #         count+=1
    #     # break
    
            
        
    if '.txt' in pdb:
        with open(r'E:\covid\covid\%s'%(pdb),'r') as fp:
            first_ = fp.read().split('\n')[0]
            # print(int(first_)+50)
        pdb = pdb.split('.')[0]
        if int(first_)+50 >= 367:
            print(pdb,first_)
            shutil.move(r'E:\covid\covid\pdb%s.ent'%(pdb),r'E:\covid\covid\unknow')
            shutil.move(r'E:\covid\covid\%s.txt'%(pdb),r'E:\covid\covid\unknow')
            shutil.move(r'E:\covid\covid\%s.aln'%(pdb),r'E:\covid\covid\unknow')
        # break
            
    

