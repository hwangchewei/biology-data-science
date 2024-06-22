from selenium import webdriver
from selenium.webdriver.common.by import By
from time import sleep
from selenium.webdriver.chrome.options import Options
import os,shutil
import json
import requests
from multiprocessing import Process, Pool
from lxml import etree
from Bio.Align.Applications import ClustalwCommandline
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import csv
from Biopython import biopython

Pcharge = ['HIS','LYS','ARG']
Ncharge = ['ASP','GLU']
Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR']
Hydrophobic = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
aa_codes = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
    'PHE':'F','GLY':'G','HIS':'H','LYS':'K',
    'ILE':'I','LEU':'L','MET':'M','ASN':'N',
    'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TYR':'Y','TRP':'W'} 
aa_code_reverse = dict(zip(aa_codes.values(),aa_codes.keys()))
user_agent = {'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.105 Safari/537.36'}


def json_download(x,download_path):
    chrome_options = Options()
    # chrome_options.add_argument('--headless')
    # chrome_options.add_argument('--disable-gpu')
    prefs = {'download.default_directory' : download_path}
    chrome_options.add_experimental_option('prefs', prefs)
    chrome_options.add_argument('log-level=3')

    bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_options)
    monomer_pdb = x
    print(monomer_pdb)
    check = True
    # if monomer_pdb == '5a5f':
    #     check = True
    if check:
        print('aaa')
        bro.get('https://www.ncbi.nlm.nih.gov/Structure/vastplus/vastplus.cgi?uid=%s'%(monomer_pdb))
        
        try:
            bro.find_element(By.XPATH,'/html/body/div[1]/div/div[3]/div[1]/div[2]/button[2]').click()
            
        except:
            with open(r"%s\no_data.txt"%(download_path),'a') as fp:
                fp.write(monomer_pdb+'\n')
            bro.quit()
            return
        # print(dimer_pdb)
        sleep(2)
        # sleep(10)
    sleep(3)
    bro.quit()
    return


def no_selection(x,download_path):
    path = r'%s\%s'%(download_path,x)
    if not os.path.isdir(path):
        os.mkdir(path)
        shutil.move(r'%s\%s_vastplus.json'%(download_path,x),path)
    if len(x) != 4:
        return
    print(x)
    with open('%s\%s_vastplus.json'%(path,x),'r') as fp:
        a = json.load(fp)
    
    if not os.path.isdir('%s\mer'%(path)):
        os.mkdir('%s\mer'%(path))
    count = 0
    list_dir = os.listdir(download_path)
    with open('./monomer.txt','r') as fp:
        monomer_pdb_list = fp.readlines()
    monomer = []
    for x in monomer_pdb_list:
        monomer.append(x.split()[1])
    with open('./homomer.txt','r') as fp:
        homomer_pdb_list = fp.readlines()
    homomer = []
    for x in homomer_pdb_list:
        homomer.append(x.split()[1])    
    print(list_dir)
    print(monomer)
    # print(homomer)
    for z in a['neighbors']:
        # print(z['pdbId'])
        # print(int(z['rmsd']))
        url = 'https://www.rcsb.org/fasta/entry/%s/display'%(z['pdbId'])
        fasta = requests.get(url=url,headers=user_agent).text
        with open('%s\mer\%s.fasta'%(path,z['pdbId']) , 'w') as fp:
            fp.write(fasta)

def monomer_fasta(x,download_path,rmsd,sequenceIdentity):
    path = r'%s\%s'%(download_path,x)
    if not os.path.isdir(path):
        os.mkdir(path)
        shutil.move(r'%s\%s_vastplus.json'%(download_path,x),path)
    if len(x) != 4:
        return
    print(x)
    with open('%s\%s_vastplus.json'%(path,x),'r') as fp:
        a = json.load(fp)
    
    if not os.path.isdir('%s\mer'%(path)):
        os.mkdir('%s\mer'%(path))
    count = 0
    list_dir = os.listdir(download_path)
    with open('./monomer.txt','r') as fp:
        monomer_pdb_list = fp.readlines()
    monomer = []
    for x in monomer_pdb_list:
        monomer.append(x.split()[1])
    with open('./homomer.txt','r') as fp:
        homomer_pdb_list = fp.readlines()
    homomer = []
    for x in homomer_pdb_list:
        homomer.append(x.split()[1])    
    print(list_dir)
    print(monomer)
    # print(homomer)
    for z in a['neighbors']:
        # print(z['pdbId'])
        # print(int(z['rmsd']))
        if z['pdbId'].lower() not in monomer and z['pdbId'].lower() in homomer and float(z['rmsd']) < 1.7 and float(z['sequenceIdentity']) > 0.2:
            url = 'https://www.rcsb.org/fasta/entry/%s/display'%(z['pdbId'])
            fasta = requests.get(url=url,headers=user_agent).text
            with open('%s\mer\%s.fasta'%(path,z['pdbId']) , 'w') as fp:
                fp.write(fasta)
            count += 1
            
        if z['pdbId'].lower() not in monomer and float(z['rmsd']) < 1.7 and float(z['sequenceIdentity']) > 0.2:
            shutil.move('%s/%s_vastplus.json'%(download_path,z['pdbId']),path)
            
    if count == 0:
        shutil.move(path,r'D:\programming_language\python\monomer_json\no_homomer')
    return
        # print(x.split('_')[0])
        # path = r'D:\programming_language\python\monomer_json\%s'%(x.split('_')[0])

def oligState(x):
    print(x)
    if x == 'oligState' or x == 'no_homomer' or x == 'no':
        return
    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    with open('%s\%s_vastplus.json'%(path,x),'r') as fp:
        a = json.load(fp)
    b = a['query']
    # print(b['oligState'])
    # return
    if b['oligState'] != 'monomeric':
        shutil.move(path,r'D:\programming_language\python\monomer_json\oligState')

def get_fasta(x): 
    # path = r'E:\monomer_json\%s\mer\dimer\really_dimer'%(x)
    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    
    if os.path.isfile('%s\%s.fasta'%(path,x)):
        return
    if x == 'oligState' or x == 'no_homomer' or x == 'no':
        return
    print(x)
    list_dir = os.listdir(path)
    for y in list_dir:
        y = y.split('_')[0]
        url = 'https://www.rcsb.org/fasta/entry/%s/display'%(y)
        fasta = requests.get(url=url,headers=user_agent).text
        with open('%s\%s.fasta'%(path,y) , 'w') as fp:
            fp.write(fasta)
            
def alignedResidues(x):
    print(x)
    path = r'D:\programming_language\python\monomer_json\%s'%(x)   
    
    if len(x) != 4:
        return  
    with open('%s\%s.fasta'%(path,x) , 'r') as fp:
        s = fp.read().split('\n')[1]
    # print(len(s))
    if not os.path.isdir(r'%s\mer\no_aligne'%(path)):
        os.mkdir(r'%s\mer\no_aligne'%(path))
    with open('%s\%s_vastplus.json'%(path,x),'r') as fp:
        a = json.load(fp)
    
    count = 0
    list_merdir = os.listdir(r'D:\programming_language\python\monomer_json\%s\mer'%(x))
    neighbors = []
    neighbors_dict = {}
    for z in a['neighbors']:
        neighbors.append(z['pdbId'])
        neighbors_dict.update({z['pdbId']:z['alignedResidues']})
    # print(neighbors)
    # print(neighbors_dict)
    # break
    for y in list_merdir:
        if y == 'no_aligne' or y == 'monomer' or 'xml' in y or y == 'false_mer':
            continue
        id = y.split('.')[0]
        Residuesaligned = ((neighbors_dict[id.split('.')[0]])/len(s))
        with open('%s\mer\%s.fasta'%(path,id) , 'r') as fp:
            try:
                c = fp.read().split('\n')[1]
            except:
                c = ''
                # shutil.move('%s\mer\%s'%(path,y),r'%s\mer\no_aligne'%(path))
                # os.remove('%s\mer\%s.xml'%(path,id))
                pass
        if c == '':
            shutil.move('%s\mer\%s'%(path,y),r'%s\mer\no_aligne'%(path))
            os.remove('%s\mer\%s.xml'%(path,id))
            continue
        # print(len(s)/len(c))
        print(Residuesaligned)
        if  Residuesaligned < 0.9 or Residuesaligned > 1.1 or len(c)/len(s) > 1.1 or len(c)/len(s) <0.9:
            shutil.move('%s\mer\%s'%(path,y),r'%s\mer\no_aligne'%(path))
            os.remove('%s\mer\%s.xml'%(path,id))
    # print(os.listdir(r'D:\programming_language\python\monomer_json\%s\mer'%(x)))
    if len(os.listdir(r'D:\programming_language\python\monomer_json\%s\mer'%(x))) == 1:
        shutil.move('%s'%(path),r'D:\programming_language\python\monomer_json\no_aligne')
        
def not_monomer(x):
    print(x)
    path = r'D:\programming_language\python\monomer_json\%s'%(x) 
    if x == 'oligState' or x == 'no_homomer' or x == 'no' or x == 'no_aligne' or x == 'monomer':
        return
    if not os.path.isdir(r'%s\mer\monomer'%(path)):
        os.mkdir(r'%s\mer\monomer'%(path))
      
    with open('%s\%s_vastplus.json'%(path,x),'r') as fp:
        a = json.load(fp)
    list_merdir = os.listdir(r'D:\programming_language\python\monomer_json\%s\mer'%(x))
    neighbors = []
    neighbors_dict = {}
    for z in a['neighbors']:
        neighbors.append(z['pdbId'])
        neighbors_dict.update({z['pdbId']:z['numberOfMolecules']})
    for y in list_merdir:
        if y == 'no_aligne' or y == 'monomer':
            continue
        id = y.split('.')[0]
        # print(Residuesaligned)
        if int(neighbors_dict[id]) == 1:
            shutil.move('%s\mer\%s'%(path,y),r'%s\mer\monomer'%(path))
    if len(os.listdir(r'D:\programming_language\python\monomer_json\%s\mer'%(x))) == 2:
        shutil.move('%s'%(path),r'D:\programming_language\python\monomer_json\monomer')
 
def pisa_download(x):
    print(x)
    # path = r'E:\monomer_json\%s\mer\dimer\really_dimer'%(x)
    path = r'D:\programming_language\python\monomer_json\%s'%(x) 
    # if x == 'oligState' or x == 'no_homomer' or x == 'no' or x == 'no_aligne' or x == 'monomer' or x == 'false_mer' or '.' in x:
    #     return
    if len(x) != 4:
        return
    list_merdir = os.listdir(path)   
    for y in list_merdir:
        if y == 'no_aligne' or y == 'monomer':
            continue
        id = y.split('_')[0]
        # print(id)
        url = 'https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?%s'%(id.lower())
        xml_ = requests.get(url=url,headers=user_agent).text
        
        with open('%s\%s.xml'%(path,id) , 'w') as fp:
            fp.write(xml_)
           
def check_mer(x):
    print(x)
    if x == 'oligState' or x == 'no_homomer' or x == 'no' or x == 'no_aligne' or x == 'monomer' or x == 'false_mer' or '.' in x:
        return
    path = r'D:\programming_language\python\monomer_json\%s'%(x) 
    list_merdir = os.listdir(r'D:\programming_language\python\monomer_json\%s\mer'%(x))  
    if not os.path.isdir(r'%s\mer\false_mer'%(path)):
        os.mkdir(r'%s\mer\false_mer'%(path))  
    for y in list_merdir:
        count = 0
        ligand = []
        if '.xml' not in y:
            continue
        tree = etree.parse('%s/mer/%s.xml'%(path,y.split('.')[0]))
        interface_list = tree.xpath('//interface')
        # print(interface_list)
        # print(len(interface_list))
        for interface_number in range(len(interface_list)):
            # print(interface_number)
            svn = tree.xpath('//interface[%d]//int_solv_en/text()'%(interface_number+1))[0]
            # print(svn)
            symop1 = tree.xpath('//interface[%d]//molecule[1]/symop/text()'%(interface_number+1))[0]
            symop2 = tree.xpath('//interface[%d]//molecule[2]/symop/text()'%(interface_number+1))[0]
            m_class1 = tree.xpath('//interface[%d]//molecule[1]/class/text()'%(interface_number+1))[0]
            m_class2 = tree.xpath('//interface[%d]//molecule[2]/class/text()'%(interface_number+1))[0]
            chain_id1 = tree.xpath('//interface[%d]//molecule[1]/chain_id/text()'%(interface_number+1))[0]
            chain_id2 = tree.xpath('//interface[%d]//molecule[2]/chain_id/text()'%(interface_number+1))[0]
        
        # print(symop1,symop2,m_class1,m_class2)
            if (symop1 == 'x,y,z' and symop2.lower() == 'x,y,z') and m_class1 == m_class2 == 'Protein' and float(svn) < -8:
                # print('aaa')
                count = 1
                break
                # count += 1
            if (symop1 == 'x,y,z' and symop2.lower() == 'x,y,z') and (m_class1 == 'Ligand' or m_class2 == 'Ligand'):
                if m_class1 == m_class2:
                    continue
                count = 2
                if m_class1 == 'Ligand' and m_class2 == 'Protein':
                    ligand.append(chain_id1)
                elif m_class2 == 'Ligand' and m_class1 == 'Protein':
                    ligand.append(chain_id2)
                
        
        if ligand != []:
            with open('%s/mer/%s_ligand.txt'%(path,y.split('.')[0]),'w') as fp:
                fp.write(y.split('.')[0]+'\n')
                for z in ligand:
                    fp.write(str(z)+' ')
                fp.write('\n')
            # print(ligand)
        if count == 2:
            print(x,y)
            shutil.move('%s/mer/%s.xml'%(path,y.split('.')[0]),'%s/mer/false_mer'%(path))
            shutil.move('%s/mer/%s.fasta'%(path,y.split('.')[0]),'%s/mer/false_mer'%(path))
            try:
                shutil.move('%s/mer/%s.pdb'%(path,y.split('.')[0]),'%s/mer/false_mer'%(path))
            except:
                pass
            try:
                shutil.move('%s/mer/%s_%s.pdb'%(path,x,y.split('.')[0]),'%s/mer/false_mer'%(path))
            except:
                pass
    for cfile in os.listdir(r'D:\programming_language\python\monomer_json\%s\mer'%(x)):
        if '.xml' in cfile:
            break
    else:
        shutil.move('%s'%(path),r'D:\programming_language\python\monomer_json\false_mer')
          
def protein_type(x):
    print(x)
    if len(x) != 4:
        return

    hydrophilic = 0
    hydrophobic = 0
    
    
    path = r'D:\programming_language\python\monomer_json\%s'%(x) 
    with open('%s\%s.fasta'%(path,x) , 'r') as fp:
        s = fp.read().split('\n')[1]
    
    for y in s:
        try:
            if aa_code_reverse[y] in Hydrophilic:
                hydrophilic += 1
            elif aa_code_reverse[y] in Hydrophobic:
                hydrophobic += 1
        except:
            pass
    i_type = []
    i_list = []
    i_type.append('{:5s}{:^16.3f}{:^16.3f}'.format(x,round(hydrophobic/len(s),3),round(hydrophilic/len(s),3)))
    list_merdir = os.listdir(r'D:\programming_language\python\monomer_json\%s\mer'%(x))  
    for y in list_merdir:
        hydrophilic = 0
        hydrophobic = 0
        if '.fasta' not in y :
            continue
        i_list.append(y.split('.')[0]+' ')
        # print(y)
        # print(y.split('.')[0])
        with open('%s\mer\%s.fasta'%(path,y.split('.')[0]) , 'r') as fp:
            try:
                s = fp.read().split('\n')[1]
            except:
                continue
            for z in s:
                try:
                    if aa_code_reverse[z] in Hydrophilic or aa_code_reverse[z] in Pcharge or aa_code_reverse[z] in Ncharge:
                        hydrophilic += 1
                    elif aa_code_reverse[z] in Hydrophobic:
                        hydrophobic += 1
                except:
                    pass
        i_type.append('{:5s}{:^16.3f}{:^16.3f}'.format(y.split('.')[0],round(hydrophobic/len(s),3),round(hydrophilic/len(s),3)))
        
        # print(mer_i_type)
        # print(residue1[0],buri1[0],energy1[0])
        # break
    with open(r'D:\programming_language\python\monomer_json\protein_type_new.txt' , 'a') as fp:
        fp.write(x)
        fp.write('\n')
        fp.writelines(i_list)
        fp.write('\n')
        fp.write('{:5s}{:^16s}{:^16s}\n'.format('ID','hydrophobic','hydrophilic'))
        for i in i_type:
            fp.write(str(i))
            fp.write('\n')
        fp.write('******************************************************************\n\n')

def interface_check(x):
    print(x)
    if len(x) != 4:
        return
      
    path = r'E:\monomer_json\%s\mer\dimer\really_dimer'%(x)
    print(path)
 
    list_merdir = os.listdir(path)
    
    print(len(list_merdir))
    if len(list_merdir) == 0:
        return
    if not os.path.isdir(r'%s\interface_error'%(path)):
        os.mkdir(r'%s\interface_error'%(path))
    if not os.path.isdir(r'%s\fasta'%(path)):
        os.mkdir(r'%s\fasta'%(path))
    with open('E:\monomer_json\%s\%s.fasta'%(x,x) , 'r') as fp:
        s = fp.read().split('\n')
    
    for y in list_merdir:
        # try:
        print(y)
        if '.fasta' in y:
            print(y)
            # y = y.split('.')[0]+'.fasta'
            # if '5UWN' not in y:
            #     continue
            with open('%s\%s.fasta'%(path,y.split('.')[0]) , 'r') as fp:
                fasta_y = fp.read().split('\n')
            with open(r'%s\fasta\%s_%s.fasta'%(path,x,y.split('.')[0]) , 'w') as fp:
                fp.write(s[0]+'\n')
                fp.write(s[1]+'\n')
                fp.write('\n')
                fp.write(fasta_y[0]+'\n')
                fp.write(fasta_y[1]+'\n')
            clustalw_cline = ClustalwCommandline("clustalw2", infile=r'%s\fasta\%s_%s.fasta'%(path,x,y.split('.')[0]))
            os.system(str(clustalw_cline))
            with open(r'%s\fasta\%s_%s.aln'%(path,x,y.split('.')[0]),'r') as fp:
                ali_seq_list = fp.read().split('\n')[3:]
            count = 0
            mono,mer = '',''
            for ali_seq in ali_seq_list:
                # print(ali_seq)
                if count % 4 == 0:
                    # print(ali_seq)
                    # print(ali_seq.split()[1])
                    mono += ali_seq.split()[1]
                elif count % 4 == 1:
                    mer += ali_seq.split()[1]
                count += 1
            # print(mono)
            # print(mer)
            mono_dic,mer_dic = {},{}
            countm = 1
            for m in range(len(mer)):
                if mer[m] == '-':
                    continue
                mer_dic.update({countm:mer[m]})
                mono_dic.update({countm:mono[m]})
                countm += 1
            # print(mono_dic)
            # print(mer_dic)
            tree = etree.parse('%s/%s.xml'%(path,y.split('.')[0]))
            # print(svn)
            all_interface = tree.xpath('//interface')
            all_right_face = {}
            # print(len(all_interface))
            for z in range(len(all_interface)):
                svn = tree.xpath('//interface[%d]//int_solv_en/text()'%(z+1))[0]
                symop1 = tree.xpath('//interface[%d]//molecule[1]/symop/text()'%(z+1))[0]
                symop2 = tree.xpath('//interface[%d]//molecule[2]/symop/text()'%(z+1))[0]
                m_class1 = tree.xpath('//interface[%d]//molecule[1]/class/text()'%(z+1))[0]
                m_class2 = tree.xpath('//interface[%d]//molecule[2]/class/text()'%(z+1))[0]
                # print(chain_id1,chain_id2)
                if (symop1 == 'x,y,z' and symop2.lower() == 'x,y,z') and m_class1 == m_class2 == 'Protein':
                    
                    # print(z+1)
                    all_right_face.update({float(svn):z+1})
                    # break
            if len(all_right_face) == 0:
                if not os.path.isdir('%s/XYZnocontact'%(path)):
                    os.mkdir('%s/XYZnocontact'%(path))
                all_file = os.listdir(r'%s'%(path))
                for z in all_file:
                    print(y)
                    if y.split('.')[0] in z:
                        shutil.move(r'%s\%s'%(path,z),r'%s\XYZnocontact'%(path))
                continue
            # print(sorted(all_right_face))
            # print(all_right_face[sorted(all_right_face)[0]],sorted(all_right_face)[0])
            interface_number = all_right_face[sorted(all_right_face)[0]]
            residue_no = tree.xpath('//interface[%d]//molecule[1]/residues/residue/seq_num/text()'%(interface_number))
            residue = tree.xpath('//interface[%d]//molecule[1]/residues/residue/name/text()'%(interface_number))
            buri = tree.xpath('//interface[%d]//molecule[1]//residues/residue/bsa/text()'%(interface_number))
            chain_id1 = tree.xpath('//interface[%d]//molecule[1]/chain_id/text()'%(interface_number))[0]
            chain_id2 = tree.xpath('//interface[%d]//molecule[2]/chain_id/text()'%(interface_number))[0]

            mono_list,mer_list,seq_num,special_amino = [],[],[],[]
            fix_no = 1

            for n in range(1,100):
                # print(n)
                try:
                    if aa_codes[residue[0]] == mer_dic[fix_no] and aa_codes[residue[1]] == mer_dic[fix_no+1] and aa_codes[residue[2]] == mer_dic[fix_no+2]:
                        break
                    fix_no += 1    
                except:
                    pass
            for z in range(len(residue)):
                if fix_no+z > len(mer_dic):
                    break
                if int(residue_no[z]) == 0:
                    continue
                # print(residue_no[z],fix_no+z)
                if residue[z] not in aa_codes:
                    # print(residue[z])
                    if float(buri[z]) != 0:
                        special_amino.append('%s %s %s'%(residue[z],fix_no+z,mono_dic[fix_no+z]))
                    continue
                
                if aa_codes[residue[z]] != mer_dic[fix_no+z]:
                    for n in range(1,100):
                        # print(n)
                        if aa_codes[residue[z]] == mer_dic[fix_no+z] and aa_codes[residue[z+1]] == mer_dic[fix_no+z+1] and aa_codes[residue[z+2]] == mer_dic[fix_no+z+2]:
                            break
                    fix_no+=1 
                    # print('aaa',fix_no)
                    # return
                if float(buri[z]) != 0:
                    # print(mer_dic[fix_no+z],residue[z],fix_no+z,z)
                    # print(mer_dic[int(residue_no[z-1])],aa_codes[residue[z]])
                    # print(residue_no[z])
                    try:

                        if aa_codes[residue[z]] == mer_dic[fix_no+z]:
                            mono_list.append(mono_dic[fix_no+z])
                            mer_list.append(mer_dic[fix_no+z])
                            seq_num.append(fix_no+z)
                            # print(mono_dic[int(residue_no[z])],mer_dic[int(residue_no[z])],int(residue_no[z]),z)
                            # print(residue[z])
                    except:
                        special_amino.append('%s %s'%(residue[z],fix_no+z))
                        # mono_list.append(mono_dic[fix_no+z])
                        # mer_list.append(mer_dic[fix_no+z])
                        # seq_num.append(fix_no+z)
            hydrophilic = 0
            hydrophobic = 0
            pcharge,ncharge,total,etotal = 0,0,0,0
            # print(mono_list,mer_list)
            for mo in mono_list:
                # print(mo)
                if mo == '-' or mo.lower() == 'x':
                    etotal+=1
                    continue
                if aa_code_reverse[mo] in Hydrophilic:
                    hydrophilic += 1
                    total += 1
                elif aa_code_reverse[mo] in Hydrophobic:
                    hydrophobic += 1
                    total += 1
                elif aa_code_reverse[mo] in Pcharge:
                    pcharge += 1
                    total += 1
                elif aa_code_reverse[mo] in Ncharge:
                    ncharge += 1
                    total += 1
            # if etotal > total or total == 0:
                # shutil.move('%s/mer/%s.xml'%(path,y.split('.')[0]),r'%s\mer\interface_error'%(path))
                # shutil.move('%s/mer/%s.pdb'%(path,y.split('.')[0]),r'%s\mer\interface_error'%(path))
                # shutil.move('%s/mer/%s.fasta'%(path,y.split('.')[0]),r'%s\mer\interface_error'%(path))
                # elist_merdir = os.listdir(r'D:/programming_language/python/monomer_json/%s/mer'%(x))
                # if '.fasta' not in elist_merdir:
                #     shutil.move(path,r'D:\programming_language\python\monomer_json\interface_error')
                # break
            # print('aaa')
            mobic = round(hydrophobic/total,3)
            mobich = hydrophobic
            molic = round(hydrophilic/total,3)
            molich = hydrophilic
            mopc = round(pcharge/total,3)
            mopch = pcharge
            monc = round(ncharge/total,3)
            monch = ncharge
            hydrophilic = 0
            hydrophobic = 0
            pcharge,ncharge,total = 0,0,0
            for mer_ in mer_list:
                if aa_code_reverse[mer_] in Hydrophilic:
                    hydrophilic += 1
                    total += 1
                elif aa_code_reverse[mer_] in Hydrophobic:
                    hydrophobic += 1
                    total += 1
                elif aa_code_reverse[mer_] in Pcharge:
                    pcharge += 1
                    total += 1
                elif aa_code_reverse[mer_] in Ncharge:
                    ncharge += 1
                    total += 1
                
            merbic = round(hydrophobic/total,3)
            merbich = hydrophobic
            merlic = round(hydrophilic/total,3)
            merlich = hydrophilic
            merpc = round(pcharge/total,3)
            merpch = pcharge
            mernc = round(ncharge/total,3)
            mernch = ncharge
            print(merpch,mernch)
            with open(r'E:\monomer_json\interface_type_test.txt' , 'a') as fp:
                fp.write(y.split('.')[0]+' '+chain_id1+' '+chain_id2+'\n')
                fp.write(mono+'\n')
                fp.write(mer+'\n')
                fp.write(str(mono_list)+'\n')
                fp.write(str(mer_list)+'\n')           
                fp.write(str(seq_num)+'\n')
                fp.write('%s  bic:%d(%.3f)  lic:%d(%.3f)  p:%d(%.3f)  n:%d(%.3f)\n'%(x,mobich,mobic,molich,molic,mopch,mopc,monch,monc))
                fp.write('%s  bic:%d(%.3f)  lic:%d(%.3f)  p:%d(%.3f)  n:%d(%.3f)\n'%(y.split('.')[0],merbich,merbic,merlich,merlic,merpch,merpc,mernch,mernc))
                fp.write(str(special_amino)+'\n')
                fp.write('****************************************\n')
        # except:
        #     # pass
        #     with open(r'D:\programming_language\python\monomer_json\interface_type_error3.txt' , 'a') as fp:
        #         fp.write(x+' '+y.split('.')[0]+'\n')
        # return
        
def dali(x):
    # print(x)
    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    if len(x) != 4:
        return
    list_merdir = os.listdir('%s\mer'%(path))
    for pdb in list_merdir:
        if '.fasta' in pdb:
            # try:
            #     with open('%s\mer\%s.pdb'%(path,pdb.split('.')[0]),'r') as fp:
            #         pdb_line = fp.read().split('\n')
            #     for cline in pdb_line:
            #         if ' B ' in cline:
            #             dali_start(x,pdb)
            
            # try:
            #     if os.path.getsize('%s\mer\%s.pdb'%(path,pdb.split('.')[0])) < os.path.getsize('%s\mer\%s_%s.pdb'%(path,x,pdb.split('.')[0])) or os.path.getsize('%s\mer\%s_%s.pdb'%(path,x,pdb.split('.')[0])) < 100000:
            #         print(os.path.getsize('%s\mer\%s_%s.pdb'%(path,x,pdb.split('.')[0])))
            #         dali_start(x,pdb)
            #     else:
            #         continue
            # except:
                
            dali_start(x,pdb)
            
            # break    
    list_c_merdir = os.listdir('%s\mer'%(path))
    for c in list_c_merdir:
        if '.' in c:
            break
    else:
        shutil.move('%s'%(path),r'D:\programming_language\python\monomer_json\monomer\monomer')            
                    
def dali_start(x,pdb,error_count = 1):
    
    print(x)
    print(pdb)
    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    pdb_id = pdb.split('.')[0]
    if error_count >= 3:
        print('error')
        with open('%s\mer\%s.txt'%(path,pdb_id) , 'w') as fp:
            fp.write('error'+str(error_count))
        return        
    url = 'https://files.rcsb.org/view/%s.pdb'%(pdb_id)
    pdb_file = requests.get(url=url,headers=user_agent).text
    with open('%s\mer\%s.pdb'%(path,pdb_id) , 'w') as fp:
        fp.write(pdb_file)
    for pdb_lines in pdb_file.split('\n'):
        if 'CRYST1' in pdb_lines:
            crystl = pdb_lines
    url = 'http://ekhidna2.biocenter.helsinki.fi/dali/'
    chrome_option = Options()
    chrome_option.add_argument('--headless')
    chrome_option.add_argument('--disable-gpu')
    ua = ' Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/107.0.0.0 Safari/537.36 Edg/107.0.1418.35'
    chrome_option.add_argument("user-agent={}".format(ua))
    bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)
    bro.get('http://ekhidna2.biocenter.helsinki.fi/dali/')
    bro.find_element(By.XPATH,'//*[@id="tabs"]/ul/li[5]').click()
    bro.find_element(By.XPATH,'//*[@id="tabs-2"]/div/form/input[2]').send_keys('%s\mer\%s.pdb'%(path,pdb_id))
    bro.find_element(By.XPATH,'//*[@id="tabs-2"]/div/form/div[4]/div[1]/input[1]').click()
    bro.find_element(By.XPATH,'//*[@id="tabs-2"]/div/form/div[4]/div[2]/input[1]').send_keys('%sA'%(x))
    bro.find_element(By.XPATH,'//*[@id="tabs-2"]/div/form/div[7]/input[2]').click()
    check = True
    count = 0
    with open('%s\mer\%s_%s.pdb'%(path,x,pdb_id),'w') as fp:
        fp.write(crystl+'\n')
        while check:
            a = bro.find_elements(By.XPATH,'/html/body/ul/table/tbody/tr')
            if count == 120:
                print('a error')
                fp.close()
                return dali_start(x,pdb,error_count = error_count+1)
            if a == []:
                sleep(5)
                count += 1
                continue
            # print(a)
            check = False
            # print(pdb_id)
            if len(a) == 1:
                break
            for chain in range(len(a)):
                print(chain)
                try:
                    bro.find_element(By.XPATH,'/html/body/ul/table/tbody/tr[%d]/td[1]/a'%(int(chain)+1)).click()
                    bro.find_element(By.XPATH,'/html/body/form/pre/a[2]').click()
                except:
                    fp.close()
                    print('error')
                    return dali_start(x,pdb,error_count = error_count+1)
                # bro.find_element(By.XPATH,'/html/body/ul/table/tbody/tr[%d]/td[1]/a'%(int(chain)+1))
                # bro.find_element(By.XPATH,'/html/body/form/pre/a[2]')
                text = bro.find_element(By.XPATH,'/html/body/pre').text.split('\n')
                        
                # print(text)
                print(chr(65+chain))
                if chr(65+chain) == 'A':
                    for pdb_text in text:
                        # if 'END' in pdb_text:
                        #     continue
                        if  ' A ' in pdb_text and 'ATOM' in pdb_text:
                            fp.write(pdb_text+'\n')
                else:
                    for pdb_text in text:
                        if  ' A ' in pdb_text and 'ATOM' in pdb_text:
                            # print(pdb_text)
                            # if 'END' in pdb_text:
                            #     continue
                            # if 'ATOM' in pdb_text:
                            pdb_text.replace(' A ',' %s '%(chr(65+chain)))
                            # print(pdb_text.replace(' A ',' %s '%(chr(65+chain)))+'\n')
                            fp.write(pdb_text.replace(' A ',' %s '%(chr(65+chain)))+'\n')
                fp.write('TER\n')
                sleep(0.3)
                bro.back()
                sleep(0.3)
                bro.back()
        fp.write('END\n')    
        
    if len(a) == 1:
        shutil.move('%s/mer/%s.xml'%(path,pdb_id),'%s/mer/monomer'%(path))
        shutil.move('%s/mer/%s.fasta'%(path,pdb_id),'%s/mer/monomer'%(path))
        shutil.move('%s/mer/%s_%s.pdb'%(path,x,pdb_id),'%s/mer/monomer'%(path))
        shutil.move('%s/mer/%s.pdb'%(path,pdb_id),'%s/mer/monomer'%(path))                 

def pisa(x):
    print(x)
    if len(x) != 4:
        return

    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    list_dir = os.listdir(r'%s\mer'%(path))
    for y in list_dir:
        if '.csv' in y:
            return
    # yes_pisa = []
    for y in list_dir:
        if '_' in y and '.xml' in y:
            os.remove(r'%s\mer\%s'%(path,y))
            # print(y)
            # yes_pisa.append(y.split('.')[0].split('_')[1])

    # return
    list_dir = os.listdir(r'%s\mer'%(path))
    if os.path.isfile(r'%s\mer\%s.csv'%(path,x)):
        os.remove(r'%s\mer\%s.csv'%(path,x))
    # return
    for y in list_dir:
        if '_' in y and '.pdb' in y:

            print(y,os.path.getsize(r'%s\mer\%s'%(path,y)))
            if os.path.getsize(r'%s\mer\%s'%(path,y)) < 10000:
                shutil.move(r'%s\mer\%s'%(path,y),r'%s\mer\false_mer'%(path))
                shutil.move(r'%s\mer\%s.xml'%(path,y.split('_')[1].split('.')[0]),r'%s\mer\false_mer'%(path))
                shutil.move(r'%s\mer\%s.pdb'%(path,y.split('_')[1].split('.')[0]),r'%s\mer\false_mer'%(path))
                continue
            
            chrome_option = Options()
            chrome_option.add_argument('--headless')
            chrome_option.add_argument('--disable-gpu')
            chrome_option.add_argument('log-level=3')
            bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)
            bro.get('https://www.ebi.ac.uk/pdbe/pisa/')
            WebDriverWait(bro, 200).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div[2]/div[2]/div/table/tbody/tr/td[2]/table/tbody/tr/td/div/div/div/form/span/span/button'))).click()   #click launch pisa
            WebDriverWait(bro, 1).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[4]/td/u/input'))).click()  #click Coordinate file
            WebDriverWait(bro, 200).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[3]/td/b/input[2]'))).send_keys(r'%s\mer\%s'%(path,y))  #upload file
            WebDriverWait(bro, 200).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[3]/td/b/input[3]'))).click()   #upload
            WebDriverWait(bro, 200).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="pdbeSubmitButton"]/span/button'))).click()
            try:
                sleep(60)
                bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/div/form/input[5]')
                shutil.move(r'%s\mer\%s'%(path,y),r'%s\mer\monomer'%(path))
                shutil.move(r'%s\mer\%s.pdb'%(path,y.split('_')[1].split('.')[0]),r'%s\mer\monomer'%(path))
                shutil.move(r'%s\mer\%s.xml'%(path,y.split('_')[1].split('.')[0]),r'%s\mer\monomer'%(path))
                return pisa(x)
            except:
                pass
            WebDriverWait(bro, 600).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="makePageHeader0"]//tbody/tr[4]/td[2]/span/span/button'))).click()
            overlapping = 'No overlapping'
            overlapping_numer = '0'
            try:
                WebDriverWait(bro, 3).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[3]/td/table/tbody/tr[2]/td/font/font/p[1]/strong')))
                overlapping = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[3]/td/table/tbody/tr[2]/td/font/font/p[1]/strong').text
                print(overlapping)
                if overlapping == 'Overlapping structures':
                    if not os.path.isdir('%s/mer/overlapping'%(path)):
                        os.mkdir(r'%s/mer/overlapping'%(path))
                    shutil.move(r'%s/mer/%s'%(path,y),r'%s/mer/overlapping'%(path))
                    shutil.move(r'%s\mer\%s.xml'%(path,y.split('_')[1].split('.')[0]),r'%s/mer/overlapping'%(path))
                    shutil.move(r'%s\mer\%s.pdb'%(path,y.split('_')[1].split('.')[0]),r'%s/mer/overlapping'%(path))
                    bro.find_element(By.XPATH,'//*[@id="makePageHeader0"]//tbody/tr[2]/td[2]/span/span/button').click()
                    overlapping_numer = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr/td[2]').text
                    bro.find_element(By.XPATH,'//*[@id="makePageHeader0"]//tbody/tr[2]/td[3]/span/span/button').click()
                # continue
            except:
                bro.find_element(By.XPATH,'//*[@id="makePageHeader0"]//tbody/tr[2]/td[2]/span/span/button').click()
            try:
                WebDriverWait(bro, 300).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[1]/td/span[3]/span/button')))
            except:
                return pisa(x)
            with open(r'%s\mer\%s.csv'%(path,x),'a',newline='') as fp:
                writer = csv.writer(fp)
                writer.writerow([y.split('_')[1].split('.')[0],overlapping,overlapping_numer])
            interface_list = bro.find_elements(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr')
            check_face = False
            for interface in range(len(interface_list)):
                
                interface = interface+1
                # print(interface)
                i = bro.find_elements(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td'%(interface))
                p1,p2 = ':',':'
                if len(i) == 20:
                    Nnumber = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[1]'%(interface)).text
                    p1 = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[3]'%(interface)).text
                    p2 = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[8]'%(interface)).text
                    s = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[9]'%(interface)).text
                    g = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[15]'%(interface)).text
                    hb = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[17]'%(interface)).text
                    sb = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[18]'%(interface)).text
                    ds = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[19]'%(interface)).text
                elif len(i) == 21:
                    Nnumber = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[2]'%(interface)).text
                    p1 = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[4]'%(interface)).text
                    p2 = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[9]'%(interface)).text
                    s = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[10]'%(interface)).text
                    g = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[16]'%(interface)).text
                    hb = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[18]'%(interface)).text
                    sb = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[19]'%(interface)).text
                    ds = bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[20]'%(interface)).text
                else:
                    continue
                if (':' not in str(p1) and ':' not in str(p2)) and '  x,y,z  ' == str(s):
                    check_face = True
                    if len(i) == 20:
                        bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[1]'%(interface)).click()
                    elif len(i) == 21:
                        bro.find_element(By.XPATH,'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[%s]/td[2]'%(interface)).click()
                    
                    
                    # bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[1]/td/span[3]/span/button').click()

                    WebDriverWait(bro, 200).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr/td/span[1]/span/button'))).click()
                    # sleep(5)
                    bro.switch_to.window(bro.window_handles[-1])
                    a = bro.find_element(By.CLASS_NAME,'pretty-print').text
                    with open(r'%s\mer\%s_%s%s.xml'%(path,y.split('.')[0],p1,p2),'w') as fp:
                        fp.write(a)
                    with open(r'%s\mer\%s.csv'%(path,x),'a',newline='') as fp:
                        writer = csv.writer(fp)
                        writer.writerow([str(Nnumber),str(p1),str(p2),str(s),str(g),str(hb),str(sb),str(ds)])
                    # sleep(3)
                    bro.close()
                    bro.switch_to.window(bro.window_handles[0])
                    bro.back()
                    # sleep(5)
                else:
                    continue
            bro.quit()
            if check_face == False:
                shutil.move(r'%s\mer\%s'%(path,y),r'%s\mer\monomer'%(path))
                shutil.move(r'%s\mer\%s.pdb'%(path,y.split('_')[1].split('.')[0]),r'%s\mer\monomer'%(path))
                shutil.move(r'%s\mer\%s.xml'%(path,y.split('_')[1].split('.')[0]),r'%s\mer\monomer'%(path))
                
    sleep(5)
    list_dir = os.listdir(r'%s\mer'%(path))
    for z in list_dir:
        if '.pdb' in z:
            break
    else:
        shutil.move(path,r'D:\programming_language\python\monomer_json\overlapping')

def check_fasta_pdb(x):
    if len(x) != 4:
        return

    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    if not os.path.isdir(r'%s\mer\overlapping'%(path)):
        return
    list_dir = os.listdir(r'%s\mer\overlapping'%(path))
    c = False
    
    for y in list_dir:
        if '.fasta' in y:
            c = True
            # print(y.split('.')[0])
            if x+'_'+y.split('.')[0]+'.pdb' not in list_dir:
                print(x,y.split('.')[0])
                shutil.move(r'%s\mer\overlapping\%s.fasta'%(path,y.split('.')[0]),r'%s\mer\alreadymer'%(path))
                # try:
                #     # print(r'%s\mer\%s_%s.pdb'%(path,x,y.split('.')[0]))
                #     shutil.move(r'%s\mer\%s_%s.pdb'%(path,x,y.split('.')[0]),r'%s\mer\monomer'%(path))     
                # except:
                #     pass  
                           
def protein_face(x):
    if len(x) != 4:
        return

    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    list_dir = os.listdir(r'%s\mer'%(path))
    i_list,i_type = [],[]
    for pdb in list_dir:
        if '_' not in pdb and '.pdb' in pdb:
            print(pdb)
            z = biopython('%s\mer\%s'%(path,pdb))
            a = z.surface_residue(4.5)
            hydrophilic = 0
            hydrophobic = 0
            pcharge,ncharge,total = 0,0,0
            # print(mono_list,mer_list)
            for b in a:
                mo = b.get_resname()
                # print(mo)
                if mo in Hydrophilic:
                    hydrophilic += 1
                    total += 1
                elif mo in Hydrophobic:
                    hydrophobic += 1
                    total += 1
                elif mo in Pcharge:
                    pcharge += 1
                    total += 1
                elif mo in Ncharge:
                    ncharge += 1
                    total += 1
            mobic = round(hydrophobic/total,3)
            mobich = hydrophobic
            molic = round(hydrophilic/total,3)
            molich = hydrophilic
            mopc = round(pcharge/total,3)
            mopch = pcharge
            monc = round(ncharge/total,3)
            monch = ncharge
            i_type.append('{:5s}{:^16.10s}{:^16.10s}{:^16.10s}{:^16.10s}'.format(pdb.split('.')[0],r'%s(%s)'%(str(mobic),str(mobich)),r'%s(%s)'%(str(molic),str(molich)),r'%s(%s)'%(str(mopc),str(mopch)),r'%s(%s)'%(str(monc),str(monch))))
            i_list.append(pdb.split('.')[0]+' ')
    with open(r'D:\programming_language\python\monomer_json\protein_face_type.txt' , 'a') as fp:
        fp.write(x)
        fp.write('\n')
        fp.writelines(i_list)
        fp.write('\n')
        fp.write('{:5s}{:^16s}{:^16s}{:^16s}{:^16s}\n'.format('ID','hydrophobic','hydrophilic','positive','negative'))
        for i in i_type:
            print(i)
            fp.write(str(i))
            fp.write('\n')
        fp.write('******************************************************************\n\n')

def pdb_not_in_file(x):
    # print(x)
    if len(x) != 4:
        return

    path = r'D:\programming_language\python\monomer_json\monomer\monomer\%s'%(x)
    try:
        list_dir = os.listdir(r'%s\mer'%(path))
    except:
        print(x)
        shutil.move(path,r'D:\programming_language\python\monomer_json\monomer')
    c = False
    for y in list_dir:
        if '.fasta' in y or '.pdb' in y:
            c = True
            return
    if not c:
        print(x)
        shutil.move(path,r'D:\programming_language\python\monomer_json\monomer')
 
def pisa_mer_check(x):
    print(x)
    if len(x) != 4:
        return

    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    list_dir = os.listdir(r'%s\mer'%(path))
    # for y in list_dir:
    #     if 'alreadymer' in y:
    #         return
    # if os.path.isdir(r'%s\mer\alreadymer'%(path)):
    #     list_dir = os.listdir(r'%s\mer\alreadymer'%(path))
    #     for y in list_dir:
    #         shutil.move(r'%s\mer\alreadymer\%s'%(path,y),r'%s\mer'%(path))
    # if os.path.isfile(r'%s\mer\%s_alreadymer.csv'%(path,x)):
    #     os.remove(r'%s\mer\%s_alreadymer.csv'%(path,x))
    #     print(x)
    # return
    for y in list_dir:
        if '.pdb' in y and '_' in y and 'del' not in y:
        # if '_del.pdb' in y:
            pisa_mer_check_(x,y)
    list_dir = os.listdir(r'%s\mer'%(path))
    for z in list_dir:
        if '.pdb' in z:
            break
    else:
        shutil.move(path,r'D:\programming_language\python\monomer_json\alreadymer')
        # shutil.move(path,r'D:\programming_language\python\monomer_json\fmer')

def pisa_mer_check_(x,y,error_count=0):
    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    chrome_option = Options()
    chrome_option.add_argument('--headless')
    chrome_option.add_argument('--disable-gpu')
    chrome_option.add_argument('log-level=3')
    bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)
    bro.get('https://www.ebi.ac.uk/pdbe/pisa/')
    WebDriverWait(bro, 200).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div[2]/div[2]/div/table/tbody/tr/td[2]/table/tbody/tr/td/div/div/div/form/span/span/button'))).click()   #click launch pisa
    WebDriverWait(bro, 1).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[4]/td/u/input'))).click()  #click Coordinate file
    WebDriverWait(bro, 10).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[3]/td/b/input[2]'))).send_keys(r'%s\mer\%s'%(path,y))  #upload file 
    WebDriverWait(bro, 10).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[3]/td/b/input[3]'))).click()   #upload
    try:
        WebDriverWait(bro, 300).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div[2]/div[2]/div/form/table/tbody//td/span[3]/span/button'))).click()
        WebDriverWait(bro, 200).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div[2]/div[2]/div/form/table[2]/tbody')))
    except:
        if error_count > 3:
            with open(r'D:\programming_language\python\monomer_json\alreadymer_timeout.txt','a') as fp:
                fp.write(x+' '+y+'\n')
            return
        error_count += 1
        print(error_count)
        print(x,y)
        bro.close()
        return pisa_mer_check_(x,y,error_count)
    all_table = bro.find_elements(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr')
    check_mer = False
    if not os.path.isdir(r'%s\mer\alreadymer'%(path)):
        os.mkdir(r'%s\mer\alreadymer'%(path))
    # if not os.path.isdir(r'%s\mer\fmer'%(path)):
    #     os.mkdir(r'%s\mer\fmer'%(path))   
    for t in range(len(all_table)):
        tt = bro.find_elements(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[%s]//tr'%(t))
        # print(len(tt))
        if len(tt) != 0:
            # return
            check_mer = True
            # print(x,y)
            # continue
            with open(r'%s\mer\%s_alreadymer.csv'%(path,x),'a') as fp:
                writer = csv.writer(fp)
                writer.writerow([y.split('_')[1].split('.')[0]])
            for tr in range(3,len(tt)+1):
                tr_xpath = '/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[%s]//tr[%s]'%(t,tr)
                mer_text = bro.find_element(By.XPATH,tr_xpath).text
                print(y,mer_text)
                with open(r'%s\mer\%s_alreadymer.csv'%(path,x),'a') as fp:
                    writer = csv.writer(fp)
                    if tr>3:
                        writer.writerow([' ']+mer_text.split())
                    else:
                        writer.writerow(mer_text.split())
        else:
            continue
    # else:
    #     print(x,y)
    #     all_file = os.listdir(r'%s\mer'%(path))
    #     for z in all_file:
    #         if y.split('_')[0] in z:
    #             shutil.move(r'%s\mer\%s'%(path,z),r'%s\mer\fmer'%(path))
    if check_mer:
        all_file = os.listdir(r'%s\mer'%(path))
        for z in all_file:
            if y.split('.')[0].split('_')[1] in z:
                shutil.move(r'%s\mer\%s'%(path,z),r'%s\mer\alreadymer'%(path))
    bro.close()
            
def pdb_del_ligand(x):
    print(x)
    if len(x) != 4:
        return

    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    list_dir = os.listdir(r'%s\mer'%(path))
    for y in list_dir:
        if '.pdb' in y and '_' not in y:
            with open(r'%s\mer\%s'%(path,y),'r') as fp:
                pdb_file = fp.read().split('\n')
                # for p in pdb_file:
                    # print(p)
            with open(r'%s\mer\%s_del.pdb'%(path,y.split('.')[0]),'w') as fp:
                for p in pdb_file:
                    if 'ATOM' in p or 'CRYST1' in p or 'TER' in p or 'END' in p:
                        fp.write(p+'\n')
                      
def move_file(x):
    print(x)
    if len(x) != 4:
        return
    path = r'D:\programming_language\python\monomer_json\%s'%(x)
    # if os.path.isfile(r'%s\mer\%s_alreadymer.csv'%(path,x)):
        # os.remove(r'%s\mer\%s_alreadymer.csv'%(path,x))
        # shutil.move(path,r'D:\programming_language\python\monomer_json\overlapping')
    #     return
    list_dir = os.listdir(path)
    for y in list_dir:
        # if '.pdb' in y:
        print(x)
        shutil.move(r'%s\mer\alreadymer\%s'%(path,y),r'%s\mer'%(path))
        return
    if os.path.isdir(r'%s\mer\alreadymer'%(path)):
        os.rmdir(r'%s\mer\alreadymer'%(path))
  
def swiss_model(x):
    path = 'E:\monomer_json\%s'%(x)
    list_dir = os.listdir(r'%s\mer'%(path))
    with open('%s/%s.fasta'%(path,x),'r') as fp:
        a = fp.read()
    for pdb in list_dir:
        if '.fasta' in pdb: 
            with open('%s/mer/%s'%(path,pdb),'r') as fp:
                b = fp.read().split('\n')
            with open('%s/model.fasta'%(path),'w') as fp:
                fp.write(a)
                fp.write('\n\n')
                fp.write(b[0])
                fp.write('\n'+b[1])
            in_file = '%s/model.fasta'%(path)
            print(in_file)
            pdb = pdb.split('.')[0]
            clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file)
            os.system(str(clustalw_cline))
            # return 
            ua = ' Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/107.0.0.0 Safari/537.36 Edg/107.0.1418.35'
            chrome_option = Options()
            # chrome_option.add_argument('--headless')
            # chrome_option.add_argument('--disable-gpu')
            chrome_option.add_argument("user-agent={}".format(ua))
            bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)
            bro.get('https://swissmodel.expasy.org/interactive#alignment')
            bro.find_element(By.ID,'id_sequence_file_upload').send_keys('%s/model.aln'%(path))
            WebDriverWait(bro, 200).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div[2]/div[1]/div[2]/form/div[3]/div/div/div[2]/div/button[2]'))).click()

            

            sleep(100)
            check = True
            while check:
                try:
                    bro.find_element(By.XPATH,'//*[@id="mdl_left_col_01"]/div[1]/div[1]/div/div[1]/button')
                    bro.get(bro.current_url+'01.pdb?display=1')
                    print(bro.current_url)
                    check = False
                    text = bro.find_element(By.XPATH,'/html/body/pre').text
                    
                    with open('%s/%s_%smodel.pdb'%(path,x,pdb),'w') as fp:
                        fp.write(text)
                    
                except:
                    print('aaa')
                    sleep(60)

def is_dimer(x,d_):
    # print(x)
    if len(x) != 4:
        return
    path = r'E:\monomer_json\%s\mer'%(x)
    dimer = []
    
    
    list_dir = os.listdir(path+r'\dimer')
    for y in list_dir:
        # print(y.split('.')[0])
        if y.split('.')[0] in d_:
            if not os.path.isdir(r'%s\dimer\really_dimer'%(path)):
                os.mkdir(r'%s\dimer\really_dimer'%(path))
            print(y.split('.'))
            shutil.move(r'%s\dimer\%s'%(path,y),r'%s\dimer\really_dimer'%(path))
    return
    with open(r'%s\%s.csv'%(path,x),newline='') as fp:
        a = csv.reader(fp)
        count = 0 
        for y in reversed(list(a)):
            count += 1
            if len(y) == 3:
                if count == 2:
                    dimer.append(y[0])
                count = 0
            
    # print(dimer)
    try:
        with open(r'%s\%s_alreadymer.csv'%(path,x),newline='') as fp:
            a = csv.reader(fp)
            for y in a:
                if len(y) == 1:
                    dimer.remove(y[0])
    except:
        pass
    # print(dimer)
    if not os.path.isdir(r'%s\dimer'%(path)):
        os.mkdir(r'%s\dimer'%(path))
    list_dir = os.listdir(path)
    for dimer_pdb in dimer:
        for y in list_dir:
            if dimer_pdb in y:
                shutil.move(r'%s\%s'%(path,y),r'%s\dimer'%(path))
        
def interface_change():
    with open(r'E:\monomer_json\interface_type_dimer1.txt','r') as fp:
        a = fp.read().split('=============================================================')
    # with open(r'E:\monomer_json\interface_type_really_dimer.txt','r') as fp:
    #     a = fp.read().split('=============================================================')
    tbic,tlic,tp,tn = 0,0,0,0
    tnbic,tnlic,tnp,tnn = 0,0,0,0 #minus
    tebic,telic,tep,ten = 0,0,0,0 #equl
    gbic,glic,gp,gn = 0,0,0,0
    gnbic,gnlic,gnp,gnn = 0,0,0,0 #minus
    gebic,gelic,gep,gen = 0,0,0,0 #equl
    count,gcount = 0,0
    tsymb,tgsymb = {},{}
    # print(a[2].split('\n')[2])
    pdb_list = []
    for x in a:
        # print(len(x.split('\n')),x.split('\n')[2])
        if len(x.split('\n')) <= 5:
            continue
        pdb_list.append(x.split('\n')[2])
        gsymb = ['0','0','0','0']
        gcount += 1
        bic,lic,p,n = 0,0,0,0
        nbic,nlic,np,nn = 0,0,0,0 #minus
        ebic,elic,ep,en = 0,0,0,0 #equl
        for y in x.split('****************************************'):
            if len(y.split('\n')) < 7:
                continue
            symb = ['0','0','0','0']
            try:
                bic1 = int(y.split('\n')[-4].split('  ')[1].split(':')[1].split('(')[0])
                bic2 = int(y.split('\n')[-3].split('  ')[1].split(':')[1].split('(')[0])
                lic1 = int(y.split('\n')[-4].split('  ')[2].split(':')[1].split('(')[0])
                lic2 = int(y.split('\n')[-3].split('  ')[2].split(':')[1].split('(')[0])
                p1 = int(y.split('\n')[-4].split('  ')[3].split(':')[1].split('(')[0])
                p2 = int(y.split('\n')[-3].split('  ')[3].split(':')[1].split('(')[0])
                n1 = int(y.split('\n')[-4].split('  ')[4].split(':')[1].split('(')[0])
                n2 = int(y.split('\n')[-3].split('  ')[4].split(':')[1].split('(')[0])
                count += 1
                if bic1 < bic2:
                    bic += 1
                    symb[0] = '+'
                elif bic1 > bic2:
                    nbic += 1
                    symb[0] = '-'
                else:
                    ebic += 1
                if lic1 < lic2:
                    lic += 1
                    symb[1] = '+'
                elif lic1 > lic2:
                    nlic += 1
                    symb[1] = '-'
                else:
                    elic += 1
                if p1 < p2:
                    p += 1
                    symb[2] = '+'
                elif p1 > p2:
                    np += 1
                    symb[2] = '-'
                else:
                    ep += 1
                if n1 < n2:
                    n += 1
                    symb[3] = '+'
                elif n1 > n2:
                    nn += 1
                    symb[3] = '-'
                else:
                    en += 1
                symb = str(symb)
                if symb not in tsymb:
                    tsymb.update({symb:1})
                else:
                    tsymb[symb] += 1
            except:
                pass
            
        tbic += bic
        tnbic += nbic
        tebic += ebic
        tlic += lic
        tnlic += nlic
        telic += elic
        tp += p
        tnp += np
        tep += ep
        tn += n
        tnn += nn
        ten += en
        if bic > nbic:
            gbic += 1
            gsymb[0] = '+'
        elif bic < nbic:
            gnbic += 1
            gsymb[0] = '-'
        else:
            gebic += 1
        if lic > nlic:
            glic += 1
            gsymb[1] = '+'
        elif lic < nlic:
            gnlic += 1
            gsymb[1] = '-'
        else:
            gelic += 1
        if p > np:
            gp += 1
            gsymb[2] = '+'
        elif p < np:
            gnp += 1
            gsymb[2] = '-'
        else:
            gep += 1
        if n > nn:
            gn += 1
            gsymb[3] = '+'
        elif n < nn:
            gnn += 1
            gsymb[3] = '-'
        else:
            gen += 1
        gsymb = str(gsymb)
        if gsymb not in tgsymb:
            tgsymb.update({gsymb:1})
        else:
            tgsymb[gsymb] += 1
    print(tbic,tlic,tp,tn)
    print(tnbic,tnlic,tnp,tnn)
    print(tebic,telic,tep,ten)
    print(count)
    print(gbic,glic,gp,gn)
    print(gnbic,gnlic,gnp,gnn)
    print(gebic,gelic,gep,gen)
    print(gcount)
    with open(r'E:\monomer_json\change_type1.csv','w',newline='') as fp:
    # with open(r'E:\monomer_json\really_dimer_change_type.csv','w',newline='') as fp:
        a = csv.writer(fp)
        a.writerow([tbic,tlic,tp,tn])
        a.writerow([-tnbic,-tnlic,-tnp,-tnn])
        a.writerow([tebic,telic,tep,ten])
        a.writerow([count])
        for x in sorted(tsymb.items(), key=lambda x:x[1],reverse=True):
            a.writerow(x)
        a.writerow([gbic,glic,gp,gn])
        a.writerow([-gnbic,-gnlic,-gnp,-gnn])
        a.writerow([gebic,gelic,gep,gen])
        a.writerow([gcount])
        for x in sorted(tgsymb.items(), key=lambda x:x[1],reverse=True):
            a.writerow(x)
    # print(sorted(tsymb.items(), key=lambda x:x[1],reverse=True))
    # print(sorted(tgsymb.items(), key=lambda x:x[1],reverse=True))
    return pdb_list

def interface_amino_change():
    with open(r'E:\monomer_json\interface_type_dimer1.txt','r') as fp:
        a = fp.read().split('=============================================================')
    # with open(r'E:\monomer_json\interface_type_really_dimer.txt','r') as fp:
    #     a = fp.read().split('=============================================================')
    amino_dict,amino_type_dict = {},{}
    same_list = [['same pdbid']]
    bic1,lic1,p1,n1,non1,special_amino1 = 0,0,0,0,0,0
    bic2,lic2,p2,n2,non2,special_amino2 = 0,0,0,0,0,0
    A,C,D,E,F,G,H,K,I,L,M,N,P,Q,R,S,T,V,Y,W = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    A2,C2,D2,E2,F2,G2,H2,K2,I2,L2,M2,N2,P2,Q2,R2,S2,T2,V2,Y2,W2 = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    for x in a:
        if len(x.split('\n')) <= 5:
            continue
        # print(x.split('\n')[2])
        for y in x.split('****************************************'):
            # print(len(y.split('\n')),y.split('\n')[0])
            if len(y.split('\n')) < 7:
                continue
            amino1 = y.split('\n')[-7].replace('[','').replace(']','').replace('\'','').replace(' ','').split(',')
            amino2 = y.split('\n')[-6].replace('[','').replace(']','').replace('\'','').replace(' ','').split(',')
            # return
            # print(amino1,amino2)
            if amino1 == amino2:
                same_list.append([y.split('\n')[-4].split('  ')[0],y.split('\n')[-3].split('  ')[0]])
                if 'same' not in amino_type_dict:
                    amino_type_dict.update({'same':1})
                else:
                    amino_type_dict['same'] += 1
                continue    
            for amino_len in range(len(amino1)):
                if amino1[amino_len] != amino2[amino_len]:
                    amino_change_type = '%s -> %s'%(amino1[amino_len],amino2[amino_len])
                    if amino_change_type not in amino_dict:
                        amino_dict.update({amino_change_type:1})
                    else:
                        amino_dict[amino_change_type] += 1
                    try:
                        a1 = aa_code_reverse[amino1[amino_len]]
                        if amino1[amino_len] == 'A':
                            A += 1
                        elif amino1[amino_len] == 'D':
                            D += 1
                        elif amino1[amino_len] == 'E':
                            E += 1
                        elif amino1[amino_len] == 'F':
                            F += 1
                        elif amino1[amino_len] == 'G':
                            G += 1
                        elif amino1[amino_len] == 'H':
                            H += 1
                        elif amino1[amino_len] == 'K':
                            K += 1
                        elif amino1[amino_len] == 'I':
                            I += 1
                        elif amino1[amino_len] == 'L':
                            L += 1
                        elif amino1[amino_len] == 'M':
                            M += 1
                        elif amino1[amino_len] == 'N':
                            N += 1
                        elif amino1[amino_len] == 'P':
                            P += 1
                        elif amino1[amino_len] == 'Q':
                            Q += 1
                        elif amino1[amino_len] == 'R':
                            R += 1
                        elif amino1[amino_len] == 'S':
                            S += 1
                        elif amino1[amino_len] == 'T':
                            T += 1
                        elif amino1[amino_len] == 'V':
                            V += 1
                        elif amino1[amino_len] == 'Y':
                            Y += 1
                        elif amino1[amino_len] == 'W':
                            W += 1
                        elif amino1[amino_len] == 'C':
                            C += 1
                        if amino2[amino_len] == 'A':
                            A2 += 1
                        elif amino2[amino_len] == 'D':
                            D2 += 1
                        elif amino2[amino_len] == 'E':
                            E2 += 1
                        elif amino2[amino_len] == 'F':
                            F2 += 1
                        elif amino2[amino_len] == 'G':
                            G2 += 1
                        elif amino2[amino_len] == 'H':
                            H2 += 1
                        elif amino2[amino_len] == 'K':
                            K2 += 1
                        elif amino2[amino_len] == 'I':
                            I2 += 1
                        elif amino2[amino_len] == 'L':
                            L2 += 1
                        elif amino2[amino_len] == 'M':
                            M2 += 1
                        elif amino2[amino_len] == 'N':
                            N2 += 1
                        elif amino2[amino_len] == 'P':
                            P2 += 1
                        elif amino2[amino_len] == 'Q':
                            Q2 += 1
                        elif amino2[amino_len] == 'R':
                            R2 += 1
                        elif amino2[amino_len] == 'S':
                            S2 += 1
                        elif amino2[amino_len] == 'T':
                            T2 += 1
                        elif amino2[amino_len] == 'V':
                            V2 += 1
                        elif amino2[amino_len] == 'Y':
                            Y2 += 1
                        elif amino2[amino_len] == 'W':
                            W2 += 1
                        elif amino2[amino_len] == 'C':
                            C2 += 1
                    except:
                        if amino1[amino_len] == '-':
                            a1 = '-'
                            non1 += 1
                        else:
                            a1 = amino1[amino_len]
                            special_amino1 += 1
                    try:    
                        a2 = aa_code_reverse[amino2[amino_len]]
                    except:
                        if amino2[amino_len] == '-':
                            a2 = '-'
                            non2 += 1
                        else:
                            a2 = amino2[amino_len]
                            special_amino2 += 1
                    if a1 in Hydrophobic:
                        a1 = 'Bic'
                        bic1 += 1
                    elif a1 in Hydrophilic:
                        a1 = 'Lic'
                        lic1 += 1
                    elif a1 in Pcharge:
                        a1 = 'P'
                        p1 += 1
                    elif a1 in Ncharge:
                        a1 = 'N'
                        n1 += 1
                    if a2 in Hydrophobic:
                        a2 = 'Bic'
                        bic2 += 1
                    elif a2 in Hydrophilic:
                        a2 = 'Lic'
                        lic2 += 1
                    elif a2 in Pcharge:
                        a2 = 'P'
                        p2 += 1
                    elif a2 in Ncharge:
                        a2 = 'N'
                        n2 += 1
                    else:
                        pass
                    
                    amino_type = '%s -> %s'%(a1,a2)
                    if amino_type not in amino_type_dict:
                        amino_type_dict.update({amino_type:1})
                    else:
                        amino_type_dict[amino_type] += 1

        with open(r'E:\monomer_json\amino_change_type1.csv','w',newline='') as fp:
        # with open(r'E:\monomer_json\really_dimer_amino_change_type.csv','w',newline='') as fp:
            a = csv.writer(fp)
            a.writerows([['A',A,A2],['C',C,C2],['D',D,D2],['E',E,E2],['F',F,F2],['G',G,G2],['H',H,H2],['K',K,K2],['I',I,I2],['L',L,L2],['M',M,M2],['N',N,N2],['P',P,P2],['Q',Q,Q2],['R',R,R2],['S',S,S2],['T',T,T2],['V',V,V2],['Y',Y,Y2],['W',W,W2]])
            for x in sorted(amino_dict.items(), key=lambda x:x[1],reverse=True):
                a.writerow(x)
            for x in sorted(amino_type_dict.items(), key=lambda x:x[1],reverse=True):
                a.writerow(x)
            a.writerows(same_list)
            a.writerow([bic1,lic1,p1,n1,non1,special_amino1])
            a.writerow([bic2,lic2,p2,n2,non2,special_amino2])


def main():
    #下載NCBI Vast
    # pisa_monomer_txt = input('輸入PDBePISA_monomer檔案位置: ')
    # download_path = input('輸入下載位置: ')
    # with open(pisa_monomer_txt,'r') as fp:
    #     for monomer_pdb in fp.readlines():
    #         print(monomer_pdb.split()[1])
    #         json_download(monomer_pdb.split()[1],download_path)
    
    
          
    '''  非必要步驟
    list_dir = os.listdir(r'D:\programming_language\python\monomer_json')
    for pdb_vastplus in list_dir:
        pdb_id = pdb_vastplus.split('_')[0]
        rmsd = input('rmsd: ')
        sequenceIdentity = input('sequenceIdentity: )
        monomer_fasta(pdb_id,download_path,rmsd,sequenceIdentity)
        oligState(pdb_id)
        get_fasta(pdb_id)
        alignedResidues(pdb_id)
        #not_monomer(pdb_id)
    '''     

    # vastplus_path = input('輸入Vast檔案位置: ')
    # vastplus_dir = os.listdir(vastplus_path)
    # for pdb_vastplus in vastplus_dir:
    #     pdb_id = pdb_vastplus.split('_')[0]
    #     if os.path.isdir('%s\%s'%(vastplus_path,pdb_id)):
    #         os.mkdir('%s\%s',%(vastplus_path,pdb_id))
            
    
    '''
    
    '''
    
    
    interface_check('6SVN')
    
    # pdb_list = interface_change()
    # print(pdb_list)
    # interface_amino_change()
    # interface_change()
    # list_dir = os.listdir(r'E:\monomer_json')
    # list_dir = os.listdir(r'D:\programming_language\python\monomer_json\alreadymer')
    # check = False
    # check = True
    # list_ = ['1EF9','1HW6','4XTL']
    # d_ = []
    # with open(r"D:\programming_language\python\dimer.txt",'r') as fp:
    #     a = fp.read().split('\n')
    # for y in a:
    #     # print(y.split('  ')[1])
    #     d_.append(y.split('  ')[1].upper())
    # count = 0
    # path = r'E:\monomer_json'
    # for x in list_dir:
        
    #     try:
    #         dir_ = os.listdir(r'%s\%s\mer\dimer'%(path,x))
    #     except:
    #         continue
    #     for y in dir_:
    #         if '.fasta' in y:
    #             pdb_list.remove(x)
    #             count += 1
    #             break
    #     # if '.fasta' in dir_:
    #     #     pdb_list.remove(x)
    #     #     count += 1
    # print(count,pdb_list)
        # pisa_download('1NP2')
        # get_fasta('1NP2')
        # break
    # swiss_model('1A2A')
    # list_ = ['1AAC', '1AB1', '1AB5', '1ACF', '1ADS', '1V0O', '1VJF', '1W24', '1W5B', '1WAM', '1WDN', '3EY6', '3F3M', '4EI2', '4FEU', '4FH1', '4FYJ', '4G0J', '4GAW', '4GPQ', '4GU5', '4GXP', '4H60', '4HDG', '4HHY', '4IED', '4IOB', '4IRX', '4J8M']
    #     pisa_mer_check(x)
        # move_file(x)
        # pdb_del_ligand(x)
        # pisa_mer_check(x)
        # break
        # if x == '5WWK':
        #     check = True
        # if check:


            
        # protein_type(x)
            # with open(r'E:\monomer_json\interface_type_dimer1.txt' , 'a') as fp:
            #     fp.write(x+'\n')
            # interface_check(x)                                                                                                                                                                                                                                     
            # #     # break  
            # with open(r'E:\monomer_json\interface_type_dimer1.txt' , 'a') as fp:
            #     fp.write('\n'+'============================================================='+'\n'+'\n')
            # break 
    # with open('./monomer.txt','r') as fp:
    #     monomer_pdb_list = fp.readlines()
    # # with open('./monomer_json/no_data.txt') as fp:
    # #     nodata = fp.readlines()
    # pdb_list = []
    # for x in monomer_pdb_list:
        
    #     pdb = x.split()[1].upper()
        
    #     print(pdb)

    # #     a = '%s_vastplus.json'%(pdb)
    # #     b = '%s\n'%(pdb.lower())
    # #     if a not in list_dir and b not in nodata:
    #     json_download(pdb)
    #     break
    #     pdb_list.append(pdb)
    # #     # break
    # # print(len(pdb_list))
    # # print(pdb_list)
    # pool = Pool(6)
    # pool.map(json_download,pdb_list)
    # pool.close()
    # pool.join()
if __name__ == '__main__':
    main()