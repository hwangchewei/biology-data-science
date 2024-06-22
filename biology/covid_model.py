from Bio.Align.Applications import ClustalwCommandline
from selenium import webdriver
from lxml import etree
from time import sleep
from selenium.webdriver.common.by import By
from lxml import etree
from selenium.webdriver.chrome.options import Options
import os
from multiprocessing import Process, Pool
from itertools import repeat
 


def covid_model(x,pdb):
    with open('E:/covid19/%s.txt'%(x),'r') as fp:
        a = fp.read()
    with open('E:/covid19/%s.fasta'%(pdb),'r') as fp:
        b = fp.read()
    with open('./%s.fasta'%(x),'w') as fp:
        fp.write(a)
        fp.write('\n\n')
        fp.write(b)
    in_file = './%s.fasta'%(x)
    print(in_file)
    
    clustalw_cline = ClustalwCommandline("clustalw2", infile=in_file)
    os.system(str(clustalw_cline))

    ua = ' Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/107.0.0.0 Safari/537.36 Edg/107.0.1418.35'
    chrome_option = Options()
    # chrome_option.add_argument('--headless')
    # chrome_option.add_argument('--disable-gpu')
    chrome_option.add_argument("user-agent={}".format(ua))
    bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)
    bro.get('https://swissmodel.expasy.org/interactive#alignment')
    bro.find_element(By.ID,'id_sequence_file_upload').send_keys('D:/programming_language/python/%s.aln'%(x))
    sleep(3)
    check = True
    while check:
        try:
            bro.find_element(By.XPATH,'/html/body/div[2]/div[1]/div[2]/form/div[3]/div/div/div[2]/div/button[2]').click()
            check = False
        except:
            pass
    sleep(300)
    check = True
    while check:
        try:
            bro.find_element(By.XPATH,'//*[@id="mdl_left_col_01"]/div[1]/div[1]/div/div[1]/button')
            bro.get(bro.current_url+'01.pdb?display=1')
            print(bro.current_url)
            check = False
            text = bro.find_element(By.XPATH,'/html/body/pre').text
            
            with open('./%s_%smodel.pdb'%(x,pdb),'w') as fp:
                fp.write(text)
            
        except:
            print('aaa')
            sleep(60)
def dali(x,pdb):
    ua = ' Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/107.0.0.0 Safari/537.36 Edg/107.0.1418.35'
    chrome_option = Options()
    chrome_option.add_argument('--headless')
    chrome_option.add_argument('--disable-gpu')
    chrome_option.add_argument("user-agent={}".format(ua))
    bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)
    bro.get('http://ekhidna2.biocenter.helsinki.fi/dali/')
    bro.find_element(By.XPATH,'//*[@id="tabs"]/ul/li[5]').click()
    bro.find_element(By.ID,'taxonomy-1').send_keys('%sA'%(pdb))
    bro.find_element(By.XPATH,'//*[@id="tabs-2"]/div/form/div[4]/div[1]/input[1]').click()
    bro.find_element(By.XPATH,'//*[@id="tabs-2"]/div/form/div[4]/div[2]/input[2]').send_keys('E:/covid19/7v8a/%s_7v8amodel.pdb'%(x))
    bro.find_element(By.XPATH,'//*[@id="tabs-2"]/div/form/div[7]/input[2]').click()
    check = True
    
    while check:
        try:
            bro.find_element(By.XPATH,'/html/body/ul/table/tbody/tr/td[1]/a').click()
            bro.find_element(By.XPATH,'/html/body/form/pre/a[2]').click()
            text = bro.find_element(By.XPATH,'/html/body/pre').text
            check = False
        except:
            sleep(60)
            print('bbb')
    with open('./%sace2.txt'%(pdb),'r') as fp:
        ace2 = fp.read()
    with open('./%s_7v8amodel_%sace2.pdb'%(x,pdb),'w') as fp:
        fp.write(text.split('END')[0])
        fp.write('\n')
        fp.write(ace2)

                
    
        

def main():
    covid_list = ['Alpha','Beta','Delta','Gamma','BA.1','BA.2','BA.2.12.1','BA.2.75','BA.4','BA.5']
    pdb_list = ['7r1a']
    # pdb = input('pdb id: ')
    
    for pdb in pdb_list:
        pool = Pool(len(covid_list))
        pool.starmap(covid_model,zip(covid_list, repeat(pdb)))
        pool.close()
        pool.join()
        sleep(300)
    # for x in covid_list:
    #     dali(x, pdb)
if __name__ == '__main__':
    main()

# from Bio import Phylo
# tree = Phylo.read("covid_dna.dnd", "newick")
# Phylo.draw_ascii(tree)

