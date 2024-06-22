from operator import concat
import requests
from selenium import webdriver
from lxml import etree
from time import sleep
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
import os

# user_agent = {'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.105 Safari/537.36'}

# chrome_option.add_argument("user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/84.0.4147.105 Safari/537.36")
def pisa(key):
    chrome_option = Options()
    chrome_option.add_argument('--headless')
    chrome_option.add_argument('--disable-gpu')
    chrome_option.add_argument("--disable-javascript") 
    bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)

    bro.get('https://www.ebi.ac.uk/pdbe/pisa/')

    bro.find_element(By.XPATH,'//*[@id="pdbeSubmitButton"]/span/button').click()

    pdb_id_input = bro.find_element(By.XPATH,'//*[@id="sform"]/tbody/tr[3]/td/i/u/input[2]')
    pdb_id_input.clear()
    pdb_id_input.send_keys(key)
    sleep(1)
    Analysis = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table/tbody/tr[6]/td/big/strong/input').get_attribute('value')
    ligand = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table/tbody/tr[8]/td').text
    bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table/tbody//button[@name="btn_submit_interfaces"]').click()
    # id = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[2]/th[1]/a/strong')
    # if id.text == ' Id ':
    #     s1 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[3]/td[4]').text
    #     s2 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[3]/td[9]').text
    #     s3 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[4]/td[3]').text
    #     s4 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[4]/td[8]').text
    # else:
    #     s1 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[3]/td[3]').text
    #     s2 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[3]/td[8]').text
    #     s3 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[4]/td[3]').text
    #     s4 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[4]/td[8]').text
        
    s1 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[3]/td[@class="data-1"]').text
    s2 = bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[3]/td[@class="data-left data-2"]').text
    # print(s1,s2)
    bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[1]/td/span[3]/span/button').click()
    try:    
        bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr/td/span[1]/span/button').click()
    except:
        if not os.path.isdir('./pisa/%s/'%(key)):
            os.mkdir('./pisa/%s/'%(key))
            print('aaa')
        with open('./pisa/%s/%s%s%s.xml'%(key,key,s1,s2),'w') as fp:
            fp.write('server error')
        a = 'no result'
        force(a,key,s1,s2,ligand,Analysis)
        return
    bro.switch_to.window(bro.window_handles[-1])
    a = bro.find_element(By.CLASS_NAME,'pretty-print').text
    sleep(1)
    bro.switch_to.window(bro.window_handles[0])
    bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[1]/tbody//tbody/tr[4]/td[2]/span/span/button').click()            
    bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr//td[@class="data-center"]/span[1]/span/button').click()
    bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[4]/td/span[2]/span/button').click()
    b = bro.find_element(By.XPATH,'/html/body/pre')
    if not os.path.isdir('./pisa/%s/'%(key)):
        os.mkdir('./pisa/%s/'%(key))
    with open('./pisa/%s/%s%s%s.xml'%(key,key,s1,s2),'w') as fp:
        fp.writelines(ligand + '\n')
        fp.writelines(Analysis + '\n')
        fp.write(a)
    with open('./pisa/%s/%s.pdb'%(key,key),'w') as fp:
        fp.write(b.text)
    
    # bro.switch_to.window(bro.window_handles[0])
    # if s1 != s3 or s2 != s4:
    #     bro.back()
    #     bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[4]/td[@class="data-0"]/input').click()
    #     bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[1]/td/span[3]/span/button').click()
    #     bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr/td/span[1]/span/button').click()
    #     bro.switch_to.window(bro.window_handles[-1])
    #     a = bro.find_element(By.CLASS_NAME,'pretty-print')
    #     with open('./pisa/%s%s%s.xml'%(key,s3,s4),'w') as fp:
    #         fp.write(a.text)
    #     force(a.text,key,s3,s4)
    # print(s1,s2,s3,s4)
    return force(a,key,s1,s2,ligand,Analysis)

def force(key_xml,key,s1,s2,ligand,Analysis):
    if key_xml == 'no result':
        with open('pisa_force.txt','a') as fp:
            fp.write(key + '\n')
            fp.write(ligand + '\n')
            fp.write(Analysis + '\n')
            fp.write(key_xml + '\n'  + '-----------------------------------------------------------' + '\n')

        return
    tree = etree.XML(key_xml)
    residue = tree.xpath('//RESIDUE2/RESIDUE/STRUCTURE/text()')
    # print(residue)
    buri = tree.xpath('//RESIDUE2/RESIDUE//BURIEDSURFACEAREA/text()')
    energy = tree.xpath('//RESIDUE2/RESIDUE//SOLVATIONENERGY/text()')
    hsdc = tree.xpath('//RESIDUE2/RESIDUE//HSDC/text()')
    # sleep(5)

    # page_text = requests.get(url=url).text
    Hydrophobic = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
    # Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR','ASP','GLU','HIS','LYS','ARG']
    Acid = ['ASP','GLU']
    Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR']
    Base = ['HIS','LYS','ARG']

    # h_bond_dic ,saltbride_dic = {} ,{}
    # for x in range(len(h_bond)):
    #     h_bond[x] = h_bond[x].strip()
    #     h_bond[x] = h_bond[x].split(' ')
    #     h_bond[x] = h_bond[x][1].split('[')
    #     if h_bond[x][0] not in h_bond_dic:
    #         h_bond_dic[h_bond[x][0]] = 1
    #     else:
    #         h_bond_dic[h_bond[x][0]] += 1
    # for x in range(len(saltbride)):
    #     saltbride[x] = saltbride[x].strip()
    #     saltbride[x] = saltbride[x].split(' ')
    #     saltbride[x] = saltbride[x][1].split('[')
    #     if saltbride[x][0] not in saltbride_dic:
    #         saltbride_dic[saltbride[x][0]] = 1
    #     else:
    #         saltbride_dic[saltbride[x][0]] += 1        
        # print(h_bond_dic)
            # h_bond_dic = dict(map(x[0],x.count(x[0])))
        # print(saltbride_dic)
    Hydrophobic_count , Hydrophilic_count , Acid_count , Base_count = 0 ,0 ,0 ,0
    Hydrophobic_energy , Hydrophilic_energy , Acid_energy , Base_energy = 0 ,0 ,0 ,0
    Hydrophobic_hsdc , Hydrophilic_hsdc , Acid_hsdc , Base_hsdc = [] ,[] ,[] ,[]
    for x in range(len(residue)):
        residue[x] = residue[x].strip()
        residue[x] = residue[x].split(' ')
        residue[x][0] = residue[x][0].split(':')
        if float(buri[x]) != 0:
            if residue[x][0][1] in Hydrophobic:
                Hydrophobic_energy -= float(energy[x])
                Hydrophobic_count += 1
                if hsdc[x] != ' ':
                    Hydrophobic_hsdc.append(str(residue[x][1])+' '+str(residue[x][0][1])+' '+str(hsdc[x]).strip())

                    
                # if hsdc[x] != ' ':
                #     Hydrophobic_hsdc.append(hsdc[x])
                # print(residue[x][0][1],residue[x][-1],energy[x])
            elif residue[x][0][1] in Hydrophilic:
                Hydrophilic_energy -= float(energy[x])
                Hydrophilic_count += 1
                # if hsdc[x] != ' ':
                #     Hydrophilic_hsdc.append(hsdc[x])
                if hsdc[x] != ' ':
                    Hydrophilic_hsdc.append(str(residue[x][1])+' '+str(residue[x][0][1])+' '+str(hsdc[x]).strip())
   
            elif residue[x][0][1] in Base:
                Base_energy -= float(energy[x])
                Base_count += 1
                if hsdc[x] != ' ':
                    Base_hsdc.append(str(residue[x][1])+' '+str(residue[x][0][1])+' '+str(hsdc[x]).strip())
                    
            elif residue[x][0][1] in Acid:
                Acid_energy -= float(energy[x])
                Acid_count += 1
                if hsdc[x] != ' ':
                    Acid_hsdc.append(str(residue[x][1])+' '+str(residue[x][0][1])+' '+str(hsdc[x]).strip())
                    
    information = 'Hydrophobic_energy:%f quantity:%d   hsdc:%s\nHydrophilic_energy:%f quantity:%d   hsdc:%s\nBase_energy:%f quantity:%d   hsdc:%s\nAcid_energy:%f quantity:%d   hsdc:%s\n' %(Hydrophobic_energy,Hydrophobic_count,str(Hydrophobic_hsdc),Hydrophilic_energy,Hydrophilic_count,str(Hydrophilic_hsdc),Base_energy,Base_count,str(Base_hsdc),Acid_energy,Acid_count,str(Acid_hsdc))
    with open('pisa_force.txt','a') as fp:
        fp.write(str(key)+str(s1)+str(s2) + '\n')
        fp.write(ligand + '\n')
        fp.write(Analysis + '\n')
        fp.write(information + '\n' + '-----------------------------------------------------------' + '\n')

    
# pisa('6xt8')

with open('dimer_pdb.txt','r') as fp:  
    pdb_list = list(fp.readlines())
count = 0
for pdb in pdb_list:
    pdb.split()[0]
    if 'HYDROLASE' in pdb.split()[-1]:
        count += 1
        print(str(count)+'/8305')
        print(pdb.split()[0])
        pisa(pdb.split()[0])
        sleep(2)
        
    # if count == 10:
    #     break
            

    


