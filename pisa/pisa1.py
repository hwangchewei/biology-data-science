from selenium import webdriver
from time import sleep
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from lxml import etree
import requests
import os

def pisa(pdb_id,protein_type):
    if not os.path.isdir('./pisa1/%s/'%(pdb_id)):
        os.mkdir('./pisa1/%s/'%(pdb_id))
        page_txt = requests.get(url='https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?%s'%(pdb_id)).text
        with open('./pisa1/%s/%s.xml'%(pdb_id,pdb_id),'w') as fp:
            fp.write(page_txt)
        sleep(3)
        print('aaa')
    chrome_option = Options()
    chrome_option.add_argument('--headless')
    chrome_option.add_argument('--disable-gpu')
    chrome_option.add_argument('--disable-javascript')
    chrome_option.add_argument('blink-settings=imagesEnabled=false')

    bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)

    bro.get('http://www.ebi.ac.uk/pdbe/pisa/cgi-bin/piserver?qa=%s'%(pdb_id))
    bro.find_element(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[2]/tbody//button[@name="btn_assembly_details"]').click()
    interface = bro.find_elements(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[3]/tbody/tr[2]/td/table/tbody/tr')
    count = 0
    atom_list ,label_list= [],[]
    for x in interface:
        count += 1
        if '×' in str(x.text):
            atom = bro.find_elements(By.XPATH,'/html/body/div[2]/div[2]/div/form/table[3]/tbody/tr[2]/td/table/tbody/tr[%d]/td/a'%(count))
            for y in atom:
                atom_list.append(y.text.split('+'))
                
    
    
    tree = etree.parse('./pisa1/%s/%s.xml'%(pdb_id,pdb_id))
    x = tree.xpath('//interface/id/text()')
    count = 0
    protein_energy = ['structure DB read fault']
    for atom in atom_list: 
        count += 1
        for y in range(len(x)):            
            m1 = tree.xpath('//interface[%d]//molecule[1]/chain_id/text()'%(y+1))
            m2 = tree.xpath('//interface[%d]//molecule[2]/chain_id/text()'%(y+1))
            css = tree.xpath('//interface[%d]/css/text()'%(y+1))
            # print(m1[0] == atom[0].strip(),m2[0] == atom[1].strip(),symop1[0] == 'x,y,z',symop2[0] == 'X,Y,Z')
            # print(m1[0],atom[0],m2[0],atom[1])
            if m1[0] == atom[0].strip() and m2[0] == atom[1].strip() and float(css[0]) != 0:
                # print(css)
                if count == 1:
                    protein_energy = tree.xpath('//interface[%d]/int_solv_en/text()'%(y+1))
                label_list.append((y+1))
                break
    print(len(atom_list),len(label_list),atom_list)
    with open('pisa2_force.txt','a') as fp:
        fp.write('-----------------------------------------------------------\n')
        fp.write(pdb_id+'@'+protein_type+'\n')
    with open('pisa2_summary.txt','a') as fp:
        fp.write('-----------------------------------------------------------\n')
        fp.write(pdb_id+'@'+protein_type+'$'+str(label_list)+'\n'+str(len(atom_list))+ ' ' +str(len(label_list))+'\n'+atom_list[0][0]+atom_list[0][1]+'\n'+protein_energy[0]+'\n')
    bro.quit()
    return force(label_list,pdb_id)
def force(all_label,pdb_id):
    for label in all_label:
        tree = etree.parse('./pisa1/%s/%s.xml'%(pdb_id,pdb_id))
        Hydrophobic = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
        Acid = ['ASP','GLU']
        Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR']
        Base = ['HIS','LYS','ARG']
        for y in range(1,3):
            chain_id = tree.xpath('//interface[%d]//molecule[%d]/chain_id/text()'%(label,y))
            residue = tree.xpath('//interface[%d]//molecule[%d]//residues/residue/name/text()'%(label,y))
            buri = tree.xpath('//interface[%d]//molecule[%d]//residues/residue/bsa/text()'%(label,y))
            energy = tree.xpath('//interface[%d]//molecule[%d]//residues/residue/solv_en/text()'%(label,y))
            ser_no = tree.xpath('//interface[%d]//molecule[%d]//residues/residue/ser_no/text()'%(label,y))
            Hydrophobic_count , Hydrophilic_count , Acid_count , Base_count = 0 ,0 ,0 ,0
            Hydrophobic_energy , Hydrophilic_energy , Acid_energy , Base_energy , else_energy = 0 ,0 ,0 ,0 ,0
            Hydrophobic_hsdc , Hydrophilic_hsdc , Acid_hsdc , Base_hsdc = [] ,[] ,[] ,[]
            for x in range(len(residue)):
                if float(buri[x]) != 0:
                    hsdc = tree.xpath('//interface[%d]//molecule[%d]//residues/residue[%d]/bonds/text()'%(label,y,x+1))
                    if residue[x] in Hydrophobic:
                        Hydrophobic_energy -= float(energy[x])
                        Hydrophobic_count += 1
                        
                        if len(hsdc) != 0:
                            Hydrophobic_hsdc.append(ser_no[x]+' '+residue[x]+' '+str(hsdc[0]))

                            
                        # if hsdc[x] != ' ':
                        #     Hydrophobic_hsdc.append(hsdc[x])
                        # print(residue[x][0][1],residue[x][-1],energy[x])
                    elif residue[x] in Hydrophilic:
                        Hydrophilic_energy -= float(energy[x])
                        Hydrophilic_count += 1
                        # if hsdc[x] != ' ':
                        #     Hydrophilic_hsdc.append(hsdc[x])
                        if len(hsdc) != 0:
                            Hydrophilic_hsdc.append(ser_no[x]+' '+residue[x]+' '+str(hsdc[0]))
        
                    elif residue[x] in Base:
                        Base_energy -= float(energy[x])
                        Base_count += 1
                        if len(hsdc) != 0:
                            Base_hsdc.append(ser_no[x]+' '+residue[x]+' '+str(hsdc[0]))
                            
                    elif residue[x] in Acid:
                        Acid_energy -= float(energy[x])
                        Acid_count += 1
                        if len(hsdc) != 0:
                            Acid_hsdc.append(ser_no[x]+' '+residue[x]+' '+str(hsdc[0]))
                    else:
                        else_energy -= float(energy[x])
            total_energy = Hydrophobic_energy + Hydrophilic_energy + Acid_energy + Base_energy
            if else_energy != 0 and total_energy == 0:
                information = 'ligand_energy:%f'%(else_energy)
            else:
                information = 'Hydrophobic_energy:%f quantity:%d   hsdc:%s\nHydrophilic_energy:%f quantity:%d   hsdc:%s\nBase_energy:%f quantity:%d   hsdc:%s\nAcid_energy:%f quantity:%d   hsdc:%s' %(Hydrophobic_energy,Hydrophobic_count,str(Hydrophobic_hsdc),Hydrophilic_energy,Hydrophilic_count,str(Hydrophilic_hsdc),Base_energy,Base_count,str(Base_hsdc),Acid_energy,Acid_count,str(Acid_hsdc))
            with open('pisa2_force.txt','a') as fp:
                # fp.write(str(key)+ '\n')
                # fp.write(ligand + '\n')
                # fp.write(Analysis + '\n')
                fp.write(chain_id[0] + '\n' + information + '\n')
        with open('pisa2_force.txt','a') as fp:
            fp.write('***********************************************************\n')
    print(all_label)    
        
            
with open('dimer_pdb.txt','r') as fp:  
    pdb_list = list(fp.readlines())
count = 0
start = 0
for pdb in pdb_list:
    # print(pdb.split()[-1])

    if '2ijj' in pdb.split():
        start = 1
        continue
    if 'notfound' not in pdb.split()[-1] and start == 1:
        # break
        # count += 1
        # print(str(count)+'/8305')
        print(pdb.split()[0])
        pisa(pdb.split()[0],pdb.split()[-1])
        sleep(2)
#     # if count == 10:
#     #     break
# force([1],'4o23')
# pisa('5bwy')
