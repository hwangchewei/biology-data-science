from posixpath import split
from time import sleep
from Biopython import biopython
from multiprocessing import Process, Pool
import multiprocessing as mp
import requests
from lxml import etree
import tqdm
import os



def protein_check(id_force):
    id = id_force[0]
    try:
        force = round(float(id_force[1]),2)
    except:
        force = id_force[1]
    if not os.path.isdir('D:/PDB/dimer/'+id+'.pdb'):
        url = 'https://files.rcsb.org/view/%s.pdb'%(id)
        response = requests.get(url = url).text.split('\n')

        
        with open('D:/PDB/dimer/'+id+'.pdb','w') as fp:
            for rLine in response:
                if 'ATOM' in rLine:
                    fp.write(rLine + '\n')
    protein = biopython('D:/PDB/dimer/'+id+'.pdb')
    result = protein.surface_area_amount(4.5,900)
    with open('./Biopython/face_dimer3.txt','a') as fp:
        fp.write(id +' '+ str(result) +' '+ str(force) + ' \n')
    sleep(10)
    return

def get_protein():
    with open('pisa2_summary.txt','r') as fp:
        id_list = fp.read().split('-----------------------------------------------------------\n')
    protein_id = []
    try:
        with open('./Biopython/face_dimer3.txt','r') as fp:
            already_check = fp.read()
    except:
        already_check = []
    for id in id_list:
        if str(id.split('@')[0]) in already_check:
            continue
        else:
            protein_id.append((id.split('@')[0],id.split('\n')[-2]))
    return protein_id
    
    

def main():
    protein = get_protein()
    cpus = int(mp.cpu_count()/2)
    pool = Pool(2)
    r = list(tqdm.tqdm(pool.imap(protein_check,protein),total=len(protein)))
    pool.close()
    pool.join()
if __name__ == '__main__':
    main()