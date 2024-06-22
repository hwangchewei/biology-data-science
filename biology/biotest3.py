from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from time import sleep



with open('pdb_id_CA.txt','r') as fp:
    pdb_id = fp.read().split(',')
    
T = False
for x in pdb_id:
    if x == '5EOI':
        T = True
        continue
    if T:
        
        url = ('https://www.ncbi.nlm.nih.gov/protein/%s_A?report=fasta'%(str(x)))
        print(url)
        chrome_option = Options()
        chrome_option.add_argument('--headless')
        chrome_option.add_argument('--disable-gpu')
        chrome_option.add_argument('--disable-javascript')
        chrome_option.add_argument('blink-settings=imagesEnabled=false')
        bro = webdriver.Chrome(executable_path='./chromedriver',chrome_options=chrome_option)
        bro.get(url=url)
        # sleep(1)
        def fasta_try(count):
            if count < 3:    
                try:
                    fasta = bro.find_element(By.XPATH,'/html/body/div[1]/div[1]/form/div[1]/div[4]/div/div[5]/div[2]/div[1]/pre').text
                except:
                    sleep(0.2)
                    return fasta_try(count+1)
            else:
                fasta = x+' time_out'
            return fasta  
        fasta = fasta_try(0)      
        print(fasta)
        bro.quit()
        if 'time_out' in fasta:
            with open('CA_fasta_retry.txt','a+') as fp:
                fp.write(fasta+'\n'+'.............'+'\n')
        else:
            with open('CA_fasta.txt','a+') as fp:
                fp.write(fasta+'\n'+'.............'+'\n')
        sleep(0.2)
