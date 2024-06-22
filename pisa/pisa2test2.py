from lxml import etree

with open('pisa2_force.txt','r') as fp:
    all_result = fp.read().split('-----------------------------------------------------------')
    all_result.pop(0)
with open('pisa2_result_lowdletaG.txt','r') as fp:
    protein_id = fp.read().strip().split('\n')

for x in range(len(protein_id)):
    protein_id[x] = protein_id[x].split('@')[0]
    
count,count1 ,cyscount,sercount,tyrcount= 0 ,0 ,0 ,0 ,0
count1_list ,elselist= [] ,[]
for result in all_result:
    result = result.strip().split('***********************************************************')
    result = result[0].strip().split('\n')
    if result[0].split('@')[0] in protein_id:

        try:
            Hydrophobic = float(result[2].split()[0].split(':')[1])+float(result[7].split()[0].split(':')[1])
            Hydrophilic = float(result[3].split()[0].split(':')[1])+float(result[4].split()[0].split(':')[1])+float(result[5].split()[0].split(':')[1])+float(result[8].split()[0].split(':')[1])+float(result[9].split()[0].split(':')[1])+float(result[10].split()[0].split(':')[1])
            count+=1
            if (Hydrophobic-Hydrophilic) / Hydrophobic >= 0.5: #疏水為親水的兩倍
                count1+=1
            else:
                # if 'CYS' in result[3].split('[')[-1]:
                #     cyscount += 1
                # elif 'TYR' in result[3].split('[')[-1]:
                #     tyrcount += 1
                # elif 'SER' in result[3].split('[')[-1]:
                #     sercount += 1
                
                # else:
                #     elselist.append(result[0]+' '+result[3])
                    
                count1_list.append(result[0].split('@')[0])
                if Hydrophobic<Hydrophilic:
                    print(result[0].split('@')[0])
            with open('lowdeltaG2_force','a') as fp:
                fp.write('%s\nHydrophobic:%f\nHydrophilic:%f\n'%(result[0].split('@')[0],Hydrophobic,Hydrophilic))
        except:
            with open('lowdeltaG2_force','a') as fp:
                fp.write('%s\nretry\n'%(result[0].split('@')[0]))
        
        # print(result[0])
residuelist ,minresiduelist1,minresiduelist2= [],[],[]
for x in count1_list:
    x = str(x)
    # print(x)
    tree = etree.parse('./pisa1/%s/%s.xml'%(x,x))
    Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR','ASP','GLU','HIS','LYS','ARG']
    for y in range(1,3):
        chain_id = tree.xpath('//interface[1]//molecule[%d]/chain_id/text()'%(y))
        residue = tree.xpath('//interface[1]//molecule[%d]//residues/residue/name/text()'%(y))
        buri = tree.xpath('//interface[1]//molecule[%d]//residues/residue/bsa/text()'%(y))
        energy = tree.xpath('//interface[1]//molecule[%d]//residues/residue/solv_en/text()'%(y))
        ser_no = tree.xpath('//interface[1]//molecule[%d]//residues/residue/ser_no/text()'%(y))
        r_name = tree.xpath('//interface[1]//molecule[%d]//residues/residue/name/text()'%(y))
        minenergy,renergy = 0 ,0
        minresidue = ''
        
        for x in range(len(residue)):
            if float(buri[x]) != 0:
                if residue[x] in Hydrophilic:
                    # print(residue[x])
                    renergy += float(energy[x])
                    if minenergy < float(energy[x]):
                        minenergy = float(energy[x])
                        minresidue = residue[x]
                            # if hsdc[x] != ' ':
                            #     Hydrophilic_hsdc.append(hsdc[x])
        # if y == 1:
        minresiduelist1.append(minresidue)
        # else:
        #     minresiduelist2.append(minresidue)
        if minresidue not in residuelist:
            residuelist.append(minresidue)
    #     print((minenergy/renergy) >= 0.4,minresidue)
    # print(str(x)+ ' ' + str(minresidue)+' '+str(minenergy))                    
for x in residuelist:
    print(x,minresiduelist1.count(x))
    # print(x,minresiduelist2.count(x))
    # break
with open('lowdeltaG2_force','a') as fp:
    fp.write(str(float(count1/count))+'\n') #+'CYS: '+str(float(cyscount/count))+'\n'+'SER: '+str(float(sercount/count))+'\n'+'TYR: '+str(float(tyrcount/count))+'\n')
    # for x in elselist:
    #     fp.write(str(x)+'\n')
    for x in count1_list:
        fp.write(str(x)+'\n')
