with open('./pisa2_result_highdletaG.txt','r') as fp:
    protein_id = fp.read().split('\n')
for x in range(len(protein_id)):
    protein_id[x] = protein_id[x].split('@')[0]
# print(protein_id)
with open('./pisa2_force.txt','r') as fp:
    protein_detail = fp.read().split('-----------------------------------------------------------')
    protein_detail.pop(0)
protein_highGList = []
plabel ,p_all_label,bic= [],[],[]
for x in range(len(protein_detail)):
    y = protein_detail[x].split('***********************************************************')
    if y[0].split('\n')[1].split('@')[0] in protein_id:
        protein_highGList.append(x)
        minenergy = 0
        minlabel = ''
        for z in y:
            z.split('ligand_energy:')[0].split('\n')[1]
            try:
                if minenergy > float(z.split('ligand_energy:')[1].split('\n')[0]):
                    minlabel = z.split('ligand_energy:')[0].split('\n')[1]
                    if len(minlabel) == 1:
                        minlabel = z.split('ligand_energy:')[0].split('\n')[-2]
                    minenergy = float(z.split('ligand_energy:')[1].split('\n')[0])
            except:
                pass
        # if 'ZN' not in minlabel and 'K'not in minlabel and 'SO4' not in minlabel:
        #     print(y[0].split('\n')[1].split('@')[0])    
        try:
            iron = minlabel.split('[')[1].split(']')[0]
        except:
            bic.append(y[0].split('\n')[1].split('@')[0])
            continue
        p_all_label.append(iron)
        if minlabel.split('[')[1].split(']')[0] not in plabel:
            plabel.append(iron)
        
        # print(minlabel.split('[')[1].split(']')[0],y[0].split('\n')[1].split('@')[0],minenergy)

for x in plabel:
    print(x,p_all_label.count(x))
print(bic)
# print(protein_highGList)

            
# for x in protein_highGList:
#     y = protein_detail[x].split('***********************************************************')
    
    
    
# count = 0
# protein_highGList = []
# for x in protein_detail:
    
#     if protein_detail[1][0].split('\n')[1].split('@')[0] in protein_id:
#         protein_highGList
#     count += 1

