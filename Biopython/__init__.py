# from operator import ne, truediv
# from select import select
# from time import sleep
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface,min_dist
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.ResidueDepth import ResidueDepth
from Bio.PDB.NeighborSearch import NeighborSearch
import Bio
from Bio.PDB import *
from Bio.PDB.vectors import calc_dihedral
import numpy
import copy
# from multiprocessing.dummy import Pool as threadpool

aa_codes = {
    'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
    'PHE':'F','GLY':'G','HIS':'H','LYS':'K',
    'ILE':'I','LEU':'L','MET':'M','ASN':'N',
    'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
    'THR':'T','VAL':'V','TYR':'Y','TRP':'W'}  

amino_acid = ['CYS','ASN','GLN','SER','THR','TYR','ASP','GLU','HIS','LYS','ARG','ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
# Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR','ASP','GLU','HIS','LYS','ARG']
# Hydrophobic = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
Charge = ['ASP','GLU','HIS','LYS','ARG']
Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR']
Hydrophobic = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
amino_weight_funtion = {
    'CYS':1.12,
    'ASN':-0.05,
    'GLN':0.16,
    'SER':-0.39,
    'THR':0.46,
    'TYR':0.01,
    'ASP':-0.28,
    'GLU':-0.3,
    'HIS':0.2,
    'LYS':0.07,
    'ARG':0.02,
    'ALA':-0.06,
    'PHE':0.43,
    'ILE':0.7,
    'LEU':1.01,
    'MET':0.86,
    'PRO':0.83,
    'VAL':0.38,
    'TRP':1.69,
    'GLY':-0.25}
amino_weight_funtion2 = {
    'CYS':1.08,
    'ASN':-0.16,
    'GLN':0.10,
    'SER':-0.45,
    'THR':0.35,
    'TYR':-0.03,
    'ASP':-0.34,
    'GLU':-0.37,
    'HIS':0.05,
    'LYS':0,
    'ARG':-0.04,
    'ALA':-0.17,
    'PHE':0.3,
    'ILE':0.49,
    'LEU':0.85,
    'MET':0.66,
    'PRO':0.59,
    'VAL':0.26,
    'TRP':1.46,
    'GLY':-0.43}

class biopython():
    def __init__(self,protein,chain_name = 'A') -> None:
        self.surface_residue_atom ,self.surface_atom ,self.surface_residue_dict ,self._asa = {},{},{},1200
        self.maxbicinterface = ''
        parser = PDBParser(QUIET = 1)
        self.structure = parser.get_structure("test",protein)
        self.model = self.structure[0]
        try:
            chain_ = self.model[chain_name].get_list()
            self.check_chain = 1
            # print(chain_)
        except:
            # print(protein+' noA')
            self.check_chain = -1
            return
        self.chain = []
        for x in chain_:
            if x.get_id()[0] == ' ':
                self.chain.append(x)
        try:        
            self.surface = get_surface(self.chain)
            self.check_surface = 1
        except:
            self.check_surface = -1
            return 
        atoms = Bio.PDB.Selection.unfold_entities(self.chain, 'A')
        self.ns = Bio.PDB.NeighborSearch(atoms)  
        self.sr = ShrakeRupley()
        self.sr.compute(self.model[chain_name], level="C")
    def residue_depth(self,residue):
        atom_list = residue.get_unpacked_list()
        length = len(atom_list)
        d = 100000
        atom_d = ''
        for atom in atom_list:
            coord = atom.get_coord()
            # print(atom)
            if d > min_dist(coord, self.surface):
                d = min_dist(coord, self.surface)
                atom_d = atom
                
            # d = d + min_dist(coord, surface)
        return (d,atom_d)
    def surface_residue(self,distance):
        total_len = len(self.chain)
        surface_atom ,surface_residue = [],[]
        r_count = 0
        for residue in self.chain: 
            # r_count+=1           // 去除首尾
            # if r_count < 5:
            #     continue
            # if r_count >= total_len - 5:
            #     break
            # if str(residue).split(' ')[1] not in amino_acid:
            #     continue
            rd = self.residue_depth(residue)
            # print(rd,residue)
            # rd = ResidueDepth.residue_depth(residue)
            # return

            # print(rd[0],residue)
            # print(rd ,str(residue).split(' ')[1],str(residue).split('=')[2].split(' ')[0])
            if rd[0] <= distance:
                # print(rd[1])
                surface_atom.append(rd[1])
                surface_residue.append(residue)
        self.surface_residue_atom = dict(zip(surface_residue,surface_atom))
        self.surface_atom = dict(enumerate(surface_atom))
        self.surface_residue_dict = dict(enumerate(surface_residue))
        
        return surface_residue
    def find_surface(self,distance):
        self.surface_residue(distance)
        self.surface_area_list ,self.surface_area_id_list = [],[]
        radius = ((round(self.total_r_asa(self.surface_residue_dict.values()),0)/6)/3.14)**0.5+2
        # print(self.total_r_asa(self.surface_residue_dict.values()),radius)
  
        for x in self.surface_atom:
            z = []
            for y in self.ns.search(self.surface_atom[x].coord,radius,'R'):
                if y in self.surface_residue_dict.values():  
                    z.append(y)
        
            # print(z)
            # print(surface_residue[x].sasa)
            surface = self.surface_area(z,self.surface_residue_dict[x],self.surface_residue_atom,self.surface_residue_dict.values())
            surface_id = []
            for y in surface:
                id = int(str(y).split('=')[2].split(' ')[0])
                surface_id.append(id)
            surface_id.sort()
            # print(surface_id)
            if surface_id not in self.surface_area_id_list:
                self.surface_area_list.append(surface)
                self.surface_area_id_list.append(surface_id)
        return self.surface_area_list
    def surface_area(self,neighbor,residue,surface_residue_atom,surface_residue):
        near = []
        for x in range(1,10):
            y = self.ns.search(surface_residue_atom[residue].coord,x,'R')
            for z in self.ns.search(surface_residue_atom[residue].coord,x,'R'):
                if z in surface_residue and z != residue:
                    if z not in near:
                        near.append(z)  
                    if len(near) == 2:
                        break
            if len(near) == 2:
                break
        else:
            return [residue]
        # print(near)
        coordp1 = surface_residue_atom[residue].get_vector()
        p1 = surface_residue_atom[residue].get_vector()
        p2 = surface_residue_atom[near[0]].get_vector()
        p3 = surface_residue_atom[near[1]].get_vector()
        # print(surface_residue_atom[residue],surface_residue_atom[near[0]],surface_residue_atom[near[1]])
        # print(residue,near)
        surface_area_ = [residue,near[0],near[1]]
        # print(surface_area_)
        non_surface_residue ,sneighbor = 0,0
        pdegrees,ndegrees = 0,0
        for x in neighbor:
            coordp4 = surface_residue_atom[x].get_coord()
            if x in surface_area_ or x not in self.surface_residue_dict.values():
                # print(x)
                continue
            # p1_p4_dis = ((coordp1[0]-coordp4[0])**2+(coordp1[1]-coordp4[1])**2+(coordp1[2]-coordp4[2])**2)**0.5
            # print(p1_p4_dis)
            # if p1_p4_dis > ((round(self.total_r_asa(self.surface_residue_dict.values()),0)/6)/3.14)**0.5:
            #     continue
            p4 = surface_residue_atom[x].get_vector()
            sneighbor += 1
            p1p4degrees = numpy.degrees(calc_dihedral(p1,p2,p3,p4))
            if abs(p1p4degrees) > 165 or abs(p1p4degrees) < 15:
                surface_area_.append(x)
            elif p1p4degrees > 0:
                pdegrees += 1
            elif p1p4degrees < 0:
                ndegrees += 1
        if pdegrees <= 5 or ndegrees <= 5:
            # print(pdegrees,ndegrees)
            return surface_area_
        else:
            return [residue,near[0],near[1]]
        # print(sneighbor,non_surface_residue,len(surface_area_))
        
    def surface_area_amount(self,distance=2,amino_weight_function=200):
        # if self.check_chain < 0 or self.check_surface < 0:
        #     return (0,)
        
        self.find_surface(distance)
        print('Protein Solvent Accessible Surface Area： %s'%(self.total_r_asa(self.surface_residue_dict.values())))
        total_P_asa = round(self.total_r_asa(self.surface_residue_dict.values()),2)
        mayx_bicasa,may_asa,maxbic_asa,maxlic_asa,max_asa,max_surface_id = 0,[0,'0'],[0,'0'],[0,'0'],0,[]
        surface_id_list ,surface_list ,liccount ,biccount ,maycount = [] ,[] ,0 ,0 ,0
        max_weight_f = 0
        check_ = False
        print('可能形成聚體的胺基酸: ')
        for x in self.surface_area_list:
            surface_id = []
            for y in x:      
                    # print(y)
                id = int(str(y).split('=')[2].split(' ')[0])
                surface_id.append(id)  
            surface_id.sort()
            # surface_list.append(x)
            surface_id_list.append(surface_id)
            biclist,liclist,chargelist,total_len,charge_asa = [],[],[],len(self.chain),[]
            amino_weight_funtion_area,amino_weight_funtion_area2,a_weight_f,a_weight_f2 = 0,0,0,0
            for y in x:
                # id = int(str(y).split('=')[2].split(' ')[0])       //不計算頭尾
                # if id < 10:
                #     continue
                # if id >= total_len - 10:
                #     break
                residue_name = str(y).split(' ')[1] 
                if residue_name in Hydrophilic:
                    liclist.append(y)
                elif residue_name in Hydrophobic:
                    biclist.append(y)
                elif residue_name in Charge:
                    chargelist.append(y)
                # a_weight_f += amino_weight_funtion[residue_name]
                a_weight_f2 += amino_weight_funtion2[residue_name]
                # print(residue_name,biclist,liclist,chargelist)
                # amino_weight_funtion_area += round(self.total_r_asa([y]),2)*amino_weight_funtion[residue_name]
                amino_weight_funtion_area2 += self.total_r_asa([y])*amino_weight_funtion2[residue_name]
            amino_weight_funtion_area2 = round(amino_weight_funtion_area2,2)
            bic_asa = round(self.total_r_asa(biclist),2)
            lic_asa = round(self.total_r_asa(liclist),2)
            charge_asa = round(self.total_r_asa(chargelist),2)
            
            surface_id = []
            # print(x)
            total_s_asa = round(self.total_r_asa(x),2)
            # if total_P_asa/6 < inputasa:
            compare_asa = total_P_asa/6
            # print(total_s_asa,compare_asa)
            # return
            # else:
            #     compare_asa = inputasa
            # print(total_P_asa)
            # if total_s_asa > max_asa:
            #     max_asa = total_s_asa
            #     # print(max_asa)
                
            #     self.maxbicinterface = x
            # if (total_s_asa >= compare_asa and (bic_asa>lic_asa)) or (bic_asa > lic_asa*2.5 and bic_asa > compare_asa/2) or (total_s_asa > 2500 and (bic_asa  > total_s_asa*0.4)): 
            # /...舊判斷條件
            # if (bic_asa + lic_asa*3/4) > (lic_asa/4 + charge_asa*3) and total_s_asa > total_P_asa/8:
            #     surface_id = [] 
            #     for y in x:      
            #         # print(y)
            #         id = int(str(y).split('=')[2].split(' ')[0])
            #         surface_id.append(id)  
            #     surface_id.sort()
            #     # print(surface_id,bic_asa,lic_asa,charge_asa,a_weight_f,a_weight_f2)
            #     print(1)
            #     # surface_list.append(x)
            #     # surface_id_list.append(surface_id)
            #     if ((bic_asa + lic_asa*3/4) - (lic_asa/4 + charge_asa*3)) > maxbic_asa[0]:
            #         maxbic_asa[0] = total_s_asa
            #         maxbic_asa[1] = 'bic:%s,lic:%s,charge:%s'%(bic_asa,lic_asa,charge_asa)
            #         max_surface_id = [surface_id]
            #     biccount += 1
            # elif (total_s_asa >= compare_asa and (lic_asa > bic_asa) ) or (lic_asa > bic_asa*2.5 and lic_asa > compare_asa/2):
            #     liccount += 1
            #     print(2)
            #     if (lic_asa - bic_asa) > maxlic_asa[0]:
            #         maxlic_asa[0] = total_s_asa
            #         maxlic_asa[1] = 'bic: %s,lic: %s'%(bic_asa,lic_asa)
            # elif total_s_asa >= total_P_asa/10 and bic_asa > charge_asa:
            #     maycount += 1
            #     surface_id = [] 
            #     print(3)
            #     for y in x:      
            #         # print(y)
            #         id = int(str(y).split('=')[2].split(' ')[0])
            #         surface_id.append(id)  
            #     surface_id.sort()
            #     # print(surface_id)
            #     if total_s_asa > may_asa[0] and bic_asa > mayx_bicasa:
            #         mayx_bicasa = bic_asa
            #         may_asa[0] = total_s_asa
            #         # print(may_asa[0])
            #         may_asa[1] = 'bic:%s,lic:%s'%(bic_asa,lic_asa)
            
                    
        # surface_id = []
        # for y in self.maxbicinterface:      
            
        #     id = int(str(y).split('=')[2].split(' ')[0])
        #     surface_id.append(id)  
        # surface_id.sort()
        # print(surface_id)
        # print(max_surface_id)
        # return ([biccount,maycount],'protein_asa:%s maxbic_asa:%s  may_asa:%s'%(total_P_asa,maxbic_asa,may_asa))
        # return ([biccount,liccount,maycount],'protein_asa:%s maxbic_asa:%s maxlic_asa:%s may_asa:%s'%(total_P_asa,maxbic_asa,maxlic_asa,may_asa))
        # 舊判斷條件.../
            # print(a_weight_f,a_weight_f2,amino_weight_funtion_area,amino_weight_funtion_area2)
            if amino_weight_funtion_area2 >= amino_weight_function:
                check_ = True
                surface_id,surface_id2 = [] ,[]
                for y in x:      
                    # print(y)
                    residue_name = str(y).split(' ')[1] 
                    id = int(str(y).split('=')[2].split(' ')[0])
                    surface_id.append([id,aa_codes[residue_name]])  
                surface_id.sort()
                for x in surface_id:
                    surface_id2.append(str(x[0])+x[1])
                # print(surface_id,round(a_weight_f,2),round(a_weight_f2,2),round(amino_weight_funtion_area,2),round(amino_weight_funtion_area2,2))
                print(surface_id2,amino_weight_funtion_area2)
                if amino_weight_funtion_area2 > max_weight_f:
                    max_weight_f = amino_weight_funtion_area2
                    max_surface_id = surface_id2
        if check_ == False:
            print('無，需調低胺基酸權重')
        return ['最有可能形成聚體的胺基酸:%s '%(max_surface_id),max_weight_f]
    
    def total_r_asa(self,residue_list):
        total_asa = 0
        for residue in residue_list:
            asa = 0
            for atom in residue:
                # print(atom.sasa)
                asa += round(atom.sasa, 2)
            total_asa += asa
        return total_asa
    def find_bic_area(self,residues):
        for residue in residues:
            atom = self.surface_residue_atom[residue]
            radius_range = int(((round(self.total_r_asa(self.surface_residue_dict.values()),0)/6)/3.14)**0.5+1)
            bic_area_list = []
            if str(residue).split(' ')[1] in Hydrophobic:
                for radius in range(1,radius_range):
                    bic_area = []
                    has_lic = False
                    for nsresidue in self.ns.search(atom.coord,radius,'R'):
                        if nsresidue in self.surface_residue_dict.values():
                           
                            bic_area.append(nsresidue)
                            if str(nsresidue).split(' ')[1] in Hydrophilic:
                                has_lic = True
                                break
                    if has_lic:   
                        remove_r = []   
                        if len(bic_area) > 3:
                            p1 = self.surface_residue_atom[bic_area[0]].get_vector()
                            p2 = self.surface_residue_atom[bic_area[1]].get_vector()
                            p3 = self.surface_residue_atom[bic_area[2]].get_vector()
                            len_area = len(bic_area)
                            for x in range(3,len_area):
                                p4 = self.surface_residue_atom[bic_area[x]].get_vector()
                                if abs(numpy.degrees(calc_dihedral(p1,p2,p3,p4))) > 170 or abs(numpy.degrees(calc_dihedral(p1,p2,p3,p4))) < 10:
                                    pass
                                else:
                                    remove_r.append(bic_area[x])
                        print(bic_area,remove_r)
                        for x in remove_r:     
                            bic_area.remove(x)
                        print(bic_area)
                        bic_area_list.append(bic_area)
                        
                        break
            else:
                continue
        return bic_area_list
    def bic_area_amount(self,distance):
        bic_area_max_len = 0
        if self.check_chain < 0 or self.check_surface < 0:
            return (0,)
        
        residue = self.surface_residue(distance)
        bic_area_list = self.find_bic_area(residue)
        print(len(bic_area_list))
        bic_area_list2 = copy.deepcopy(bic_area_list)
        new_bic_area = []
        for bic_area in bic_area_list:
            id_list = []
            # print(len(x))
            for y in bic_area:
                y = int(str(y).split('=')[2].split(' ')[0])
                # print(y,id_list)
                if y in id_list:
                    continue
                id_list.append(y)
            id_list.sort()
            print(id_list)
            # bic_area_list2.remove(bic_area)
            for bic_area2 in bic_area_list2:
                if len(bic_area) > 3:
                    p1 = self.surface_residue_atom[bic_area[0]].get_vector()
                    p2 = self.surface_residue_atom[bic_area[1]].get_vector()
                    p3 = self.surface_residue_atom[bic_area[2]].get_vector()
                    if len(bic_area2) < 3:
                        p4 = self.surface_residue_atom[bic_area2[0]].get_vector()
                        if abs(numpy.degrees(calc_dihedral(p1,p2,p3,p4))) > 170 or abs(numpy.degrees(calc_dihedral(p1,p2,p3,p4))) < 10:
                            bic_area += bic_area2
                    if len(bic_area2) >= 3:
                        p4 = self.surface_residue_atom[bic_area2[0]].get_vector()
                        p5 = self.surface_residue_atom[bic_area2[int(len(bic_area2)/2)]].get_vector()
                        p6 = self.surface_residue_atom[bic_area2[-1]].get_vector()
                        ave_degree = (abs(numpy.degrees(calc_dihedral(p1,p2,p3,p4)))+abs(numpy.degrees(calc_dihedral(p1,p2,p3,p5)))+abs(numpy.degrees(calc_dihedral(p1,p2,p3,p6))))/3
                        if ave_degree > 170 or ave_degree < 10:
                            bic_area += bic_area2
            
                        
                else:
                    break
            new_bic_area.append(bic_area)
            # print(len(bic_area))
        # print(bic_area_list)
        # print(new_bic_area)
        bic_area_id_list = []
        for x in new_bic_area:
            id_list = []
            # print(len(x))
            for y in x:
                y = int(str(y).split('=')[2].split(' ')[0])
                # print(y,id_list)
                if y in id_list:
                    continue
                id_list.append(y)
            id_list.sort()
            # print(id_list)
            if id_list not in bic_area_id_list:
                bic_area_id_list.append(id_list)
        for x in bic_area_id_list:
            # print(x)
            pass
        

                

        return bic_area_id_list
            
               
                

 
#          def bic_area(surface_a):
#             z = [surface_a]
#             lic = False
#             for x in range(1,radius):
#                 for y in self.ns.search(surface_a.coord,radius,'R'):
#                     if y in self.surface_residue_dict.values():
#                         if y in Hydrophobic:
#                             z.append(y)
#                         else:
#                             lic = True
#                 if lic:
#                     return z         
# thread_pool = threadpool(3)
#         thread_pool.map(bic_area,self.surface_atom)