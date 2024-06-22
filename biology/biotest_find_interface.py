from operator import ne
from time import sleep
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface,residue_depth
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.NeighborSearch import NeighborSearch
import Bio
from Bio.PDB.vectors import calc_angle,calc_dihedral
from biotest7 import asa_total_dict
import numpy
import copy







# def min_dist(coord, surface):
#     """Return minimum distance between coord and surface."""
#     d = surface - coord
#     d2 = numpy.sum(d * d, 1)
#     return numpy.sqrt(min(d2))
# def residue_depth(residue, surface):
#     """Residue depth as average depth of all its atoms.

#     Return average distance to surface for all atoms in a residue,
#     ie. the residue depth.
#     """
#     atom_list = residue.get_unpacked_list()
#     length = len(atom_list)
#     d = 100000
#     atom_d = ''
#     for atom in atom_list:
#         coord = atom.get_coord()
#         if d > min_dist(coord, surface) and str(atom) != '<Atom O>':
#             d = min_dist(coord, surface)
#             atom_d = atom
            
#         # d = d + min_dist(coord, surface)
#     return (d,atom_d)
def vector_centersite(r_list):
    x,y,z = 0 ,0 ,0
    for r in r_list:
        r = r['CA'].get_vector()
        x += r[0]
        y += r[1]
        z += r[2]
    return (x,y,z)
def sdistance(surface1,surface2):
    # print(surface1)
    x = surface1[0] - surface2[0]
    y = surface1[1] - surface2[1]
    z = surface1[2] - surface2[2]
    return (x**2+y**2+z**2)**0.5
def surface_distance(surface_dict,site_dict):
    surface_distance_dict = {}
    for x in range(len(surface_dict)):
        surface1 = site_dict[x]
        surface_dict.pop(x)
        site_dict.pop(x)
        d_list = []
        # print(total_surface_dict)

        for y in surface_dict:
            surface2 = site_dict[y]
            d = sdistance(surface1,surface2)
            d_list.append((surface_dict[y],d))   
        surface_distance_dict.update({x:d_list})    
    return surface_distance_dict
def surface_angle(total_surface_dict):
    angle_list = []
    for x in range(len(total_surface_dict)-2):
        p1 = (total_surface_dict[x][1][0],total_surface_dict[x][1][1],total_surface_dict[x][1][2])
        p2 = (total_surface_dict[x+1][1][0],total_surface_dict[x+1][1][1],total_surface_dict[x+1][1][2])
        p3 = (total_surface_dict[x+2][1][0],total_surface_dict[x+2][1][1],total_surface_dict[x+2][1][2])
        v1 = numpy.array((p1[0] - p2[0],p1[1] - p2[1],p1[2] - p2[2]))
        v2 = numpy.array((p3[0] - p2[0],p3[1] - p2[1],p3[2] - p2[2]))
        cosangle = numpy.dot(v1,v2)/(numpy.sqrt(v1.dot(v1))*numpy.sqrt(v2.dot(v2)))
        angle_list.append(numpy.degrees(numpy.arccos(cosangle)))
        print(numpy.degrees(numpy.arccos(cosangle)),total_surface_dict[x])
    
    return angle_list
def surface_area(neighbor,residue,surface_residue_atom,surface_residue):
    near = []
    for x in range(1,7):
        y = ns.search(surface_residue_atom[residue].coord,x,'R')
        for z in ns.search(surface_residue_atom[residue].coord,x,'R'):
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
    p1 = surface_residue_atom[residue].get_vector()
    p2 = surface_residue_atom[near[0]].get_vector()
    p3 = surface_residue_atom[near[1]].get_vector()
    # print(surface_residue_atom[residue],surface_residue_atom[near[0]],surface_residue_atom[near[1]])
    # print(residue,near)
    surface_area_ = [residue,near[0],near[1]]
    # print(surface_area_)
    for x in neighbor:
        if x in surface_area_:
            # print(x)
            continue
        p4 = surface_residue_atom[x].get_vector()
        if abs(numpy.degrees(calc_dihedral(p1,p2,p3,p4))) > 170 or abs(numpy.degrees(calc_dihedral(p1,p2,p3,p4))) < 10:
            surface_area_.append(x)
    return surface_area_
def total_r_asa(residue_list):
    total_asa = 0
    for residue in residue_list:
        asa = 0
        for atom in residue:
            asa += round(atom.sasa, 2)
        total_asa += asa
    return total_asa

# with open('CA_fasta_dimer.txt','r') as fp:
#     protein_list = fp.read().split('\n')
bicface ,mayface = 0,0
protein = input('root :')
protein = protein+'.pdb'

# url = 'https://files.rcsb.org/view/%s.pdb'%(protein)
# response = requests.get(url = url).text
# tree = etree.HTML(response)
# with open('D:/PDB/CA/'+protein,'w') as fp:
#     fp.write(response)
amino_acid = ['CYS','ASN','GLN','SER','THR','TYR','ASP','GLU','HIS','LYS','ARG','ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
Hydrophilic = ['CYS','ASN','GLN','SER','THR','TYR','ASP','GLU','HIS','LYS','ARG']
Hydrophobic = ['ALA','PHE','ILE','LEU','MET','PRO','VAL','TRP','GLY']
parser = PDBParser(QUIET = 1)
structure = parser.get_structure("test",protein)
model = structure[0]
try:
    chain_ = model['A'].get_list()
except:
    print(protein+' noA')

chain = []
for x in chain_:
    if x.get_id()[0] == ' ':
         chain.append(x)

surface = get_surface(chain)
atoms = Bio.PDB.Selection.unfold_entities(model['A'], 'A')
ns = Bio.PDB.NeighborSearch(atoms)  
sr = ShrakeRupley(n_points=1000)
sr.compute(model['A'], level="C")
# sr = ShrakeRupley(probe_radius=1.5, n_points=1000)
# sr.compute(chain, level="C")
# asa_dict = asa_total_dict(chain)
surface_atom,surface_residue = [],[]
total_len = len(chain)
r_count = 0
for residue in chain: 
    r_count+=1
    if r_count < 10:
        continue
    if r_count >= total_len - 10:
        break
    if str(residue).split(' ')[1] not in amino_acid:
        continue
    rd = residue_depth(residue, surface)
    # print(rd ,str(residue).split(' ')[1],str(residue).split('=')[2].split(' ')[0])
    if rd[0] < 4.5:
        surface_atom.append(rd[1])
        surface_residue.append(residue)
        # ns.search(residue['CA'].coord,5,'R')
# print(surface_residue)
surface_residue_atom = dict(zip(surface_residue,surface_atom))
surface_atom = dict(enumerate(surface_atom))
surface_residue = dict(enumerate(surface_residue))

surface_area_list ,surface_id_list = [] ,[]
            
for x in surface_atom:
    z = []
    for y in ns.search(surface_atom[x].coord,30,'R'):
        if y in surface_residue.values():  
            z.append(y)
    # print(z)
    # print(surface_residue[x].sasa)
    surface = surface_area(z,surface_residue[x],surface_residue_atom,surface_residue.values())
    surface_id = []
    for y in surface:
        id = int(str(y).split('=')[2].split(' ')[0])
        surface_id.append(id)
    surface_id.sort()
    # print(surface_id)
    if surface_id not in surface_id_list:
        surface_area_list.append(surface)
        surface_id_list.append(surface_id)
surface_list, surface_id_list,biccount,liccount,maycount= [],[],0,0,0
for x in surface_area_list:
    surface_id = []
    for y in x:      
            # print(y)
        id = int(str(y).split('=')[2].split(' ')[0])
        surface_id.append(id)  
    surface_id.sort()
    surface_list.append(x)
    surface_id_list.append(surface_id)
    biclist,liclist,total_len = [],[],len(chain)
    for y in x:
        # id = int(str(y).split('=')[2].split(' ')[0])
        # if id < 10:
        #     continue
        # if id >= total_len - 10:
        #     break
        residue_name = str(y).split(' ')[1] 
        if residue_name in Hydrophilic:
            liclist.append(y)
        elif residue_name in Hydrophobic:
            biclist.append(y)
    bic_asa = total_r_asa(biclist)
    lic_asa = total_r_asa(liclist)
    surface_id = []
    # print(x)
    total_asa = total_r_asa(x)
    if total_asa >= 950 and (bic_asa - lic_asa) > 100:
        # print(x)
        surface_id = [] 
        for y in x:      
            # print(y)
            id = int(str(y).split('=')[2].split(' ')[0])
            surface_id.append(id)  
        surface_id.sort()
        # print(surface_id,bic_asa,lic_asa)
        # surface_list.append(x)
        # surface_id_list.append(surface_id)
        biccount += 1
    elif total_asa >= 950 and (lic_asa - bic_asa) > 100:
        liccount += 1
        surface_id = [] 
        for y in x:      
            # print(y)
            id = int(str(y).split('=')[2].split(' ')[0])
            surface_id.append(id)  
        surface_id.sort()
        # surface_list.append(x)
        # surface_id_list.append(surface_id)
        # print(surface_id,bic_asa,lic_asa,total_asa)
    elif total_asa >= 800 and abs(lic_asa - bic_asa) < 100:
        maycount += 1
        surface_id = [] 
        for y in x:      
            # print(y)
            id = int(str(y).split('=')[2].split(' ')[0])
            surface_id.append(id)  
        surface_id.sort()
        # surface_list.append(x)
        # surface_id_list.append(surface_id)
        # print(surface_id,bic_asa,lic_asa,total_asa)
    
print(biccount,liccount,maycount,' '+protein)
if biccount > 0:
    bicface += 1
elif maycount > 0:
    mayface += 1




'''
for x in range(len(surface_list)):
    # print('a')

    surface1 = surface_list[x]

    for y in surface_list:
        if y != surface1:
            surface2 = y
        else:
            continue
        if len(y) < 3:
            break
        p1 = surface_residue_atom[surface1[0]].get_vector()
        p2 = surface_residue_atom[surface1[int(len(surface1)/2)]].get_vector()
        p3 = surface_residue_atom[surface1[-1]].get_vector()
        p4 = surface_residue_atom[surface2[0]].get_vector()
        p5 = surface_residue_atom[surface2[int(len(surface2)/2)]].get_vector()            
        p6 = surface_residue_atom[surface2[-1]].get_vector()
        # print(p1,p2,p3,p4,p5,p6)

        a1 = abs(numpy.degrees(calc_dihedral(p1,p2,p3,p4)))
        a2 = abs(numpy.degrees(calc_dihedral(p1,p2,p3,p5)))
        a3 = abs(numpy.degrees(calc_dihedral(p1,p2,p3,p6)))

        if (a1+a2+a3) / 3 > 170 or (a1+a2+a3) / 3 < 10:
            for z in surface2:
                id = int(str(z).split('=')[2].split(' ')[0])
                
                if z not in surface_list[x]:
                    # print(z,id)

                    surface_list[x].append(z)
                    surface_id_list[x].append(id)
            # print(surface_id_list[x])
    # print(surface_id_list[x])
    # print(len(surface_list[x]),len(surface_id_list[x]))
biccount,liccount,maycount = 0,0,0            
for x in surface_list:
    biclist,liclist,total_len = [],[],len(chain)
    for y in x:
        # id = int(str(y).split('=')[2].split(' ')[0])
        # if id < 10:
        #     continue
        # if id >= total_len - 10:
        #     break
        residue_name = str(y).split(' ')[1] 
        if residue_name in Hydrophilic:
            liclist.append(y)
        elif residue_name in Hydrophobic:
            biclist.append(y)
    bic_asa = total_r_asa(biclist)
    lic_asa = total_r_asa(liclist)
    surface_id = []
    # print(x)
    total_asa = total_r_asa(x)

    if total_asa >= 1000 and (bic_asa - lic_asa) > 100:
        # print(x)
        surface_id = [] 
        for y in x:      
            # print(y)
            id = int(str(y).split('=')[2].split(' ')[0])
            surface_id.append(id)  
        surface_id.sort()
        # print(surface_id,bic_asa,lic_asa)
        # surface_list.append(x)
        # surface_id_list.append(surface_id)
        biccount += 1
    elif total_asa >= 1000 and (lic_asa - bic_asa) > 100:
        liccount += 1
        surface_id = [] 
        for y in x:      
            # print(y)
            id = int(str(y).split('=')[2].split(' ')[0])
            surface_id.append(id)  
        surface_id.sort()
        # surface_list.append(x)
        # surface_id_list.append(surface_id)
        # print(surface_id,bic_asa,lic_asa,total_asa)
    elif total_asa >= 1000 and abs(lic_asa - bic_asa) < 100:
        maycount += 1
        surface_id = [] 
        for y in x:      
            # print(y)
            id = int(str(y).split('=')[2].split(' ')[0])
            surface_id.append(id)  
        surface_id.sort()
        # surface_list.append(x)
        # surface_id_list.append(surface_id)
        # print(surface_id,bic_asa,lic_asa,total_asa)
print(biccount,liccount,maycount)
'''
