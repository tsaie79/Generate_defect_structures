__author__ = 'dchakraborty'


import sys, os, glob, shutil
from pymatgen.analysis.defects.core import Vacancy, Substitution, Interstitial
from pymatgen.analysis.defects.generators import VacancyGenerator, SubstitutionGenerator, \
                                         SimpleChargeGenerator,VoronoiInterstitialGenerator
from pymatgen.io.vasp import Poscar
from pymatgen.ext.matproj import MPRester
from pymatgen.core import PeriodicSite
from pymatgen.core import periodic_table
from string import digits


basepath=os.getcwd()
print(basepath)

#remove previous directories
for i in os.listdir():
	if os.path.isdir(i): shutil.rmtree(i) 

#sub_element_list=['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl']


'''Generate undefected files'''

with MPRester() as mp:
	ids = mp.get_all_substrates()
	print(ids)
	structures=[]
	for i in ids:

#Creating directories with mp-ids

		str_file=os.path.join(os.path.join(basepath,'%s')%(i))
		str_file_exists=os.path.exists(str_file)
		if not str_file_exists: os.mkdir(str_file)
		os.chdir(str_file)

#Creating bulk_str files with the structure_formula

		no_defect_data = open('undefected_data.dat','w')				
		bulk_structure = mp.get_structure_by_material_id(i) 
		no_defect_data.writelines(str(bulk_structure))
		no_defect_data.close()


#Finding elements and formula of the structure 

		elements=[]
		tot_elem=len(bulk_structure.types_of_specie)
		for ele in range(tot_elem):
			elements.append(str(bulk_structure.types_of_specie[ele]).split()[-1])
		print(elements)

	
#		for ele in elements[:]:
#			if ele in sub_element_list:
#				print(ele)
#				sub_element_list.remove(ele)
#		print(sub_element_list)


		count=0		
		with open('undefected_data.dat') as file:
			for line in file:
				count += 1
				if count == 2 : str_formula=line.split()[2]
		print(str(str_formula))
		os.rename('undefected_data.dat','%s'%(str_formula))
		structures.append(str(str_formula))


#Creating vacacies with the bulk structure formula

		defects_to_run=[]
		print('Vacancy Generator:')
		for vac_defect in VacancyGenerator( bulk_structure):
			print("\tCreated Defect {} at site {}".format( vac_defect.name, vac_defect.site))
			for charged_defect in SimpleChargeGenerator( vac_defect):
				print("\tcreated defect with charge {}".format( charged_defect.charge))
				defects_to_run.append( charged_defect.copy())
				charged_vac_data= open('%s-chg(%s).dat'%(vac_defect.name,charged_defect.charge),'w')
				charged_defect_str=charged_defect.generate_defect_structure( supercell=(2, 1, 1))
				charged_vac_data.writelines(str(charged_defect_str))
				charged_vac_data.close()
				#print(charged_defect.generate_defect_structure( supercell=(2, 1, 1)))


#Generating substitutions and antisites with the bulk structure formula

		defects_to_run=[]
		print('Substitution Generator:')
		sub_element_list=['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl']
		for ele in elements[:]:
			if ele in sub_element_list:
				print(ele)
				sub_element_list.remove(ele)
		print(sub_element_list)

		for element in sub_element_list:
			for sub_defect in SubstitutionGenerator( bulk_structure, element):
				print("\tCreated Defect {} at site {}".format( sub_defect.name, sub_defect.site))
				for charged_defect in SimpleChargeGenerator( sub_defect):
					print("\tcreated defect with charge {}".format( charged_defect.charge))
					defects_to_run.append( charged_defect.copy())
					charged_sub_data= open('%s-chg(%s).dat'%(sub_defect.name,charged_defect.charge),'w')
					charged_defect_str=charged_defect.generate_defect_structure( supercell=(2, 1, 1))
					charged_sub_data.writelines(str(charged_defect_str))
					charged_sub_data.close()
					#print(charged_defect.generate_defect_structure( supercell=(2, 1, 1)))


		defects_to_run=[]
		print('Antisite Generator:')
		for element in elements:
			for sub_defect in SubstitutionGenerator( bulk_structure, element):
				print("\tCreated Defect {} at site {}".format( sub_defect.name, sub_defect.site))
				for charged_defect in SimpleChargeGenerator( sub_defect):
					print("\tcreated defect with charge {}".format( charged_defect.charge))
					defects_to_run.append( charged_defect.copy())
					charged_sub_data= open('antisite-%s-chg(%s).dat'%(sub_defect.name,charged_defect.charge),'w')
					charged_defect_str=charged_defect.generate_defect_structure( supercell=(2, 1, 1))
					charged_sub_data.writelines(str(charged_defect_str))
					charged_sub_data.close()
					#print(charged_defect.generate_defect_structure( supercell=(2, 1, 1)))



#Creating Extrinsic and Intrinsic interstitials with the bulk structure formula
		defects_to_run=[]
		print('Extrinsic Interstitial Generator:')
		sub_element_list=['H', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl']
		for ele in elements[:]:
			if ele in sub_element_list:
				print(ele)
				sub_element_list.remove(ele)
		print(sub_element_list)
		
		for element in sub_element_list:
			for inter_defect in VoronoiInterstitialGenerator( bulk_structure, element):
				print("\tCreated Defect {} at site {}".format( inter_defect.name, inter_defect.site))
				for charged_defect in SimpleChargeGenerator( inter_defect):
					print("\tcreated defect with charge {}".format( charged_defect.charge))
					defects_to_run.append( charged_defect.copy())
					charged_inter_data= open('ext-%s-chg(%s).dat'%(inter_defect.name,charged_defect.charge),'w')
					charged_defect_str=charged_defect.generate_defect_structure( supercell=(2, 1, 1))
					charged_inter_data.writelines(str(charged_defect_str))
					charged_inter_data.close()
					#print(charged_defect.generate_defect_structure( supercell=(2, 1, 1)))


		defects_to_run=[]
		print('Intrinsic Interstitial Generator:')
		for element in elements:
			for inter_defect in VoronoiInterstitialGenerator( bulk_structure, element):
				print("\tCreated Defect {} at site {}".format( inter_defect.name, inter_defect.site))
				for charged_defect in SimpleChargeGenerator( inter_defect):
					print("\tcreated defect with charge {}".format( charged_defect.charge))
					defects_to_run.append( charged_defect.copy())
					charged_inter_data= open('int-%s-chg(%s).dat'%(inter_defect.name,charged_defect.charge),'w')
					charged_defect_str=charged_defect.generate_defect_structure( supercell=(2, 1, 1))
					charged_inter_data.writelines(str(charged_defect_str))
					charged_inter_data.close()
					#print(charged_defect.generate_defect_structure( supercell=(2, 1, 1)))


#Renaming the folders according to their structural formula
		os.rename(os.path.join(os.path.join(basepath,'%s')%(i)),os.path.join(basepath,'%s'+'-'+'%s')%(i,str_formula))
	print(structures)
