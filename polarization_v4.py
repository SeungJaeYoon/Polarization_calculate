import numpy as np
import numpy.linalg as lin
import os
import math
import sys
import copy
f_1=open("1_reference/POSCAR","r") 
f_2=open("1_reference/OUTCAR","r")
f_3=open("2_original/POSCAR","r")
f_4=open("2_original/OUTCAR","r")
line_1=f_1.readlines()
line_2=f_2.readlines()
line_3=f_3.readlines()
line_4=f_4.readlines()


def cartesian_transformation(lattice_vector_array,atom_coordinate_array):
	temp_cartesian=[]
	for i in range(3):
		temp_cartesian.append(lattice_vector_array[0][i]*atom_coordinate_array[0]+
			lattice_vector_array[1][i]*atom_coordinate_array[1]+
			lattice_vector_array[2][i]*atom_coordinate_array[2])
	temp_cartesian=np.array(temp_cartesian,dtype=np.float64)
	return temp_cartesian
def read_poscar(x,y): #x는 오픈한 것, y는 라인 넘버
	temp=x
	temp=temp[y-1]
	temp=temp.replace('\n','') #n을 날린다
	temp=temp.split() #공백을 기준으로 리스트 만든거임
	temp=np.array(temp,dtype=np.float64)
	return temp
def read_poscar_coor(x,y):
	temp=x
	temp=temp[y-1]
	temp=temp.replace('\n','') #n을 날린다
	temp=temp.split()
	for i in range(0,3):
		if float(temp[i])>0.85:
			temp[i]=float(temp[i])-1
	temp=np.array(temp,dtype=np.float64)
	return temp
def find_born(x):
	line_num=0
	start=0
	end=0
	for i in x:
		line_num+=1
		if 'BORN EFFECTIVE CHARGES' in i:
			start=line_num
		elif 'INTERNAL STRAIN TENSOR FOR ION    1' in i:
			end=line_num
	return [start,end]
def read_vol(x):
        for line in x:
                if "volume of cell" in line:
                        line=line.split()
                        vol=float(line[4])
                        vol=vol*10**(-30)
        return vol

sys.stdout=open('Polarization.txt','w')
atom_list=line_3[5]
atom_list=atom_list.replace('\n','')
atom_list=atom_list.split()
atom_num_list=read_poscar(line_1,7)
atom_num=int(atom_num_list.sum())
new_atom_list=[]
for i in range(0,len(atom_list)):
	num=atom_num_list[i]
	while num>0:
		num=num-1
		num_1=int(atom_num_list[i]-num)
		new_atom_list.append('%s%s'%(atom_list[i],num_1))
		
reference_coor=[]
original_coor=[]
for i in range(9,9+atom_num): #range는 끝 숫자 포함 x
	reference_coor.append(read_poscar_coor(line_1,i))
	original_coor.append(read_poscar_coor(line_3,i))
reference_coor=np.array(reference_coor)
original_coor=np.array(original_coor)

lattice_vector_reference=[]
lattice_vector_original=[]
for i in range(3,6):
	lattice_vector_reference.append(read_poscar(line_1,i))
	lattice_vector_original.append(read_poscar(line_3,i))
lattice_vector_reference=np.array(lattice_vector_reference)
lattice_vector_reference=lattice_vector_reference.reshape(3,3)
lattice_vector_original=np.array(lattice_vector_original)
lattice_vector_original=lattice_vector_original.reshape(3,3)

displacement_coor_dir=original_coor-reference_coor
displacement_coor=[]


for i in range(0,atom_num):
	x=displacement_coor_dir[i]
	x=cartesian_transformation(lattice_vector_original,x)
	displacement_coor.append(x)
displacement_coor=np.array(displacement_coor)

print('\n Displacement value \n')
for i in range(0,atom_num):
	print('%s %s %s %s'%(new_atom_list[i],displacement_coor[i][0],displacement_coor[i][1],displacement_coor[i][2]))

num_start=find_born(line_2)[0]
num_end=find_born(line_2)[1]
born_effective_per_ion=[]

for i in range(num_start+1,num_end-2):
	born_effective_per_ion.append(line_2[i])
'''
new_born_=copy.deepcopy(born_effective_per_ion)

print('\n Born effective charge of each ions(only xx,yy,zz) \n')
num=0
for i in range(0,len(new_born_)):
	temp_born=new_born_[i]
	temp_born=temp_born.replace('\n','')
	temp_born=temp_born.split()
	temp_born=np.array(temp_born)
	if 'ion' in new_born_[i]:
		temp_born_=[]
	elif (i%4)==1:
		temp_born_x=float(temp_born[1])
		temp_born_.append(temp_born_x)
	elif (i%4)==2:
		temp_born_y=float(temp_born[2])
		temp_born_.append(temp_born_y)
	else:
		temp_born_z=float(temp_born[3])
		temp_born_.append(temp_born_z)
		print('%s %s %s %s'%(new_atom_list[num],temp_born_[0],temp_born_[1],temp_born_[2]))
		num+=1
'''
total_charge=[]
def charge_component(x):
	charge=born_effective_per_ion[x].split()
	charge=[charge[1],charge[2],charge[3]]
	return charge
for i in range(0,atom_num):
	charge_x=charge_component(4*i+1)
	charge_x[1]=0
	charge_x[2]=0
	charge_x=np.array(charge_x,dtype=np.float64)
	charge_y=charge_component(4*i+2)
	charge_y[0]=0
	charge_y[2]=0
	charge_y=np.array(charge_y,dtype=np.float64)
	charge_z=charge_component(4*i+3)
	charge_z[0]=0
	charge_z[1]=0
	charge_z=np.array(charge_z,dtype=np.float64)
	charge_xyz=[charge_x,charge_y,charge_z]
	charge_xyz=np.array(charge_xyz)
	charge_xyz=charge_xyz.reshape(3,3)
	temp_dis=displacement_coor[i]
	temp_dis=temp_dis.reshape(3,1)
	mul_charge=charge_xyz@temp_dis
	mul_charge=mul_charge.reshape(3,)
	total_charge.append(mul_charge)
total_charge=np.array(total_charge)
total_charge=total_charge*1.60217646*(10**(-19))*0.0000000001/read_vol(line_4)
zeros=[0,0,0]
total_charge_sum=np.array(zeros)

print('\n Reference | Polarizations of each atoms \n')
for i in range(0,atom_num):
	print('%s %s %s %s'%(new_atom_list[i],total_charge[i][0],total_charge[i][1],total_charge[i][2]))
for i in total_charge:
	total_charge_sum=total_charge_sum+i


reference_pol_x=total_charge_sum[0]
reference_pol_y=total_charge_sum[1]
reference_pol_z=total_charge_sum[2]
reference_pol=(math.sqrt(reference_pol_x**2+reference_pol_y**2+reference_pol_z**2))


num_start=find_born(line_4)[0]
num_end=find_born(line_4)[1]
born_effective_per_ion=[]

for i in range(num_start+1,num_end-2):
	born_effective_per_ion.append(line_4[i])

new_born_=copy.deepcopy(born_effective_per_ion)
print('\n Born effective charge of each ions(only xx,yy,zz) \n')
num=0
for i in range(0,len(new_born_)):
        temp_born=new_born_[i]
        temp_born=temp_born.replace('\n','')
        temp_born=temp_born.split()
        temp_born=np.array(temp_born)
        if 'ion' in new_born_[i]:
                temp_born_=[]
        elif (i%4)==1:
                temp_born_x=float(temp_born[1])
                temp_born_.append(temp_born_x)
        elif (i%4)==2:
                temp_born_y=float(temp_born[2])
                temp_born_.append(temp_born_y)
        else:
                temp_born_z=float(temp_born[3])
                temp_born_.append(temp_born_z)
                print('%s %s %s %s'%(new_atom_list[num],temp_born_[0],temp_born_[1],temp_born_[2]))
                num+=1
	

total_charge=[]

for i in range(0,atom_num):
	charge_x=charge_component(4*i+1)
	charge_x[1]=0
	charge_x[2]=0
	charge_x=np.array(charge_x,dtype=np.float64)
	charge_y=charge_component(4*i+2)
	charge_y[0]=0
	charge_y[2]=0
	charge_y=np.array(charge_y,dtype=np.float64)
	charge_z=charge_component(4*i+3)
	charge_z[0]=0
	charge_z[1]=0
	charge_z=np.array(charge_z,dtype=np.float64)
	charge_xyz=[charge_x,charge_y,charge_z]
	charge_xyz=np.array(charge_xyz)
	charge_xyz=charge_xyz.reshape(3,3)
	temp_dis=displacement_coor[i]
	temp_dis=temp_dis.reshape(3,1)
	mul_charge=charge_xyz@temp_dis
	mul_charge=mul_charge.reshape(3,)
	total_charge.append(mul_charge)
total_charge=np.array(total_charge)
total_charge=total_charge*1.60217646*(10**(-19))*0.0000000001/read_vol(line_4)
zeros=[0,0,0]
total_charge_sum=np.array(zeros)

print('\n Original | Polarizations of each atoms \n')
for i in range(0,atom_num):
	print('%s %s %s %s'%(new_atom_list[i],total_charge[i][0],total_charge[i][1],total_charge[i][2]))
for i in total_charge:
	total_charge_sum=total_charge_sum+i

original_pol_x=total_charge_sum[0]
original_pol_y=total_charge_sum[1]
original_pol_z=total_charge_sum[2]
original_pol=(math.sqrt(original_pol_x**2+original_pol_y**2+original_pol_z**2))
average_pol_x=(reference_pol_x+original_pol_x)/2
average_pol_y=(reference_pol_y+original_pol_y)/2
average_pol_z=(reference_pol_z+original_pol_z)/2
average_pol=(reference_pol+original_pol)/2

#sys.stdout=open('Polarization.txt','w')
print("\nReference: Px= %s Py= %s Pz= %s"%(reference_pol_x,reference_pol_y,reference_pol_z))

print("Original:  Px= %s Py= %s Pz= %s"%(original_pol_x,original_pol_y,original_pol_z))

print("Average:   Px= %s Py= %s Pz= %s"%(average_pol_x,average_pol_y,average_pol_z))

print("Reference Polarization: %s"%reference_pol)

print("Original Polarization:  %s"%original_pol)

print("Average Polarization:   %s"%average_pol)
sys.stdout.close()
