from roteasy import align_frame
from math import sqrt

pos_sun=(6,6,0.)
dist_sun=sqrt(pos_sun[0]*pos_sun[0]+pos_sun[1]*pos_sun[1]+pos_sun[2]*pos_sun[2])
x=(0,8.5,pos_sun[0])
y=(0,0,pos_sun[1])
z=(0,0,pos_sun[2])



cord_n=align_frame(cord=(x,y,z),pos_vec=pos_sun,ax='x',xoff=dist_sun,unpacked=True,reference='l', change_reference='x')

print(cord_n)