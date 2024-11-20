#!/usr.bin/env python

# Calculate cluster size distribution for randomly placed 'atoms' in 18x18x35 lattice
# create 18x18x35 array
# randomly add atom to previously unoccupied voxel
# on each addition, check for already occupied neighboring voxel
# identify added atom by cluster: new value or neighbor
# if multiple cluster numbers, choose first in list
# update number of atoms in that cluster, transfer all other atoms to that cluster
# write cluster statistics every 25 cycles after 100
# in a results file - it will be written where-ever terminal is pointed
#
# R.B. Stephens October 2023
#
# Written for Python3
import sys, math, os, datetime
from numpy import *


# print entire array
set_printoptions(threshold=sys.maxsize)

# atoms added
a_no = 0
max_a_no = 2000

# clusters added
c_no = 0
max_c_no = 400
max_c_print = 32

# cluster size distribution
c_size_dist = zeros(100, dtype=int)
# output parameters
min_out_no = 100
out_int = 25

#simulation volume (usable indices are 1-18 and 1-35)
x_dim = 18 # use 1 to 18, p.b.c
y_dim = 18 # use 1 to 18, p.b.c
z_dim = 34 # use 1 to 33, no neighbors at 0 or 34
vox_no = (x_dim)*(y_dim)*(z_dim-1)

# atom cluster id - = 0 if unoccupied
a_list  = zeros((x_dim+1,y_dim+1,z_dim+1), dtype=int)
c_atm_no = zeros(800, dtype=int)
#cluster id of neighbors
cc = zeros(26, dtype=int)
# neighbor cluster  list
na = zeros(9, dtype=int)

# create output file and headers
c_dist_file = "cluster_distrib.txt"
osout = open(c_dist_file, 'w')
osout.write('Cluster size distribution from randomly placed additions \n')
osout.write('Simulation size %i x %i x %i with %i voxels  \n' % (x_dim, y_dim, (z_dim-1), vox_no))
osout.write('atom_tot, clust_tot, max_cluster_size, \n')
osout.write(' , , , , 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, > 30 \n')


while a_no <= max_a_no:
  ran_x = random.randint(1,x_dim+1)
  ran_y = random.randint(1,y_dim+1)
  ran_z = random.randint(1,z_dim)
  if (ran_x*ran_y*ran_z == 0): continue  # can't go out of bounds
  if (a_list[ran_x,ran_y,ran_z] >= 1): continue #can't select already occupied voxel
  a_no += 1
  a_list[ran_x,ran_y,ran_z] = c_no
  # what are neighbor coordinates (p.b.c on x and y axes)
  x0 = (ran_x - 1 + x_dim - 1) % x_dim + 1
  x1 =  ran_x
  x2 = (ran_x + 1 + x_dim - 1) % x_dim + 1
  y0 = (ran_y - 1 + y_dim - 1) % y_dim + 1
  y1 =  ran_y
  y2 = (ran_y + 1 + y_dim - 1) % y_dim + 1
  z0 =  ran_z - 1
  z1 =  ran_z
  z2 =  ran_z + 1
  cc[0] = a_list[x0,y0,z0]
  cc[1] = a_list[x0,y0,z1]
  cc[2] = a_list[x0,y0,z2]
  cc[3] = a_list[x0,y1,z0]
  cc[4] = a_list[x0,y1,z1]
  cc[5] = a_list[x0,y1,z2]
  cc[6] = a_list[x0,y2,z0]
  cc[7] = a_list[x0,y2,z1]
  cc[8] = a_list[x0,y2,z2]
  cc[9] = a_list[x1,y0,z0]
  cc[10] = a_list[x1,y0,z1]
  cc[11] = a_list[x1,y0,z2]
  cc[12] = a_list[x1,y1,z0]
  #cc[4] = a_list[x1,y1,z1]
  cc[13] = a_list[x1,y1,z2]
  cc[14] = a_list[x1,y2,z0]
  cc[15] = a_list[x1,y2,z1]
  cc[16] = a_list[x1,y2,z2]
  cc[17] = a_list[x2,y0,z0]
  cc[18] = a_list[x2,y0,z1]
  cc[19] = a_list[x2,y0,z2]
  cc[20] = a_list[x2,y1,z0]
  cc[21] = a_list[x2,y1,z1]
  cc[22] = a_list[x2,y1,z2]
  cc[23] = a_list[x2,y2,z0]
  cc[24] = a_list[x2,y2,z1]
  cc[25] = a_list[x2,y2,z2]

  nbr_flag = (cc[0]>0) +(cc[1]>0) +(cc[2]>0) +(cc[3]>0) +(cc[4]>0) +(cc[5]>0) +(cc[6]>0) +(cc[7]>0) +(cc[8]>0)\
            +(cc[9]>0) +(cc[10]>0) +(cc[11]>0) +(cc[12]>0) +(cc[13]>0) +(cc[14]>0) +(cc[15]>0) +(cc[16]>0)\
            +(cc[17]>0) +(cc[18]>0) +(cc[19]>0) +(cc[20]>0) +(cc[21]>0) +(cc[22]>0) +(cc[23]>0) +(cc[24]>0) +(cc[25])>0
            
  nbr_id = cc[0] +cc[1] +cc[2] +cc[3] +cc[4] +cc[5] +cc[6] +cc[7] +cc[8] +cc[9] +cc[10] +cc[11] +cc[12]\
            +cc[13] +cc[14] +cc[15] +cc[16] +cc[17] +cc[18] +cc[19] +cc[20] +cc[21] +cc[22] +cc[23] +cc[24] +cc[25]
  if nbr_flag == 0: # no neighbors
    c_no += 1
    a_list[ran_x,ran_y,ran_z] = c_no
    print ('a_no = %i, (%i,%i,%i), c_no = %i, cluster id = %i; no neighbor' % (a_no, ran_x, ran_y, ran_z, c_no, a_list[ran_x,ran_y,ran_z]))
# one neighbor - set cluster id for new atom to that of neighbor
  elif nbr_flag == 1:
    print ('a_no = %i, (%i,%i,%i), c_no = %i, cluster id = %i; one neighbor (%i)' % (a_no, ran_x, ran_y, ran_z, c_no, a_list[ran_x,ran_y,ran_z],nbr_id))
    a_list[ran_x,ran_y,ran_z] = nbr_id
# multi-neighbors - set cluster id for new atom to that of 1st neighbor in list
#                   and reset id of other neighboring clusters to the same
#                   (that will leave one or more clusters numbers containing no atoms)
  elif nbr_flag >= 1:
    print ('a_no = %i, (%i,%i,%i), c_no = %i, cluster id = %i; multi-neighbor (%i)' % (a_no, c_no, a_list[ran_x,ran_y,ran_z],na[0]))
    n_num = 0
    # make list of  the index of neighbors in cc[*]
    for i in range(len(cc)):
      if (cc[i]>0):
        na[n_num] = i
        n_num += 1
    # set new atom cluster id to that of the first cluster in the list
    a_list[ran_x,ran_y,ran_z] =  cc[na[0]]
    # relabel all the atoms in the neighbor clusters by stepping through the entire atom array
    # stop at the end of the neighbor-cluster id  list.
    for i in range(1,len(cc)):
     if (cc[i]>0) == 0: break
     for ix in range(x_dim):
       for iy in range(y_dim):
         for iz in range(z_dim):
           if a_list[ran_x,ran_y,ran_z] == na[i]:
              a_list[ran_x,ran_y,ran_z] = cc[na[0]]
  
  # output cluster data at intervals of 25 atom additions starting at min_out_no ending at 2000
  if ((a_no % out_int) == 0) & (a_no >= min_out_no):
    # zero the cluster list of atom sizes
    for i in range(len(c_atm_no)): c_atm_no[i]=0
    # tally size of each cluster
    for ix in range(x_dim):
       for iy in range(y_dim):
         for iz in range(z_dim):
           index = a_list[ix,iy,iz]
           if index >= len(c_atm_no): index = len(c_atm_no)-1
           c_atm_no[index] += 1
    print ('===========c_atm_no========')
    print (c_atm_no)
    # calculate cluster size distribution
    for i in range(len(c_size_dist)): c_size_dist[i]=0
    top_c_size = 0
    for i in range(1,len(c_atm_no)):
      if c_atm_no[i] == 0: continue # some clusters numbers unused because of mergers
      # track the largest cluster in the simulation
      if c_atm_no[i]>top_c_size: top_c_size = c_atm_no[i]
      if c_atm_no[i]<(len(c_size_dist)):
        indx = c_atm_no[i]
      else: # last index in list counts all cluster sizes larger than array size
        indx = len(c_size_dist)-1
      c_size_dist[indx] += 1
    c_size_dist[(max_c_print-1)]=sum(c_size_dist[(max_c_print-1):])
    print ('===========c_size_dist[:max_c] ============= ')
    print (c_size_dist)
    osout.write(' %i, %i, %i ,'% (a_no, c_no, top_c_size))
    # output clusters size distribution
    # <<<< max cluster size is a little larger than total number of voxels
    for i in range(max_c_print):
      print ('there are %i cluster(s) of size %i' % (c_size_dist[i],i))
      osout.write(', %i' % (c_size_dist[i]))
    osout.write(' \n')
osout.close()
print ('++++++++++++++++++++++')
print ('Simulation Finished ')
print ('++++++++++++++++++++++')
