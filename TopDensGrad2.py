#!/usr.bin/env python

# Analyze density and size distribution for swapping in >>previously grown<<< film
#         The swap region is swept through the film at 0.1 layer per cycle
#         starting at height xx on cycle 1
# Inputs: atomfile describing the set of atoms in this simulation
#         dump file created during growth of deposited film recording time series of atom positions
#         have to grow films with dump file
# Calculations:
# Average density vs distance from growing surface
# Atom type distribution vs distance from growing surface
#
#

# Output: type_dist.txt. density and type distribution vs depth below surface averaged over all data-dumps:

#
# R.B. Stephens August 2022
#
# Written for Python3
import sys, math, os
from numpy import *

try:
  inatom = sys.argv[1]; indump = sys.argv[2];
except:
  print("Usage:", sys.argv[0], "(1) inatom (2) indump ")
  sys.exit(1)
  
# for swap region with height 6, buried 1 below the surface rising by 0.1 on each cycle
sw_ht = 6
sw_bur = 1
del_ht = 0.1

# open files for reading (the data file and the dump file)
iatom = open( inatom, 'r') # contains specs for atom types
idump = open( indump, 'r') # contains positions of all atoms at each dump time

# create tag for output
infile = os.path.basename(inatom)
outtag = infile.split('.')
print (outtag[0])

box_x = zeros(2)
box_y = zeros(2)
box_z = zeros(2)

# read first lines of atomfile
# determine number of types and sizes of each
line = iatom.readline()
if not line:
  print("missing line in atomfile"); exit()
iwords = line.split()
if not iwords[0]=='LAMMPS':
  print("missing atomfile"); exit()
readatm = 1
while readatm:
    line = iatom.readline()
    iwords = line.split()
    if len(iwords)>2:
        if iwords[2]=='types':
          # get simulation dimensions
          ntype = int(iwords[0])+1 # add one; type list starts with 1, not 0
          atm_dia  = zeros(ntype+1)
          atm_vol  = zeros(ntype+1)
          vol_coef = 3.1415926535/6.0
          vol_cor  = sqrt(2) # correcting for distance to potential minimum
        elif iwords[2]=='xlo':
          box_x[0] = float(iwords[0])
          box_x[1] = float(iwords[1])
        elif iwords[2]=='ylo':
          box_y[0] = float(iwords[0])
          box_y[1] = float(iwords[1])
        elif iwords[0]=='PairIJ':
          atm_cnt = 1
          while 1:
            line = iatom.readline()
            iwords = line.split()
            if len(iwords)>0:
               if (int(iwords[0])==atm_cnt)&(int(iwords[1])==atm_cnt):
                 atm_dia[atm_cnt] = float(iwords[3])
                 atm_vol[atm_cnt] = vol_coef*atm_dia[atm_cnt]*atm_dia[atm_cnt]*atm_dia[atm_cnt]*vol_cor
                 atm_cnt += 1
                 if atm_cnt==ntype: readatm = 0; break
iatom.close();
# calculate type number distribution and density when growing surface is over
char_start = 10.0

# range and height bins for density and type distribution
# zero at one above the top of the swap region - increasing with depth
depth_top = 0.0
depth_bot = 9.0
depth_bin = 0.5
bin_vol = depth_bin*(box_x[1]-box_x[0])*(box_y[1]-box_y[0])
base = 3.0
ndepth   = int((depth_bot-depth_top)/depth_bin+0.5)

# read idump line by line
# each block of data is headed by 'ITEM:' followed by an identifier
# and then a fixed number of data lines

# height starts at swap height = 6 and increases by 0.1 on each cycle

# read every atom dump from dump file - starting when film is thicker than depth + base.
# for every atom in every dump increment type distribution array

# repeat read/calculation/output film_top for each dump until the end of file is reached
dump_no = 0
# one line for each atom dump
while 1:
    line = idump.readline()
    if not line: break # jump out of loop at the end of the file
    in_words = line.split()
    if in_words[0] =="ITEM:":
            item_id = in_words[1];
            # "TIMESTEP" is the keyword in the first line of each atom dump- count all TIMESTEPs
            if item_id == "TIMESTEP":
                dump_no += 1;
            # Assume the simulation space is a box; then three pairs of dimensions follow
            elif ((item_id == "BOX") and (dump_no == 1)):
                line = idump.readline()
                if not line: break
                in_words = line.split()
                minx = float(in_words[0]); maxx = float(in_words[1])
                line = idump.readline()
                if not line: break
                in_words = line.split()
                miny = float(in_words[0]); maxy = float(in_words[1])
                line = idump.readline()
                if not line: break
                in_words = line.split()
                minz = float(in_words[0]); maxz = float(in_words[1])
            elif ((item_id == "ATOMS") and (dump_no == 1)):
                # the rest of this keyword line gives variable names - we identify columns containing atom_id, type, z coord
                for i in range(0, len(in_words), 1):
                    if (in_words[i] == 'id'):
                        i_id = i-2
                    if (in_words[i] == 'z'):
                        i_z = i-2
                    if (in_words[i] == 'type'):
                        i_typ = i-2
f_tot = dump_no
print ('Finished %i atom dumps' % (f_tot))
idump.close()
    
first_dump = 80 # start calculating average after swap region moved up 8 (bottom is above base)
last_dump = 320 # stop calculating average when swap region is within 2 of top surface

avg_top = 25 # bounds within which final density average is taken
avg_bot = 10
avg_vol = (avg_top-avg_bot)*(maxx-minx)*(maxy-miny)
avg1_atm_no = zeros(ntype)
avg2_atm_no = zeros(ntype)
avg_count_beg = 0
avg_count_end = 0

# go through dump file again - increment distribution files
idump = open( indump, 'r') # contains positions of all atoms at each dump time
atom_dist = zeros([ndepth+1,ntype+1]) # have to calculate atom type number at each depth bin
dump_sum_no = 0
print ('f_tot = %i' % (f_tot))
for findx in range(f_tot):
    item_id = "NEXTLINE"
    while not(item_id == "TIMESTEP"):
        line = idump.readline()
        if not line: break # jump out of loop at the end of the file
        in_words = line.split()
        if in_words[0] =="ITEM:":
            item_id = in_words[1]
    if (findx > first_dump)and(findx < last_dump): # only increment lists after swap region moved up a distance
        while not(item_id == "NUMBER"): #determine number of atoms in this dump
          line = idump.readline()
          if not line: break # jump out of loop at the end of the file
          in_words = line.split()
          if in_words[0] =="ITEM:":
             item_id = in_words[1]
        line = idump.readline()
        if not line: break
        in_words = line.split()
        atom_no = int(in_words[0])
        while not(item_id == "ATOMS"):
          line = idump.readline()
          if not line: break # jump out of loop at the end of the file
          in_words = line.split()
          if in_words[0] =="ITEM:":
             item_id = in_words[1]
        dens_no = 0
        max_atom_id = 14
        dump_sum_no += 1
        top_height = sw_ht + del_ht*findx + sw_bur # top reference for this dump
        print('analyzing dump_no %i, height = %g' %(findx, top_height))
        for i in range(atom_no):
            line = idump.readline()
            if not line: break
            in_words = line.split()
            at_id  = int(in_words[i_id])
            at_typ = int(in_words[i_typ])
            atom_z = top_height-float(in_words[i_z])
            z_bin = int(atom_z/depth_bin)
            # print('dump[%i]: id = %i, type = %i, height = %g, z_bin = %i' % (findx, at_id, at_typ, atom_z, z_bin))
            if ((z_bin < ndepth) and (z_bin >= 0) and (at_typ <= max_atom_id)):
               atom_dist[z_bin,at_typ] +=1
    # determine average density at beginning of sweep
    if (findx>1)and(findx<12):
        print ('beginning average %i \n' % (findx))
        while not(item_id == "NUMBER"): #determine number of atoms in this dump
          line = idump.readline()
          if not line: break # jump out of loop at the end of the file
          in_words = line.split()
          if in_words[0] =="ITEM:":
             item_id = in_words[1]
        line = idump.readline()
        if not line: break
        in_words = line.split()
        atom_no = int(in_words[0])
        while not(item_id == "ATOMS"):
          line = idump.readline()
          if not line: break # jump out of loop at the end of the file
          in_words = line.split()
          if in_words[0] =="ITEM:":
             item_id = in_words[1]
        avg_count_beg += 1
        for i in range(atom_no):
            line = idump.readline()
            if not line: break
            in_words = line.split()
            at_id  = int(in_words[i_id])
            at_typ = int(in_words[i_typ])
            atom_z = float(in_words[i_z])
            if (atom_z<avg_top)and(atom_z>avg_bot):
               avg1_atm_no[at_typ] += 1
    #determine average density at end of sweep
    if (findx>f_tot-12)and(findx<f_tot-1):
        print ('end average %i \n' % (findx))
        while not(item_id == "NUMBER"): #determine number of atoms in this dump
          line = idump.readline()
          if not line: break # jump out of loop at the end of the file
          in_words = line.split()
          if in_words[0] =="ITEM:":
             item_id = in_words[1]
        line = idump.readline()
        if not line: break
        in_words = line.split()
        atom_no = int(in_words[0])
        while not(item_id == "ATOMS"):
          line = idump.readline()
          if not line: break # jump out of loop at the end of the file
          in_words = line.split()
          if in_words[0] =="ITEM:":
             item_id = in_words[1]
        avg_count_end += 1
        for i in range(atom_no):
            line = idump.readline()
            if not line: break
            in_words = line.split()
            at_id  = int(in_words[i_id])
            at_typ = int(in_words[i_typ])
            atom_z = float(in_words[i_z])
            if (atom_z<avg_top)and(atom_z>avg_bot):
               avg2_atm_no[at_typ] += 1
      
print ('Finished collecting distribution for %i dumps' % (f_tot))

# after collecting info calculate distributions
density_depth = zeros(ndepth+1)
density_value  = zeros(ndepth+1)
type_dist      = zeros([ndepth+1,ntype+1])

# calc density data - normalize
print ('dump_sum_no = %i, bin_vol = %g, atm_vol[4] = %g' % (dump_sum_no, bin_vol, atm_vol[4]))
for i in range(ndepth):
   density_depth[i] = depth_bin*i
   for j in range(ntype):
     density_value[i] += atm_vol[j]*atom_dist[i,j]/dump_sum_no/bin_vol
     
# show density vs depth
for i in range(ndepth):
    print('At film depth[%i], density = %g' % (i, density_value[i]))
  
type_dist_name = ("type_dist.txt")
otype = open(type_dist_name,'w')

#identify data source
out_string = "type distribution from " + outtag[0] + "\n"
otype.write(out_string)

# write first set of labels
out_string = "volume of each type, "
for i in range(ntype):
  next_string = ( '%i, ' % (i))
  out_string = out_string + next_string
out_string = out_string + (" \n")
otype.write(out_string)

out_string = "--, "
for i in range(ntype):
  next_string = ( '%g, ' % (atm_vol[i]))
  out_string = out_string + next_string
out_string = out_string + (" \n")
otype.write(out_string)

# write second set of labels
out_string = "depth, number of given of atom types \n"
otype.write(out_string)

out_string = "--, "
for i in range(ntype):
  next_string = ( '%i, ' % (i))
  out_string = out_string + next_string
out_string = out_string + (", , density \n")
otype.write(out_string)

for i in range(ndepth):
   out_string = ('%g, ' % (density_depth[i]))
   for j in range(ntype):
     next_string = ('%i, ' %(atom_dist[i,j]))
     out_string = out_string + next_string
   out_string = out_string + (', , %g \n' % (density_value[i]))
   otype.write(out_string)

# calculate average density within the usual bounds -- 10 to 25
tot_vol1 = 0
tot_vol2 = 0
for i in range(ntype):
   tot_vol1 += avg1_atm_no[i]*atm_vol[i]
   tot_vol2 += avg2_atm_no[i]*atm_vol[i]
avg_dens1 = tot_vol1/avg_vol/avg_count_beg
avg_dens2 = tot_vol2/avg_vol/avg_count_end

out_string = ('average density from %i to %i =, %g, (start), %g, (end) \n' % (avg_bot, avg_top, avg_dens1, avg_dens2))
otype.write(out_string)
otype.close()
print ('avg_count_beg/end = %i, %i' % (avg_count_beg,avg_count_end))
print ('average density from %i to %i =, %g, (start), %g, (end) \n' % (avg_bot, avg_top, avg_dens1, avg_dens2))
print ('Finished type dist file')
