#!/usr.bin/env python
# note calcdisp5 can't be used for any data series before qb (202110xx) - when the temp/dens and dump outputs were synchronized
# This version is meant for TT series for which the film is held at a temp above Tg to allow transformation
#
# Inputs: atomfile describing the set of atoms in this simulation
#         dump file created during thermocycle of deposited film recording time series of atom positions
#         temp file recording the addatom temperature at each step
#         (the dump and temp files steps must be the same - except the dump file starts at step 0)
#
# Average displacement of atoms in offset files, bin by height
# Calculate Density of atoms within a specified height range
# Collate temperature at each step with density and mobility
#
# Output: mobility.txt: scratch pad for organizing data from dump-file
# Output: mobility_vs_temp.png: contour plot of displacement vs height in simulation vs time step
#           the mobility data in the time-step range specified in the previous graph
#           is binned by temperature and plotted over a specified temperature range
# Output: out_tmd.txt. For each data-dump: frame-no, time-step, temperature of add-atoms averaged over
#           every time step between data-dumps, density and log(mobility) averaged over the middle of the film
#
# R.B. Stephens Jun 2022
#
# Written for Python3
import sys, math
from numpy import *

try:
  inatom = sys.argv[1]; indump = sys.argv[2]; intemp = sys.argv[3];
except:
  print("Usage:", sys.argv[0], "(1) indump (2) inatom (3) intemp")
  sys.exit(1)

iatom = open( inatom, 'r') # open files for reading
idump = open( indump, 'r')
scratchfile  = "mobility.txt"
oscratch = open(scratchfile, 'w') # create & open file for writing into scratchpad

box_x = zeros(2)
box_y = zeros(2)
box_z = zeros(2)

# needed if atoms have been lost - then max index no is larger than total atom no.
# this allows for 200 lost atoms
evap_no = 200

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
          at_type_no = int(iwords[0])
          atm_dia  = zeros(at_type_no+1)
          atm_vol  = zeros(at_type_no+1)
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
                 if atm_cnt==at_type_no: readatm = 0; break
iatom.close();

# read idump line by line
# each block of data is headed by 'ITEM:' followed by an identifier
# and then a fixed number of data lines

# max count - can limit length of cycle that is graphed
# there are 10 frames for the initial hold
# at ramp rate 65, from 0.27 to 0.21 is reached in 78 frames
#                  from 0.26                    in 65
#                  from 0.25                    in 52
#                  from 0.24                    in 39
#                  from 0.23                    in 26
#                  from 0.22                    in 13
#                  from 0.21                    in 0
#                  from 0.20                    in 13
# and from 0.21 to 0.30, in 117 frames
# and from 0.21 to 0.31, in 130 frames
# and from 0.21 to 0.32, in 143 frames
# and from 0.21 to 0.33  in 156 frames
# and frin 0.21 ti 0.34  in 169 frames
# and from 0.21 to 0.35  in 182 frames
# and from 0.21 to 0.36  in 195 frames
# the first 10 frames are not usable because mobility depends on offset
# there are 1000 frames for the high temp hold
# so one initial hold and initial cool, one heat, one hold, one cool
# first ramp up   88 to 205  then hold until 1205 and down to 1322 starting from 0.27 with max 0.30
#                 88 to 231                  1231             1374               0.27 with max 0.32
#                 88 to 257                  1257             1426               0.27 with max 0.34
#                 16 to 185    485 785       1185             1360               0.215 w/  max 0.34  
#                 10 to 127    427 727       1127             1244               0.21 with max 0.30
#                 10 to 140    440 740       1140             1270               0.21 with max 0.31
#                 10 to 153    453 753       1153             1296               0.21 with max 0.32
#                 10 to 166    466 766       1166             1322               0.21 with max 0.33
#                 10 to 179    479 779       1179             1348               0.21 with max 0.34
#                 10 to 192    492 792       1192             1374               0.21 with max 0.35
#                 10 to 205    505 792       1205             1395               0.21 with max 0.36
#                 29 to 198    498 798       1198             1367               0.225 w/  max 0.36
#                   
graph_min_step = 198
graph_max_step = 1198
# uncomment plt.xlim for this to have effect.

# average mobility transition during high temp hold - the front moves downward at a constant velocity
# pick the middle third of the high temp hold
avg_start_step = graph_min_step + 300
avg_length = 300
avg_stop_step = avg_start_step + avg_length
# down motion per time step (fraction of box height per dump step)
front_vel = -0.00001
# put profile into mob_profile
mob_profile = zeros(avg_length)
mob_n = zeros(avg_length)
mobn_profile = zeros(avg_length)

# assume simulation box x,y width = 10 (could read from file)
# also assume apparent motion more than 7 is caused by p.b.c
width = 10.0
ch_width = width/2.0

# present pixels as average(log(mobility)) or log(average(mobility))
# geometric or arithmetic average of mobility
arith_flag = 1 # =1 for arithmetic =0 for geometric
# lower limit on motion
disp_min = 0.1

# bin atom heights into sections
ht_div = 40
# define mobility using position difference on separate frames
offset = 4
# accumulate frame numbers for calculating mobility start at 0 - go to size of offset
f_no = -1
# frame number given by indx
indx = 0
# initialize final arrays
tot_dist = zeros(ht_div)
tot_num  = zeros(ht_div)
ht_label = zeros(ht_div)
i_id = 0
i_xs = 0
i_ys = 0
i_zs = 0
i_typ = 0
# create two lists for each frame
# 1) average offset as a fn of binned height
# 2) number of atome in each height bin
# also density for the frame
# and frame number, time_step
# and get temperature from the tempdens file
#
# create labels for scratchpad file
for i in range(0, ht_div, 1):
    ht_label[i] = (i+0.496)/float(ht_div)
if arith_flag==1:
  oscratch.write('Average density, log(average(displacement)), and atom no. as a function of binned height.\n')
  oscratch.write('Dump_step, Time_step, Time_step, Temp, Density, Density no., Mob at height =  ')
else:
  oscratch.write('Average density, average(log(displacement)), and atom no. as a function of binned height.\n')
  oscratch.write('Dump_step, Time_step, Time_step, Temp, Mob at height = ,')
for i in range(0, ht_div, 1):
   oscratch.write('%4.2g , ' %(ht_label[i]))
oscratch.write(' -- , Atom no at height = ,')
for i in range(0, ht_div, 1):
   oscratch.write(' %5.2g ,' %(ht_label[i]))
oscratch.write('\n')

# open time-step, temperature file - verify type and measure size
itemp = open( intemp, 'r')
# read first lines
# determine number of types and sizes of each
line = itemp.readline()
if not line:
  print("missing line in time-temp file"); exit()
iwords = line.split(',')
if not iwords[0]=="steptempdens":
  print("missing time-temp file"); exit()
# find columns for each variable
istp  = -1
itmp  = -1
idns  = -1
iatno = -1
itmob = -1
ittop = -1
itbot = -1
# read second line - with column headings
line = itemp.readline()
iwords = line.split()
if (len(iwords) > 1):
 for wndx in range(len(iwords)):
   text = iwords[wndx]
   if text == "Step":
     istp = wndx
   elif text == "Temp":
     itmp = wndx
   elif text == "Dens":
     idns = wndx
   elif text == "AtomNo":
     iatno = wndx
   elif text == "Temp_mob":
     itmob = wndx
   elif text == "Temp_top":
     ittop = wndx
   elif text == "Temp_bot":
     itbot = wndx
if (istp<0 or itmp<0 or idns<0 or iatno<0 or itmob<0 or ittop<0 or itbot<0):
    print("missing column in time-temp file"); exit()
fndx = 0
while itemp.readline():
  fndx += 1
tfilelineno = fndx

mob_step = zeros(tfilelineno+7)
mob_temp = zeros(tfilelineno+7)
mob_dens = zeros(tfilelineno+7)
itemp.close()

# open time-step, temperature file - read out data
itemp = open( intemp, 'r')
#skip over file title and column headers
line = itemp.readline()
line = itemp.readline()

# read step-temp-dens data into array
mob_step[0] = 0
mob_temp[0] = 0.0
mob_dens[0] = 0.0
rndx = 1
for rndx in range(1, tfilelineno, 1):
  line   = itemp.readline()
  iwords = line.split()

  mob_step[rndx] = iwords[istp]
  mob_temp[rndx] = iwords[itmp]
  mob_dens[rndx] = iwords[idns]
itemp.close()

# repeat read/calculation/output for each frame until the end of file is reached
# one line for each frame
while 1:
    line = idump.readline()
    if not line: break # jump out of loop at the end of the file
    in_words = line.split()
    if in_words[0] =="ITEM:":
            item_id = in_words[1];
            # "TIMESTEP" is the keyword in the first line of each frame
            if item_id == "TIMESTEP":
                indx += 1
                f_no += 1;
                # only keep 'offset' number of frames
                # with current frame no = offset
                if (f_no >= offset):
                    f_no = offset;
                line = idump.readline()
                if not line: break
                in_words = line.split()
                step_no = int(in_words[0])
            # "NUMBER" = the number of atoms; therefore the number of lines in the ATOMS block
            elif item_id == "NUMBER":
                line = idump.readline()
                if not line: break
                in_words = line.split()
                atom_no = int(in_words[0])
                padded_no = atom_no+evap_no
                pos = zeros(padded_no)
                reg_vol = 0 # zero volume of atoms inside density measuring region
                # create array to hold atom position info for offset number of frames
                # then the atom_id is needed to ensure displacements are for same atom
                if f_no == 0:
                    fpos  = zeros([offset+2,padded_no,3])
                # move atom positions from previous frames down one in array
                for i in range(offset):
                    for j in range(padded_no):
                        for k in range(3):
                            fpos[i,j,k] = fpos[i+1,j,k]
                # initialize output arrays
                tot_num = zeros(ht_div)
                tot_dist = zeros(ht_div)
            # Assume the simulation space is a box; then three pairs of dimensions follow
            elif item_id == "BOX":
                line = idump.readline()
                if not line: break
                in_words = line.split()
                minx = float(in_words[0]); maxx = float(in_words[1])
                htx  = maxx - minx
                line = idump.readline()
                if not line: break
                in_words = line.split()
                miny = float(in_words[0]); maxy = float(in_words[1])
                hty  = maxy - miny
                line = idump.readline()
                if not line: break
                in_words = line.split()
                minz = float(in_words[0]); maxz = float(in_words[1])
                htz  = maxz - minz
            elif item_id == "ATOMS":
                # the rest of this keyword line gives variable names - we identify columns containing atom_id and x,y,z coords
                for i in range(0, len(in_words), 1):
                    if (in_words[i] == 'id'): # some ids may be missing - atoms evaporated
                        i_id = i-2
                    if (in_words[i] == 'xs'):
                        i_xs = i-2
                    if (in_words[i] == 'ys'):
                        i_ys = i-2
                    if (in_words[i] == 'zs'):
                        i_zs = i-2
                    if (in_words[i] == 'type'):
                        i_typ = i-2
                # put list of id, type, and position of each atom into list for current frame
                max_atom_id = 0
                for i in range(atom_no):
                    line = idump.readline()
                    if not line: break
                    in_words = line.split()
                    at_id  = int(in_words[i_id])
                    if (at_id>max_atom_id):
                      max_atom_id = at_id
                    at_typ = int(in_words[i_typ])
                    x_cpos = float(in_words[i_xs])*htx
                    y_cpos = float(in_words[i_ys])*hty
                    z_cpos = float(in_words[i_zs])*htz
                    at_pos = [x_cpos,y_cpos,z_cpos]
                    fpos[f_no][at_id] = at_pos
                # calculate atom displacement (the coordinates xs, ys, zs are scaled to simulation box size
                this_frame = fpos[0]
                prev_frame = fpos[offset]
                dispx = zeros((padded_no))
                dispy = zeros((padded_no))
                dispz = zeros((padded_no))
                disp  = zeros((padded_no,2))
                offsetx = 0.0
                offsety = 0.0
                for iatm in range(padded_no):
                    dispx[iatm] = this_frame[iatm][0] - prev_frame[iatm][0]
                    dispy[iatm] = this_frame[iatm][1] - prev_frame[iatm][1]
                    dispz[iatm] = this_frame[iatm][2] - prev_frame[iatm][2]
                    # correct displacements for p.b.c errors
                    if abs(dispx[iatm])>ch_width:
                        if dispx[iatm]>0.0:
                            dispx[iatm] = dispx[iatm]-width
                        else:
                            dispx[iatm] = dispx[iatm]+width
                    if abs(dispy[iatm])>ch_width:
                        if dispy[iatm]>0.0:
                            dispy[iatm] = dispy[iatm]-width
                        else:
                            dispy[iatm] = dispy[iatm]+width
                # correct displacement for mass motion
                offsetx = sum(dispx)/len(dispx)
                offsety = sum(dispy)/len(dispy)
                for iatm in range(padded_no):
                    if iatm==0:continue
                    disp[iatm,0] = sqrt((dispx[iatm]-offsetx)**2 + (dispy[iatm]-offsety)**2 + dispz[iatm]**2)
                    disp[iatm,1] = (this_frame[iatm][2] + prev_frame[iatm][2])/2.0
                # create array with average displacement binned by height
                for iatm in range(padded_no): # sum displacements (or log displacements) at each binned height
                    h_indx = int(disp[iatm,1]/htz*float(ht_div)) # convert calculation to integer - ensure >0
                    if (h_indx >= ht_div):
                      h_indx = ht_div-1
                    tot_num[h_indx] += 1
                    if arith_flag==1:
                      tot_dist[h_indx] += abs(disp[iatm,0])
                    else:
                      if abs(disp[iatm,0])>disp_min:
                        tot_dist[h_indx] += log(abs(disp[iatm,0]))
                      else:
                        tot_dist[h_indx] += log(disp_min)
                        
                # output frame no, time step for this frame
                if indx==0:
                   out_string = (' %i , %i, %i, %g, ' %(indx, step_no, mob_step[0], mob_temp[1]))
                elif indx-1 < len(mob_temp):
                   out_string = (' %i , %i, %i, %g, ' %(indx, step_no, mob_step[indx-1], mob_temp[indx-1]))
                else :
                   out_string = (' %i , %i, nan, nan ' %(indx, step_no))
                # output log(average(mobility)) or average(log(mobility) for each bin
                #   (summed displacement of all relevant atoms divided by their number)
                for i in range(ht_div):
                    mob = ' nan ,'
                    if tot_num[i]>0:
                        val = tot_dist[i]
                        mob = (' % g ,' % (val))
                    out_string = out_string + mob
                out_string = out_string + ' -- ,'
                # output atom number in each bin
                for i in range(ht_div):
                    at_no = (' % i ,' % (tot_num[i]))
                    out_string = out_string + at_no
                out_string = out_string + ' \n'
                oscratch.write(out_string)
                out_string = ('frame %i ' % (indx))
                print (out_string)
            else: break
f_tot = indx
print ('Finished %i frames' % (indx))
idump.close(); oscratch.close()

# graph the mobility data as contour plot of height in simulation box vs frame no
# overlay temperature data
import matplotlib.pyplot as plt
from matplotlib import cm

iscratch = open(scratchfile, 'r') # read entire array and create 3-d graph

# read titles
line = iscratch.readline()
# read data labels
line = iscratch.readline()

# data resolution
data_time_no = f_tot
data_height_no = ht_div

delta = 1.0/float(data_height_no)

# arrays used to calculate value to put into each graph point
z = zeros((data_height_no,data_time_no), dtype=float)
n = zeros((data_height_no,data_time_no), dtype=float)

flm_dens    = zeros(data_time_no)
flm_dens_no = zeros(data_time_no)
time_step   = zeros(data_time_no)
frame_no    = zeros(data_time_no)
flm_temp    = zeros(data_time_no)
        
# read log(mobility) vs height bin and time step from scratchpad (delimited by commas)
for i in range(data_time_no):
    line = iscratch.readline()
    if not line:
        print('unexpected end of file')
        exit()
    in_words = line.split(',')
    flm_temp[i] = float(in_words[3])
    out_string = ('flm_temp[%i] = %g' % (i, flm_temp[i]))
    print(out_string)
    time_step[i] = in_words[1]
    frame_no[i] = in_words[0]
    for j in range(data_height_no):
        z[j,i] = float(in_words[j+4])
        n[j,i] = float(in_words[j+data_height_no+5])
        if (i - avg_start_step) in range(avg_length):
           ki = i - avg_start_step
           kj = int(j - ki*front_vel)
           if (kj > 0):
              mob_profile[kj] = + z[j,i]
              mob_n[kj] = + n[j,i]

# close input file
iscratch.close()

# normalize mobility profile
for i in range(data_height_no):
    mobn_profile[i] = nan
    if mob_n[i] > 0:
       mobn_profile[i] = log(mob_profile[i]/mob_n[i])
       
# output mobility profile
profile_file = "mobility_profile.txt"
oprofile = open(profile_file, 'w')
oprofile.write('Mobility vs height across the melt front \n')
oprofile.write('index, height, number, total motion, log_mobility \n')
for i in range (data_height_no):
   prof_height = i*htz/ht_div
   oprofile.write('%i, %g, %g, %g, %g \n' % (i,prof_height, mob_n[i], mob_profile[i], mobn_profile[i]))

print ("Starting mobility & temp vs time-step plot")

# graph resolution
bin_step = 3 # no of time-step observations to average
graph_x_no = int(data_time_no/bin_step + 2)
graph_y_no = ht_div
graph_min_bin_step = int(graph_min_step/bin_step)
graph_max_bin_step = int(graph_max_step/bin_step)

min_step = 0
max_step = graph_x_no
bin_step = 3 # no of observations to average
max_x_axis = max_step/bin_step
z_mob = zeros((graph_y_no,graph_x_no), dtype=float)
z_temp = zeros((graph_y_no,graph_x_no), dtype=float)
xx = arange(min_step, max_step)
yy = arange(0.0, 1.0, delta)
xx_1, yy_1 = meshgrid(xx,yy)
# bin temps by time step
# sum mobilities, atom nos, and temps
tot_disp = zeros((graph_y_no,graph_x_no), dtype=float)
tot_no = zeros((graph_y_no,graph_x_no), dtype=float)
avg_temp = zeros((graph_x_no), dtype=float)

stop = data_time_no
for i in range(stop):
  bin_indx = int(i/bin_step)
  for j in range(graph_y_no):
      tot_disp[j,bin_indx] += z[j,i]/bin_step
      tot_no[j,bin_indx] += n[j,i]/bin_step
      
# calculate average
print('calculate average \n')
print('data_time_no = %i, graph_x_no = %i \n' % (data_time_no, graph_x_no))
for i in range(graph_x_no):
   for j in range(graph_y_no):
    z_mob[j,i] = nan
    if tot_no[j,i]>0:
      if arith_flag == 1:
        z_mob[j,i]=log(tot_disp[j,i]/tot_no[j,i])
        if z_mob[j,i] < -3: z_mob[j,i] = -3
      else:
        z_mob[j,i]=tot_disp[j,i]/tot_no[j,i]
        if z_mob[j,i] < -3: z_mob[j,i] = -3

fig = plt.figure(figsize=(13,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
# comment this out to get temp span matched to data
plt.xlim(graph_min_bin_step,graph_max_bin_step)

# define plot
levels = list(range(-3,2,1))
contour = plt.contour(xx_1, yy_1, z_mob, levels, colors='darkgrey')
plt.clabel(contour,colors = 'k',fmt = '%2.1f', fontsize=12)
cs = plt.contourf(xx_1, yy_1, z_mob, alpha=1.0, cmap='inferno', vmax=2, vmin=-3, levels=25)
plt.colorbar(cs)

# plot label
if arith_flag==1:
  plt.title('Thin film log(average(mobility)) plot during thermocycle')
else:
  plt.title('Thin film average(log(mobility)) plot during thermocycle')
file_words = indump.split('/')
name_length = len(file_words)
print_name = file_words[name_length-1]
print (print_name)
text_kwargs = dict(ha='left', va='center', fontsize=12, color='black', backgroundcolor='white')
x_loc = min_step + (max_step-min_step)/40
plt.text(x_loc,0.9,print_name,**text_kwargs)
plt.xlabel('time steps')
plt.ylabel('height in simulation box -- temperature')

plt.savefig('mobility_vs_time.png', dpi=500)
plt.show()
plt.close()

print ("Starting mobility vs time-step plot")
# data from within min_step, max_step
# temp value of x-axis shown by modifying appropriate mobility pixel
min_temp = 0.15
max_temp = 0.35
