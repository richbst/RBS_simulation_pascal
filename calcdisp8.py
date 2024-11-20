#!/usr.bin/env python
# note calcdisp8 can't be used for any data series before qb (202110xx) - when the temp/dens and dump outputs were synchronized
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
#           one graph for initial heat up
#           one graph for final cool-down
# Output: out_tmd.txt. For each data-dump: frame-no, time-step, temperature of add-atoms averaged over
#           every time step between data-dumps, density and log(mobility) averaged over the middle of the film
#
# R.B. Stephens Dec 2021
# replace displacement with square of displacement to more nearly be diffusivity
# R.B. Stephens Oct 2023
#
# R.B. Stephens July 2024
# Written for Python3
import sys, math
from numpy import *

try:
  inatom = sys.argv[1]; indump = sys.argv[2]; intemp = sys.argv[3];
except:
  print("Usage:", sys.argv[0], "(1) inatom (2) indump (3) intemp")
  sys.exit(1)

print("start")
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
# calculate average mobility and density between height limits (assume film goes to 36)
mob_top = 30.0
mob_bot = 10.0

print('read atom file, type no = %i' % (at_type_no))
# density averaging volume - initial volume only - can change to maintain 0 lateral pressure
den_vol = (box_x[1]-box_x[0])*(box_y[1]-box_y[0])*(mob_top-mob_bot)

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
#                  from 0.127                   in 108
#                  from 0.125                   in 110
#                  from 0.123                   in 113
#                  from 0.120                   in 117
# and from 0.21 to 0.41, in 260 frames
# and from 0.21 to 0.45, in 312 frames
# and from 0.21 to 0.47, in 338 frames
# and from 0.21 to 0.48  in 351 frames
# the first 10 frames are not usable because mobility depends on offset
# there are 10 frames for the high temp hold
# so one initial hold and initial cool, two heats, two cools, two holds
# give first ramp from  88 to 348 and last cooldown from 888 to 1148 frames starting from 0.27 (early 0.27 runs didn't have the initial hold - so 878 to 1138)
#                   88 to 348
#                   75 to 335        warming         cooling
#                   62 to 322         592  to  862 to 1122 frames starting from 0.25 with max 0.41
#                   49 to 309         579  to  849 to 1109                 from 0.24 with max 0.41
#                   36 to 296         566  to  836 to 1096                 from 0.23 with max 0.41
#                   23 to 283         553  to  823 to 1083                 from 0.22 with max 0.41
#
#                   36 to 348         670  to  992 to 1304                 from 0.23 with max 0.45
#                   10 to 322         644  to  966 to 1278                 from 0.21 with max 0.45
#
#                   88 to 426         774  to 1121 to 1460                 from 0.27 with max 0.47
#                   81 to 419         767  to 1115 to 1453                 from 0,265 w/  max 0.47
#                   75 to 413         761  to 1109 to 1447                 from 0.26 with max 0.47
#                   68 to 406         754  to 1102 to 1440                 from 0.255 w/  max 0.47
#                   62 to 400         748  to 1096 to 1434                 from 0.25 with max 0.47
#                   55 to 393         741  to 1089 to 1427                 from 0.245 w/  max 0.47
#                   49 to 387         735  to 1083 to 1421                 from 0.24 with max 0.47
#                   42 to 380         728  to 1076 to 1414                 from 0.235 w/  max 0.47
#                   36 to 374         722  to 1070 to 1408                 from 0.23 with max 0.47
#                   29 to 367         715  to 1063 to 1401                 from 0.225 w/  max 0.47
#                   23 to 361         709  to 1057 to 1395                 from 0.22 with max 0.47
#                   10 to 348         696  to 1044 to 1382                 from 0.21 with max 0.47
#                   23 to 361         709  to 1057 to 1395                 from 0.20 with max 0.47
#
#
#                   36 to 387         750  to 1105 to 1460                 from 0.23 with max 0.48
#                   30 to 380         744  to 1106 to 1456                 from 0.225 w/  max 0.48
#                   23 to 374         740  to 1101 to 1447                 from 0.22 with max 0.48
#                   16 to 368         xxx  to xxxx to xxxx                 from 0.215 w/  max 0.48
#                   10 to 361         698  to 1057 to 1395                 from 0.21 with max 0.48
#
#                   127 to 463        801  to 1160 to 1498                 from 0.120 with max 0.48
# range for time_step graph

dep_temp = 0.120
start_cyc_temp = 0.21
fin_cyc_temp = 0.47
rest_steps= 10
steps_per_T = 1300
min_indx = int(rest_steps + abs(start_cyc_temp - dep_temp)*steps_per_T)
max_indx = int(min_indx + (fin_cyc_temp - start_cyc_temp)*steps_per_T)
max_f_indx = int(max_indx + rest_steps + 2*(fin_cyc_temp - start_cyc_temp)*steps_per_T)
min_f_indx = int(max_f_indx + (fin_cyc_temp - start_cyc_temp)*steps_per_T)
# max_indx = 450
# min_indx = 110
# second heating
# max_f_indx = 798
# min_f_indx = 1134
bin_indx = 1

print ('dep_temp = %f, start_cyc_temp = %f, fin_cyc_temp = %f ' % (dep_temp, start_cyc_temp, fin_cyc_temp ))
print ('rest_steps = %i, steps_per_T = %i' % (rest_steps, steps_per_T))
print ('calc: min_indx = %i, max_indx = %i, max_f_indx = %i, min_f_indx = %i' % (min_indx, max_indx, max_f_indx, min_f_indx))
# print ('prev: min_indx = %i, max_indx = %i, max_f_indx = %i, min_f_indx = %i' % (min_indx, max_indx, max_f_indx, min_f_indx))
# uncomment plt.xlim for this to force chart limits rather than match data.

# range for display in temperature graph
max_temp = 0.5
min_temp = 0.2
bin_temp = 0.002
ntemps   = int((max_temp-min_temp)/bin_temp+0.5)
# assume simulation box x,y width = 10 (could read from file)
# also assume apparent motion more than 7 is caused by p.b.c
width = 10.0
ch_width = width/2.0

# present pixes as average(log(mobility)) or log(average(mobility))
# geometric or arithmetic average of mobility
arith_flag = 1 # =1 for arithmetic =0 for geometric
# lower limit on motion
disp_min = 0.1

# bin atom heights into sections
ht_div = 20
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
  oscratch.write('Dump_step, Time_step, Time_step, Temp, Density, Density no., Mob at height =  ')
for i in range(0, ht_div, 1):
   oscratch.write('%4.2g , ' %(ht_label[i]))
oscratch.write(' -- , Atom no at height =')
for i in range(0, ht_div, 1):
   oscratch.write('%5.2g , ' %(ht_label[i]))
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

# create temp-mob-dens file
out_tmd = open("temp-mob-dens.txt", 'w')
# create labels for mob-dens file
if arith_flag==1:
  out_tmd.write('Average temperature, atom-no density, and log(average(displacement^2)) at each time step. \n')
  out_tmd.write('Atom_No, density, and displacement in the height range between %5.2f and %5.2f \n' %(mob_bot,mob_top))
  out_tmd.write('dump_step, time_step, mob_temp, mob_dens, dump_density, dump_density-no, log(average(displacement)) \n')
else:
  out_tmd.write('Average temperature, atom-no density, and average(log(displacement^2)) at each time step. \n')
  out_tmd.write('Atom_No, density, and displacement in the height range between %5.2f and %5.2f \n' %(mob_bot,mob_top))
  out_tmd.write('dump_step, time_step, mob_temp, mob_dens, dump_density, dump_density-no, average(log(displacement)) \n')

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
                col_no = len(in_words)-2
                # put list of id, type, and position of each atom into list for current frame
                dens_no = 0
                max_atom_id = 0
                for i in range(atom_no):
                    line = idump.readline()
                    if not line: break
                    in_words = line.split()
                    if len(in_words)<col_no: break
                    at_id  = int(in_words[i_id])
                    if (at_id>max_atom_id):
                      max_atom_id = at_id
                    at_typ = int(in_words[i_typ])
                    x_cpos = float(in_words[i_xs])*htx
                    y_cpos = float(in_words[i_ys])*hty
                    z_cpos = float(in_words[i_zs])*htz
                    at_pos = [x_cpos,y_cpos,z_cpos]
                    fpos[f_no][at_id] = at_pos
                    # do volume sum for density
                    if (z_cpos >= mob_bot)& (z_cpos <= mob_top):
                       reg_vol += atm_vol[at_typ]
                       dens_no += 1
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
                for iatm in range(padded_no): # sum displacements^2 (or log displacements^2) at each binned height
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
                        
                # output frame no, time step, density, and no of atoms used in density calc for this frame
                # use htx, hty because box size can change from frame to frame
                dens = reg_vol/(htx*hty*(mob_top-mob_bot))
                if indx==0:
                   out_string = (' %i , %i, %i, %g, %g, %i, ' %(indx, step_no, mob_step[0], mob_temp[1], dens, dens_no))
                elif indx-1 < len(mob_temp):
                   out_string = (' %i , %i, %i, %g, %g, %i, ' %(indx, step_no, mob_step[indx-1], mob_temp[indx-1], dens, dens_no))
                else :
                   out_string = (' %i , %i, nan, nan, %g, %i, ' %(indx, step_no, dens, dens_no))
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
import matplotlib.pyplot as plt
from matplotlib import cm

iscratch = open(scratchfile, 'r') # read entire array and create 3-d graph

# read titles
line = iscratch.readline()
# read data labels
line = iscratch.readline()

nrows = f_tot
ncols = ht_div

delta = 1.0/float(ncols)

z = zeros((ncols,nrows), dtype=float)
n = zeros((ncols,nrows), dtype=float)
if f_tot>max_indx:
  array_no = f_tot
else:
  array_no = max_indx
flm_dens    = zeros(array_no)
flm_dens_no = zeros(array_no)
time_step   = zeros(array_no)
frame_no    = zeros(array_no)
flm_temp    = zeros(array_no)

print('array_no = %i, f_tot = %i, tfilelineno = %i'%(array_no, f_tot, tfilelineno))
# read log(mobility) vs height bin and time step from scratchpad (delimited by commas)
if (tfilelineno<array_no): array_no = tfilelineno
for i in range(array_no):
    line = iscratch.readline()
    if not line:
        print('unexpected end of file')
        exit()
    in_words = line.split(',')
    flm_dens_no[i] = in_words[5]
    flm_dens[i] = in_words[4]
    flm_temp[i] = float(in_words[3])
    time_step[i] = in_words[1]
    frame_no[i] = in_words[0]
    for j in range(ncols):
        z[j,i] = float(in_words[j+6])
    for j in range(ncols):
        n[j,i] = float(in_words[j+ncols+7])

# close input file
iscratch.close()

# convert height limits to indices
ntop = int(mob_top/htz*float(ht_div))
nbot = int(mob_bot/htz*float(ht_div))

if (mob_step[3]!=time_step[3]):
  print('mob/filmtop time-steps not synchronized')
  exit()

for i in range(array_no):
  mob_avg = nan
  # get average mobility values within index limits for each time step
  mob_sum = 0
  mob_num = 0
  if i<(array_no-15):
    for j in range(nbot,ntop,1):
      mob_sum += z[j,i]
      mob_num += n[j,i]
    if mob_num > 0:
      if arith_flag == 1:
        mob_avg = log(mob_sum/mob_num)
      else:
        mob_avg = mob_sum/mob_num
    # write data line to output - include density calculated from dump file and step, temp, dens from std file
    out_string = ('%i, %i, %6.4f, %6.4f, %6.4f, %i, %6.4f \n' % (frame_no[i], time_step[i], mob_temp[i], mob_dens[i], flm_dens[i], flm_dens_no[i], mob_avg))
    out_tmd.write(out_string)
out_tmd.close()

print ("Starting mobility vs temp plot - 1")
# create plot matrix
# data from within min_indx, max_indx
# temp range of x-axis can be separately defined
zz = zeros((ncols,ntemps), dtype=float)
xx = arange(min_temp, max_temp, bin_temp)
yy = arange(0.0, 1.0, delta)
xx_1, yy_1 = meshgrid(xx,yy)
# bin time_steps by temp
# sum mobilities and atom nos
tot_disp = zeros((ncols,ntemps), dtype=float)
tot_no = zeros((ncols,ntemps), dtype=int)
for i in range(min_indx, min(f_tot, max_indx), 1):
  tndx = int(ntemps*(flm_temp[i]-min_temp)/(max_temp-min_temp)+0.5)
  bin_str = ( 'indx = %i tndx = %g, film_temp = %g, max_temp = %g, min_temp = %g' % (i, tndx, flm_temp[i], max_temp, min_temp))
  print (bin_str)
  if tndx>=ntemps: continue
  for j in range(ncols):
      tot_disp[j,tndx] += z[j,i]
      tot_no[j,tndx] += n[j,i]
# calculate average
min_color = -3
max_color = 2
for i in range(ntemps):
  for j in range(ncols):
    zz[j,i] = nan
    if tot_no[j,i]>0:
      if arith_flag == 1:
        print ('j = %i, i = %i, tot_disp = %g, tot_no = %i' % (j,i,tot_disp[j,i],tot_no[j,i]))
        zz[j,i]=log(tot_disp[j,i]/tot_no[j,i])
        if zz[j,i]<min_color: zz[j,i]=min_color
      else:
        zz[j,i]=tot_disp[j,i]/tot_no[j,i]
        if zz[j,i]<min_color: zz[j,i]=min_color

fig = plt.figure(figsize=(13,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
# comment this out to get temp span matched to data
# plt.xlim(min_temp,max_temp)

# define plot
levels = list(range(min_color,max_color,1))
contour = plt.contour(xx_1, yy_1, zz, levels, colors='darkgrey')
plt.clabel(contour,colors = 'k',fmt = '%2.1f', fontsize=12)
cs = plt.contourf(xx_1,yy_1, zz, alpha=1.0, cmap='inferno', vmax=max_color, vmin=min_color, levels=25)
plt.colorbar(cs)

# plot label
if arith_flag==1:
  plt.title('Thin film log(average(mobility)) plot during 1st heating')
else:
  plt.title('Thin film average(log(mobility)) plot during 1st heating')
file_words = indump.split('/')
name_length = len(file_words)
print_name = file_words[name_length-1]
print (print_name)
text_kwargs = dict(ha='left', va='center', fontsize=8, color='black', font='Helvetica', backgroundcolor='white')
x_loc = min_temp + (max_temp-min_temp)/40
plt.text(x_loc,-0.1,print_name,**text_kwargs)
plt.xlabel('Temperature')
plt.ylabel('Fractional height in simulation box')
plt.savefig('mobility_vs_temp_heat.png', dpi=500)
plt.show()
plt.close()
print ("Starting mobility vs temp plot - 2")
# create plot matrix
# data from within max_f_indx, min_f_indx
# in the same temperature range as the first chart
# temp range of x-axis can be separately defined
if min_f_indx > f_tot: min_f_indx = f_tot-1
zz = zeros((ncols,ntemps), dtype=float)
xx = arange(min_temp, max_temp, bin_temp)
yy = arange(0.0, 1.0, delta)
xx_1, yy_1 = meshgrid(xx,yy)
# bin time_steps by temp
# sum mobilities and atom nos
tot_disp = zeros((ncols,ntemps), dtype=float)
tot_no = zeros((ncols,ntemps), dtype=int)
print ('i = %i, ')
for i in range(min_f_indx, min(f_tot, max_f_indx), -1):
  tndx = int(ntemps*(flm_temp[i]-min_temp)/(max_temp-min_temp)+0.5)
  bin_str = ( 'indx = %i tndx = %g, film_temp = %g, max_temp = %g, min_temp = %g' % (i, tndx, flm_temp[i], max_temp, min_temp))
  print (bin_str)
  if tndx>=ntemps: continue
  for j in range(ncols):
      tot_disp[j,tndx] += z[j,i]
      tot_no[j,tndx] += n[j,i]
# calculate average
for i in range(ntemps):
  for j in range(ncols):
    zz[j,i] = nan
    if tot_no[j,i]>0:
      if arith_flag == 1:
        print ('j = %i, i = %i, tot_disp = %g, tot_no = %i' % (j,i,tot_disp[j,i],tot_no[j,i]))
        zz[j,i]=log(tot_disp[j,i]/tot_no[j,i])
        if zz[j,i]<-3: zz[j,i]=-3
      else:
        zz[j,i]=tot_disp[j,i]/tot_no[j,i]
        if zz[j,i]<-3: zz[j,i]=-3

fig = plt.figure(figsize=(13,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
# comment this out to get temp span matched to data
plt.xlim(min_temp,max_temp)

# define plot
levels = list(range(-3,2,1))
contour = plt.contour(xx_1, yy_1, zz, levels, colors='darkgrey')
plt.clabel(contour,colors = 'k',fmt = '%2.1f', fontsize=12)
cs = plt.contourf(xx_1,yy_1, zz, alpha=1.0, cmap='inferno', vmax=max_color, vmin=min_color, levels=25)
plt.colorbar(cs)

# plot label
if arith_flag==1:
  plt.title('Thin film log(average(mobility)) plot during 2nd heat')
else:
  plt.title('Thin film average(log(mobility)) plot during 2nd heat')
file_words = indump.split('/')
name_length = len(file_words)
print_name = file_words[name_length-1]
print (print_name)
text_kwargs = dict(ha='left', va='center', fontsize=8, color='black', font='Helvetica', backgroundcolor='white')
x_loc = min_temp + (max_temp-min_temp)/40
plt.text(x_loc,-0.1,print_name,**text_kwargs)
plt.xlabel('Temperature')
plt.ylabel('Fractional height in simulation box')
plt.savefig('mobility_vs_temp-2heat.png', dpi=500)
plt.show()
plt.close()
