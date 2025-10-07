#!/usr.bin/env python
# Inputs: atomfile describing the set of atoms in this simulation
#         first dump file created during anneal of deposited film recording time series of atom positions
#         program will find and open the rest of atom position files
#
# Calculate mean square xy offset of atoms in middle of film vs log time (z offset is affected by expansion)
# calculate distribution of offsets for each atom dump relative to first observation
# Collate with delay time
# Plot distribution vs temperature
# ------ potentially sort by atom type (or at least large/med/small)
#
# Output: mobility.txt: scratch pad for organizing data from dump-file
# Output: mob_dist_vs_time.txt: line plots of displacement spectrum for delay times
#           one line for each delay time
# Output: out_tmd.txt. displacement spectrum for each temperature bin
#
# R.B. Stephens Aug 2025
#
# Written for Python3
# command: Python3 intermedscatt.py infile#1 infile#2 infile#3
#
import os, sys, math
from numpy import *

# the dump file input is any in the series of dump files
try:
  inatom = sys.argv[1]; indump = sys.argv[2];
except:
  print("Usage:", sys.argv[0], "(1) inatom (2) indump")
  sys.exit(1)

print("start")

if os.path.exists(indump):
  mydirect = os.path.dirname(indump)
  myfilename = os.path.basename(indump)
else:
  print("Usage:", sys.argv[0], "(1) inatom (2) indump")
  sys.exit(1)
mystring = myfilename.split('.')
if len(mystring)>0:
  basename = mystring[0]
# now clip the last 9 digits to get name stem
bas_len = len(basename) - 9
print("File name:", basename, "length = ",len(basename))
print("File basename:", basename[0:(len(basename)-9)])
if bas_len>0:
  basestem = basename[0:bas_len]
  print("File basename:", basename)
else:
  print("File name:", basename, " is too short")
  sys.exit(1)
# then find the files by adding suffix to the stem - continue until end of file_indx or no more files
file_indx = ['000000000', '000000010', '000000012', '000000015', '000000017', '000000020', '000000027', '000000035', '000000042', '000000050', '000000062', '000000075', '000000087', '000000100', '000000120', '000000150', '000000170', '000000200', '000000270', '000000350', '000000420', '000000500', '000000620', '000000750', '000000870', '000001000', '000001200', '000001500', '000001700', '000002000', '000002700', '000003500', '000004200', '000005000', '000006200', '000007500', '000008700', '000010000', '000012500', '000015000', '000017500', '000020000', '000027500', '000035000', '000042500', '000050000', '000062500', '000075000', '000087500', '000100000', '000125000', '000150000', '000175000', '000200000', '000275000', '000350000', '000425000', '000500000', '000625000', '000750000', '000875000', '001000000', '001250000', '001500000', '001750000', '002000000', '002750000', '003500000', '004250000', '005000000', '006250000', '007500000', '008750000', '010000000', '012500000', '015000000', '017500000', '020000000', '027500000', '035000000', '042500000', '050000000', '062500000', '075000000', '087500000', '100000000', '125000000', '150000000', '175000000', '200000000']
ntime = len(file_indx)
file_end = ".film"

# displacment range for calculation (geometrically spaced intervals)
max_log_disp = 1
min_log_disp = -1.66667
bin_log_disp = 0.1
ndisp   = int((max_log_disp-min_log_disp)/bin_log_disp+0.5)

# needed if atoms have been lost - then max index no is larger than total atom no.
evap_no = 200
# calculate average mobility and density between height limits (assume film goes to 36)
mob_top = 30.0
mob_bot = 10.0

# open atom file for reading simulation parameters
# =======================================
iatom = open( inatom, 'r')

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

print('finished reading atom file, type no = %i' % (at_type_no))

# -- record calcularted density at each time
dens = zeros((ntime+1), dtype=float)
# from first dump file record location, type and id of each atom in array
# -- that data will be stored in fpos[0] and will be the reference for all subsequent dump files - each stored in fpos[1]

disp_time = zeros((ntime+1, ndisp+1), dtype=int)
typ_time  = zeros((ntime+1, ndisp+1), dtype=int)
tot_time  = zeros((ntime+1), dtype=int)
# axis values
time_step  = zeros((ntime+1), dtype=float)
displac    = zeros((ndisp+1), dtype=float)
log_disp   = zeros((ndisp+1), dtype=float)

for dndx in range(ndisp):
  log_disp[dndx]=dndx*bin_log_disp+min_log_disp
  displac[dndx]=10**(dndx*bin_log_disp+min_log_disp)
for tndx in range(ntime):
  time_step[tndx]=float(file_indx[tndx])
time_min = time_step[0]
time_max = time_step[ntime]

# record displacement distribution
# - number vs log10 displacement
# - sum typeID vs log10 displacement
# - and density at each time

f_no = -1
# time delay index
f_indx = 0

# create temp-mob-dens file
out_file_name = mydirect + "/" + "time-dens-disp.txt"
out_tmd = open(out_file_name, 'w')
# =================================
# create labels for mob-dens file
out_tmd.write('Density Displacement spectrum and average type at each time step from %5.2f to %5.2f from files %s <nnnnnnn> in directory %s \n' % (time_min,time_max, basestem, mydirect))
out_tmd.write('For log10 displacement values between %5.2f and %5.2f \n' % (min_log_disp,max_log_disp))
out_tmd.write('And for atoms in the height range between %5.2f and %5.2f \n' %(mob_bot,mob_top))
out_tmd.write('The displacement spectrum has %i time steps geometrically spaced \n ' % (ndisp))


out_tmd.write('time, density, disp tot, log10(displacement), \n')
out_tmd.write('-- , -- , ')
for indx in range(ndisp):
  out_tmd.write(', %6.4f' % (displac[indx]))
out_tmd.write(',   ,   ,  ')
for indx in range(ndisp):
  out_tmd.write(', %6.4f' % (displac[indx]))
out_tmd.write(' \n')
out_tmd.write('-- , -- , ')
for indx in range(ndisp):
  out_tmd.write(', %5.2f' % (log_disp[indx]))
out_tmd.write(',   ,   ,  ')
for indx in range(ndisp):
  out_tmd.write(', %5.2f' % (log_disp[indx]))
out_tmd.write(' \n')
out_tmd.close


# find atom file to open - the suffixes in order of file_indx
for fndx in range(ntime):
  indump = mydirect + "/" + basestem + file_indx[fndx] + file_end
  print (indump)
  if not(os.path.exists(indump)):
    print("reached end of time sequence")
    break
  print ("opening file no %i: %s" % (fndx, indump))
  idump = open(indump, 'r')
  # read all the lines - break at end of file
  while 1:
    line = idump.readline()
    if not line: break # jump out of loop at the end of the file
    in_words = line.split()
    if in_words[0] =="ITEM:":
      item_id = in_words[1];
      # "TIMESTEP" is the keyword in the first line of each frame
      if item_id == "TIMESTEP":
        f_indx += 1
        f_no += 1
        if f_no>1: f_no = 1
        line = idump.readline()
        if not line: break
        in_words = line.split()
        step_no = int(in_words[0])
      # "NUMBER" = the number of atoms; therefore the number of lines in the ATOMS block
      elif item_id == "NUMBER":
        line = idump.readline()
        if not line: break
        in_words = line.split()
        if f_no==0:
          atom_no = int(in_words[0])
          padded_atm_no = atom_no+evap_no
          pos = zeros(padded_atm_no)
          # atom position and type array for reference file [0] and current delay file [1]
          fpos  = zeros([2,padded_atm_no,4])
          # initialize output arrays
        tot_num1 = zeros(ndisp+1)
        reg_vol = 0 # zero volume of atoms inside density measuring region
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
        width = htx
        ch_width = width/2.0
      elif item_id == "ATOMS":
        # the rest of this keyword line gives variable names - we identify columns containing atom_id and x,y,z coords
        for i in range(0, len(in_words), 1):
          if (in_words[i] == 'id'): # note some ids may be missing - atoms evaporated
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
        # put list of id, type, and position of each atom into list for current time_step
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
          # changes value to simulation units from fraction of width
          x_cpos = float(in_words[i_xs])*htx
          y_cpos = float(in_words[i_ys])*hty
          z_cpos = float(in_words[i_zs])*htz
          at_pos = [x_cpos,y_cpos,z_cpos,at_typ]
          fpos[f_no][at_id] = at_pos
          # do volume sum for density
          if (z_cpos >= mob_bot)& (z_cpos <= mob_top):
            reg_vol += atm_vol[at_typ]
            dens_no += 1
          # calculate displacements after the first file
          if fndx>0:
            # calculate atom displacement (in the data file, the coordinates xs, ys, zs are scaled to simulation box size
            this_frame = fpos[1]
            prev_frame = fpos[0]
            dispx = zeros((padded_atm_no))
            dispy = zeros((padded_atm_no))
            dispz = zeros((padded_atm_no))
            disp  = zeros((padded_atm_no,2))
            offsetx = 0.0
            offsety = 0.0
            for iatm in range(padded_atm_no):
              dispx[iatm] = this_frame[iatm][0] - prev_frame[iatm][0]
              dispy[iatm] = this_frame[iatm][1] - prev_frame[iatm][1]
              dispz[iatm] = this_frame[iatm][2] - prev_frame[iatm][2]
              # correct displacements for periodic boundary condition offsets
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
            for iatm in range(padded_atm_no):
              if iatm==0:continue
              disp[iatm,0] = sqrt((dispx[iatm]-offsetx)**2 + (dispy[iatm]-offsety)**2)
              disp[iatm,1] = (this_frame[iatm][2] + prev_frame[iatm][2])/2.0
            # sum occurrence binned by displacement
            for iatm in range(padded_atm_no):
              if (this_frame[iatm][2]>mob_bot) and (this_frame[iatm][2]<mob_top): # atom is within measuring height
                disp_bin = int((log10(abs(disp[iatm,0]))-min_log_disp)/(max_log_disp-min_log_disp)*ndisp+0.5) # calculate displacement index
                if (disp_bin<0): disp_bin = 0           # - collapse extremes into selected range
                if (disp_bin>ndisp): disp_bin = ndisp
                disp_time[fndx][disp_bin] += 1
                tot_time[fndx] += 1
                typ_time[fndx][disp_bin]  += this_frame[iatm][3]
        # use htx, hty because box size can change from frame to frame
        dens[fndx] = reg_vol/(htx*hty*(mob_top-mob_bot))
        # write the displacement spectra and atom type for this atom file
        if fndx>0:
          print ('Writing data for step = %i, time = %i' % (fndx, time_step[fndx]))
          out_tmd = open(out_file_name, 'a')
          out_tmd.write(' %5.4f ' % (time_step[fndx]))
          out_tmd.write(', %6.5f ' % (dens[fndx]))
          out_tmd.write(', %i' %(tot_time[fndx]))
          for dndx in range(ndisp):
            out_tmd.write(', %5i' % (disp_time[fndx][dndx]))
          out_tmd.write(',   ,   ,  ')
          for dndx in range(ndisp):
            out_tmd.write(', %5i' % (typ_time[fndx][dndx]))
          out_tmd.write(' \n')
          out_tmd.close
      else: break
print ('Finished %i times' % (f_indx))
