#!/local/usr/bin/python3

'''
Used for Bootstrapping analysis with divided blocks on WT-MetaD-HREX simulation results.
1. Divide the COLVAR_ave_reweight into n blocks
2. Calculate the free energies of different substates in each block
3. Random picking up x blocks for t times
4. Calculate the standard deviations of the free energies in the four substates

Example of the COLVAR_ave_reweightx file:
#! FIELDS time s12 s34 h12 h34 metad.bias
 0.000000 18.282737 18.926259 87.418201 92.196632 221.898893
 1.000000 18.793214 18.730255 91.201423 90.731802 221.138819
 2.000000 19.298136 18.103352 94.998325 86.102117 220.458616

Here we take h12 and h34 (or h13 and h24) as the CV space
Keep the header and this script will read the header and extract the column based on it.
#
Usage: python3 ./blocks_bs.py <block num> <x pick> <CV1> <CV2> <bias file>
#
Note: repeating time is set 200 (iteration)
Note: Specifically used for HJ 2D MetaD results. Should be modified for general use.
Interesting observables: free energies of opened, closed, and two intermediate states.
Upper limit of opened state: 5; lower limit of closed state: 50.
Date: 17-02-2022, Author: Zhengyue Zhang
'''

import os
import sys
import numpy as np
import math

def is_number(string):
    # Check if string is a valid number (float)
    try:
        float(string)
        return True
    except ValueError:
        return False


def rand_pick(timescale,bnum,pick):
    # Generate trajectory index list for one iteration
    traj = []
    bsize = round(timescale/bnum)

    blockpick = np.random.choice(bnum,pick).tolist()
    for sample in blockpick:
        if sample != bnum - 1:
            traj = traj + [i for i in range(sample*bsize,(sample+1)*bsize)]
        else:
            traj = traj + [i for i in range(sample*bsize,timescale)]

    return traj


def cal_fe(colvar,traj):
    # Calculate free energies in different substates
    # In each iteration where there is different traj combination
    boundary_close = 50
    boundary_open  = 5
    kbt = 2.478957

    tot_w = 0
    open_w = 0
    close_w = 0
    interi_w = 0
    interii_w = 0

    for frame in traj:
        tot_w += math.exp(colvar[frame,2]/kbt)
        if colvar[frame,0] > boundary_close and colvar[frame,1] > boundary_close:
            open_w += math.exp(colvar[frame,2]/kbt)
        elif colvar[frame,0] < boundary_open and colvar[frame,1] < boundary_open:
            close_w += math.exp(colvar[frame,2]/kbt)
        elif colvar[frame,0] > boundary_close and colvar[frame,1] < boundary_open:
            interi_w += math.exp(colvar[frame,2]/kbt)
        elif colvar[frame,0] < boundary_open and colvar[frame,1] > boundary_close:
            interii_w += math.exp(colvar[frame,2]/kbt)

    opened = -kbt * np.log(open_w/tot_w)
    closed = -kbt * np.log(close_w/tot_w)
    interi = -kbt * np.log(interi_w/tot_w)
    interii = -kbt * np.log(interii_w/tot_w)

    return closed - opened, interi - opened, interii - opened


if __name__ == "__main__":
    
    # First import arguments and check
    if len(sys.argv) != 6:
        print("Usage: python3 ./blocks_bs.py <block num> <x pick> <CV1> <CV2> <bias file>")
        quit()
    else:
        if not is_number(sys.argv[1]) or not is_number(sys.argv[2]):
            print("Block number should be number")
            quit()
        if not os.path.isfile(sys.argv[5]):
            print("Bias file not exists!")
            quit()
        bnum  = int(sys.argv[1])
        pick  = int(sys.argv[2])
        cv1   = sys.argv[3]
        cv2   = sys.argv[4]
        fname = sys.argv[5]

    # Get values from bias file
    colvar = np.loadtxt(fname, comments = '#')
    # Get the header
    with open(fname, 'r') as inp:
        header = inp.readline()
        if header.startswith('#'):
            header = header.replace("#! FIELDS ","").split()
            if not cv1 in header or not cv2 in header:
                print("Input CVs not in the file!")
                quit()
            ind_cv1 = header.index(cv1)
            ind_cv2 = header.index(cv2)
            ind_bias = header.index("metad.bias")
    
    # Only keep our interesting CVs and the bias
    colvar = colvar[:,[ind_cv1,ind_cv2,ind_bias]]

    # Generate the index of each block, then do 200 random pick on the blocks
    # The picked blocks (x blocks) are represented by the index list
    # Extracting CVs and bias from colvar based on index list, then calculate free energies
    iteration = 200
    timescale = colvar.shape[0]

    all_opened = []
    all_closed = []
    all_interi = []
    all_interii = []
    for i in range(0,iteration):
        traj = rand_pick(timescale,bnum,pick)
        closed, interi, interii = cal_fe(colvar,traj)
        all_closed.append(closed)
        all_interi.append(interi)
        all_interii.append(interii)
    
    with open('boot.dat','w') as op:
        op.write("# Closed interi interii\n")
        op.write(str(np.std(all_closed)) + "\t" + str(np.std(all_interi)) + "\t" + str(np.std(all_interii)) + "\t" + fname + "\n")














