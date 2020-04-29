#script to run tumor spatial growth simulation

import itertools
import operator
import copy
import numpy as np
import pandas as pd
import math
import operator
from matplotlib import pyplot as plt
import matplotlib as mpl
import baltic as bt
from statsmodels.nonparametric.api import KernelReg
import pickle
import random
import matplotlib.cm as cm
import os
import random

os.chdir("/Users/mayalewinsohn/Documents/PhD/Bedford_lab/hcc-tumor-evol")

def adjacent_cell(lattice,i,j): #function to check if adjacent cells are empty
    free_cells = []
    change = [1,-1]
    for c in change:
        checki = lattice[i+c,j]
        if checki == None:
            free_cells.append([i+c,j])
        checkj = lattice[i, j+c]
        if checkj == None:
            free_cells.append([i,j+c])
        for h in change:
            checkd = lattice[i+h,j+c]
            if checkd == None:
                free_cells.append([i+h,j+c])
    return free_cells

def cancer_sim(founder_cells = 1,deltat=(1/float(24)),CSC = True,max_cells = 1000,proliferationCSC = 1, proliferationCC = 2, motility = 15, pmax = 10, mortality = 0.1, mutation = 0.1, selection = False,ps = 0.05):    

    #create classes
    class Cell:
        def __init__(self):
            
            self.parent_index = 0 #index of cell it was derived from (for df conversion)
            self.parent = None #cell it was derived from
            self.children = [] #set of daughter cells
            self.locx = None #x location in lattice
            self.locy = None #y location in lattice
            self.birthdate = None
            self.deathdate = None
            self.pmax = pmax #proliferation potential
            self.cellnum = None #unique cell number to identify
            self.alpha = mortality #probability of sponaneous death
            self.mutation_rate = mutation #probability of mutation per cell division
            self.motility = motility #probability of migration per day
            self.proliferation_rate = proliferationCC #average cell divisions per time unit
            self.clone = None
            self.mutations = []
            self.ps = 0
            self.curnode = None
            self.new_mut = False
            self.group = None
            self.alive = True
            
        def as_dict(self):
            return {'index': self.cellnum,
                    'parent_index': self.parent_index,
                    'birthdate': self.birthdate,
                    'locx': self.locx,
                    'locy': self.locy,
                    'deathdate': self.deathdate,
                    'mutations': self.mutations,
                    'proliferation_rate': self.proliferation_rate,
                    'alpha': self.alpha,
                    'mutation_rate': self.mutation_rate}


    class StemCell(Cell): #define stem cell class
        def __init__(self):
            Cell.__init__(self)
            self.pmax = float('inf') #infinite proliferation potential
            self.alpha = 0 #immortal
            self.ps = ps #probability of symmetric division
            self.proliferation_rate = proliferationCSC
            self.alive = True
    
    class Group():
        def __init__(self):
            self.alive_cells = []
            self.node = None
            self.id = None


    #create lattice array
    N = max_cells # starting demensions of lattice
    lattice = np.empty( (N,N), dtype=object)

    #parameters
    proliferation_rate_CC = proliferationCC #proliferation rate of clonal cell
    proliferation_rate_CSC = proliferationCSC #proliferation rate of cancer stem cell
    pmax_CC = pmax
    migration_potential = motility
    alpha_CC = mortality
    cur_cellnum = 0
    cur_mutnum = 0
    cur_clonenum = 1
    cur_index = 1
    mutation_rate = mutation
    time = 0

    alive_cells = []
    cells = []
    mutations = []
#     clones = ['a','b','c','d']

    i = int(N/2) #to place founder cell in center of lattice
    modi = [0,0,1,1]
    modj = [0,1,0,1]
    if CSC == True:
        for x in range(founder_cells):
            Fcell = StemCell() #initiate founder cells, stem-cells
            Fcell.ps = ps
            Fcell.proliferation_rate = proliferation_rate_CSC #average cell divisions per day
            Fcell.locx = i+modi[x] #x location in lattice
            Fcell.locy = i+modj[x] #y location in lattice
            Fcell.birthdate = 0
            Fcell.cellnum = cur_cellnum #unique cell number to identify
            cur_cellnum += 1
            Fcell.mutation_rate = mutation_rate #probability of mutation per cell division
            Fcell.motility = migration_potential #probability of migration per day

            cells.append(Fcell)
            alive_cells.append(Fcell)

            lattice[Fcell.locx,Fcell.locy] = Fcell
    else:
        for x in range(founder_cells):
            Fcell = StemCell() #initiate founder cells, non-stem cells
            Fcell.proliferation_rate = proliferation_rate_CC
            Fcell.locx = i+modi[x] #x location in lattice
            Fcell.locy = i+modj[x] #y location in lattice
            Fcell.alpha = mortality
            Fcell.ps = 0
            Fcell.birthdate = 0
            Fcell.cellnum = cur_cellnum  #unique cell number to identify
            cur_cellnum += 1
            Fcell.mutation_rate = mutation_rate #probability of mutation per cell division
            Fcell.motility = migration_potential #probability of migration per day

            cells.append(Fcell)
            alive_cells.append(Fcell)
            lattice[Fcell.locx,Fcell.locy] = Fcell


    #time parameters
    dt = deltat # time is equilavent to 1/24 of a day or 1 hour
#     stop_time =stime#stop simulation after this many days

    
    #while time < stop_time:
    while len(alive_cells) < max_cells:
        if len(alive_cells) < 1:
            print('no cells alive')
            break
        time += dt
        cell_stack = random.sample(alive_cells,len(alive_cells)) #random order of cells
        for cell in cell_stack:
            alive = True
            r = random.uniform(0, 1)
            pd = cell.proliferation_rate * dt #probability of proliferation in time dt
            free_cells = adjacent_cell(lattice,cell.locx,cell.locy)
            if r < pd: # Does cell attempt to divide? 
                r = random.uniform(0,1)
                pdie = cell.alpha * pd #probability of spontaneously dying
                if r < pdie:
                    alive = False
                elif len(free_cells) > 0: #is there any space to divide?
                    if cell.pmax > 0: #is cell proliferation capacity exhausted?
                        cell.pmax -= 1
                        r = random.uniform(0,1)
                        new_cell = copy.deepcopy(cell)
                        new_cell_2 = copy.deepcopy(cell)
                        
                        if r <= cell.ps: #does cell divide asymmetrically?
                            new_cell.pmax = pmax
                            new_cell.proliferation_rate = proliferation_rate_CC #average cell divisions per day
                            new_cell.alpha = mortality
                            new_cell.ps = 0

                        new_cell.parent = cell
                        new_cell_2.parent = cell
                        newloc = random.choice(free_cells)
                        new_cell.locx = newloc[0] #x location in lattice
                        new_cell.locy = newloc[1] #y location in lattice
                        new_cell_2.locx = cell.locx #x location in lattice
                        new_cell_2.locy = cell.locy #y location in lattice
                        new_cell.birthdate = time
                        new_cell_2.birthdate = time
                        new_cell.cellnum = cur_cellnum #unique cell number to identify
                        cur_cellnum += 1
                        new_cell_2.cellnum = cur_cellnum
                        new_cell.parent_index = new_cell.parent.cellnum
                        new_cell_2.parent_index = new_cell_2.parent.cellnum
                        cur_cellnum += 1
                        new_cell.children = []
                        new_cell_2.children = []


                        r1 = random.uniform(0,1)
                        if r1 < cell.mutation_rate: #does cell gain a mutation?
                            #new_mut = Mutation()
                            #new_mut.mutnum = cur_mutnum
                            new_mut = cur_mutnum
                            cur_mutnum +=1
                            new_cell.mutations.append(new_mut)
                            cell.mutations.append(new_mut)
                            new_cell.new_mut = True
                            mutations.append(new_mut)
                        
                        r2 = random.uniform(0,1)
                        if r2 < cell.mutation_rate: #does cell gain a mutation?
                            #new_mut = Mutation()
                            #new_mut.mutnum = cur_mutnum
                            new_mut = cur_mutnum
                            cur_mutnum +=1
                            new_cell_2.mutations.append(new_mut)
                            cell.mutations.append(new_mut)
                            new_cell_2.new_mut = True
                            mutations.append(new_mut)



                        cell.children.append(new_cell)
                        cell.children.append(new_cell_2)
                        alive_cells.append(new_cell)
                        alive_cells.append(new_cell_2)
                        lattice[new_cell.locx,new_cell.locy] = new_cell
                        lattice[cell.locx,cell.locy] = new_cell_2
                        cells.append(new_cell)
                        cells.append(new_cell_2)
                        alive_cells.remove(cell)
                        cell.alive = False
                        cell.deathdate = time
                    else:
                        lattice[cell.locx,cell.locy] = None
                        alive_cells.remove(cell)
                        cell.alive = False
                        cell.deathdate = time



    for cell in alive_cells:

        cell.deathdate = time


    return cells, alive_cells, lattice
#run simulation

random.seed(1029)
cells_CC, alive_cells_CC, lattice_CC = cancer_sim(founder_cells = 1,max_cells = 1000,proliferationCSC = 1,proliferationCC = 1,CSC = False,pmax = 10,mutation = 0.1,mortality = 0.5)

random.seed(3872)
cells_CSC, alive_cells_CSC, lattice_CSC = cancer_sim(founder_cells = 1,max_cells = 1000,proliferationCSC = 1,proliferationCC = 1,CSC = True,pmax = 10,mutation = 0.1,ps = 0.1,mortality = 0.5)

#convert results to dataframe

CC_df = pd.DataFrame([x.as_dict() for x in cells_CC])
CSC_df = pd.DataFrame([x.as_dict() for x in cells_CSC])

#save simulation output to csv
CC_df.to_csv('outputs/sim_CC_cells.csv')
CSC_df.to_csv('outputs/sim_CSC_cells.csv')
