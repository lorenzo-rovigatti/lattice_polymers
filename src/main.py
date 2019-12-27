#!/usr/bin/env python3

import sys
import math
import numpy as np

L_COORD = 6
L_DIRS = [
    np.array(( 1,  0,  0)),
    np.array(( 0,  1,  0)),
    np.array(( 0,  0,  1)),
    np.array((-1,  0,  0)),
    np.array(( 0, -1,  0)),
    np.array(( 0,  0, -1))
    ]
L_DIM = 100
LATTICE = np.zeros((L_DIM, L_DIM, L_DIM))
TRIES = 1000

def get_lattice_value(coordinate):
    c = coordinate % L_DIM
    return LATTICE[tuple(c)]

def set_lattice_value(coordinate, value):
    c = coordinate % L_DIM
    LATTICE[tuple(c)] = value

class Chain():
    def __init__(self):
        self.monomers = []
        self.r_ee_sqr = None
        
    def add_monomer(self, m):
        self.monomers.append(m)
        set_lattice_value(m, 1)
        if len(self.monomers) > 1:
            dist = self.monomers[-1] - self.monomers[0]
            self.r_ee_sqr = np.dot(dist, dist)
            
    def new_monomer_energy(self, m):
        energy = 0.
        if get_lattice_value(m):
            return np.inf
        return energy
            
    def size(self):
        return len(self.monomers)
    
    
def delete_chain(chain):
    for m in chain.monomers:
        set_lattice_value(m, 0)
        
        
def build_chain_rosenbluth(N):
    new_chain = Chain()
    last_inserted = np.random.randint(0, L_DIM, 3)
    new_chain.add_monomer(last_inserted)
    new_chain.rosenbluth_weight = float(L_COORD)
    
    while new_chain.size() < N:
        weights = []
        for dir in L_DIRS:
            attempt = last_inserted + dir
            energy = new_chain.new_monomer_energy(attempt)
            if energy == np.inf:
                weights.append(0.)
            else:
                weights.append(math.exp(-energy))
        tot_weight = sum(weights)
        if tot_weight == 0.:
            # dead end, we start over
            delete_chain(new_chain)
            new_chain = Chain()
            last_inserted = np.random.randint(0, L_DIM, 3)
            new_chain.add_monomer(last_inserted)
            new_chain.rosenbluth_weight = float(L_COORD)
        else:
            # choose a direction
            rnd = np.random.random()
            tot = 0.
            for i, w in enumerate(weights):
                tot += w / tot_weight
                if tot > rnd:
                    break
                
            last_inserted = last_inserted + L_DIRS[i]
            new_chain.add_monomer(last_inserted)
            new_chain.rosenbluth_weight *= tot_weight
            
    return new_chain
        

for N in range(40, 100, 5):
    r_ee_sqr = 0.
    tot_weight = 0.
    t = 0
    last_chain = build_chain_rosenbluth(N)
    while t < TRIES:
        chain = build_chain_rosenbluth(N)
        if np.random.random() < chain.rosenbluth_weight / last_chain.rosenbluth_weight:
            last_chain = chain
        t += 1
        r_ee_sqr += chain.r_ee_sqr
        delete_chain(chain)
    r_ee_sqr /= TRIES - 1
    print(N - 1, file=sys.stderr)
    print(N - 1, r_ee_sqr)
