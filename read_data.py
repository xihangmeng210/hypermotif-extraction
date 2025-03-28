import numpy as np
import pickle
import os
import scipy.io

def relabel_edges(uni_nodes,edges):      #relabel nodes and hyperedges
    node_id={}
    node_list=sorted(list(uni_nodes))
    id=0
    for v in node_list:
        id+=1
        node_id[v]=id
    cooked_edges=[[node_id[v] for v in e] for e in edges]
    return cooked_edges

def read_data(name):
    print('#reading hypergraph')
    uni_edges = set()
    uni_nodes = set()
    max_of_e = 0
    uni_edges_list = []

    if name in ['iAF1260b','iJO1366']:
        mat_data = scipy.io.loadmat('data/Metabolism/' + name + '.mat')
        CSC = mat_data[name]['S']
        incidence_matrix = np.array(CSC[0, 0], dtype='int')
        incidence_matrix_T = incidence_matrix.T

        for e in incidence_matrix_T:
            this_e=[]
            for index, x in enumerate(e):
                if x != 0:
                    this_e.append(index)

            this_e_size = len(this_e)
            if (this_e_size == 1) or (this_e_size > 25):
                continue

            max_of_e = max(max_of_e, this_e_size)
            for v in this_e:
                uni_nodes.add(v)

            this_e_tuple = tuple(this_e)
            if this_e_tuple not in uni_edges:
                uni_edges.add(this_e_tuple)
                uni_edges_list.append(this_e)

    elif name in ['ibm03','ibm07','ibm09','ibm10','ibm11','ibm12',
                  'ibm13','ibm14','ibm15','ibm16','ibm17','ibm18']:

        with open('data/VLSI/ISPD98_'+name+'.hgr', 'r') as file:
            lines = file.readlines()
        raw_edges = []
        for line in lines:
            raw_edges.append([int(v) for v in line.strip().split(' ')])

        for this_e in raw_edges:

            this_e_size = len(this_e)
            if (this_e_size == 1) or (this_e_size > 25):
                continue

            max_of_e = max(max_of_e, this_e_size)
            for v in this_e:
                uni_nodes.add(v)

            this_e.sort()
            this_e_tuple = tuple(this_e)
            if this_e_tuple not in uni_edges:
                uni_edges.add(this_e_tuple)
                uni_edges_list.append(this_e)

    N=len(uni_nodes)
    cooked_edges=relabel_edges(uni_nodes,uni_edges_list)
    print("dataname:",name,", vertex number:",N,", hyperedge number:",len(cooked_edges))
    return N,cooked_edges


