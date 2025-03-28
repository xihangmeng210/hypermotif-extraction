import random
from itertools import permutations
import numpy as np
import pickle
from multiprocessing import Pool
from  tqdm import  tqdm
import time
import os
import psutil

if not os.path.exists('subprocess_result'):
    os.mkdir('subprocess_result')

class hypergraph:

    def __init__(self, N, edges,order):
        self.N = N  # number of vertices
        self.edges = edges     # edge set
        self.E = len(self.edges)    # number of edges

        self.Ego,self.S = self.get_Ego() # Ego: edges containing each v in H, S: edge sizes
        self.Snb,self.Lnb,self.Bro = self.get_neibors()  #neibor information of edges

        self.order = order   # order of motifs
        self.W = [2**(self.order-1-i) for i in range(self.order)]    # initial edge weights
        self.perms = list(permutations(range(self.order)))     # permutations
        self.NewIndexMap = get_new_index_maps(self.order,self.W,self.perms)

    def get_Ego(self):
        Ego = [set() for _ in range(self.N+1)]
        S = []
        for i, e in enumerate(self.edges):
            for v in e:
                Ego[v].add(i)
            S.append(len(e))
        return Ego, S

    def get_neibors(self):
        Snb,Lnb,Bro=[],[],[]  #Snb: neibor sets of e; Lnb: neibor lists of e; Bro: brothers of e
        for i,e in enumerate(self.edges):
            this_neibor_set = set()
            this_brother_list = []
            for v in e:
                for ego_e in self.Ego[v]:
                    if ego_e == i:
                        continue
                    if ego_e not in this_neibor_set:
                        this_neibor_set.add(ego_e)
                        if ego_e > i:
                            this_brother_list.append(ego_e)
            this_neibor_list = sorted(this_neibor_set, reverse=True)  # descending order
            Snb.append(this_neibor_set)
            Lnb.append(this_neibor_list)
            Bro.append(this_brother_list)
        return Snb, Lnb, Bro


    def add(self, e,Esub, t, index_map, w):  #add e to index_map
        Esub.append(e)
        t[0]+=self.S[e]
        for v in self.edges[e]:
            if v in index_map:
                index_map[v] += w
            else:
                index_map[v] = w

    def back(self,e, Esub,t,index_map, w):  # remove e from index_map
        Esub.pop()
        t[0]-=self.S[e]
        if self.S[e]>t[0]:
            index_map.clear()
            t[0]=0
            for i,ei in enumerate(Esub):
                self.add(ei,[],t,index_map, self.W[i])
        else:
            for v in self.edges[e]:
                index_map[v] -= w
                if index_map[v] == 0:
                    del index_map[v]


    def update(self,index_map, C): # calculate encoding vector T based on index_map
        T = [0] * (2 ** self.order - 1)
        for v in index_map:
            T[index_map[v] - 1] += 1 # -1 is due to the removal of the first entry
        S = t_to_s(T)     # convert T to string form to reduce memory usage
        if S in C:           #update its occurrences in Count C
            C[S] += 1
        else:
            C[S] = 1

    def codingtransform(self, T, pi):
        new_T = [0] * (2 ** self.order - 1)
        for i in range(2 ** self.order - 1):
            if T[i] == 0:
                continue
            new_index = self.NewIndexMap[pi][i + 1]
            new_T[new_index - 1] = T[i]
        return tuple(new_T)

    def Isomorphism_Resolver(self,raw_C):        # aggregate the occurrences of each nonisomorphic pattern in C
        cooked_C = {}
        motifs = [S for S in raw_C]
        for S in motifs:
            if S not in raw_C:
                continue

            this_value = raw_C[S]
            del raw_C[S]

            T = s_to_t(S)
            max_T = tuple(T)

            for pi in range(1, len(self.perms)):
                new_T = self.codingtransform(T, pi)
                new_S = t_to_s(new_T)
                if new_S in raw_C:
                    this_value += raw_C[new_S]     # aggregate the occurrences
                    del raw_C[new_S]
                max_T = max(max_T, new_T)

            cooked_C[t_to_s(max_T)] = this_value  #The lexicographically largest one serves as the index for the isomorphism class.
        return cooked_C

    def ternary_hmotif_census(self,taskid,slice):
        P=P_global
        C = {}
        for e0 in slice:
            if random.random() >= P[0]:
                continue
            E_sub,t,index_map =[],[0],{}
            self.add(e0,E_sub, t, index_map, self.W[0])

            for j in range(len(self.Bro[e0])):
                if random.random() >= P[1]:
                    continue
                e1 = self.Bro[e0][j]
                self.add(e1, E_sub,t, index_map, self.W[1])

                for e2 in self.Lnb[e1]:
                    if e2 > e0:
                        if e2 not in self.Snb[e0]:
                            if random.random() >= P[2]:
                                continue
                            self.add(e2, E_sub,t, index_map, self.W[2])
                            self.update(index_map, C)                      #Case 3
                            self.back(e2,E_sub,t, index_map, self.W[2])
                    else:
                        break

                for k in range(j + 1, len(self.Bro[e0])):
                    if random.random() >= P[2]:
                        continue
                    e2 = self.Bro[e0][k]
                    self.add(e2, E_sub, t, index_map, self.W[2])
                    self.update(index_map, C)                    # Case 1&2
                    self.back(e2, E_sub, t, index_map, self.W[2])

                self.back(e1, E_sub, t, index_map, self.W[1])

        C = self.Isomorphism_Resolver(C)
        with open('subprocess_result/' + str(taskid), 'wb') as file:
            pickle.dump(C, file)

    def quaternary_hmotif_census(self, taskid, slice):
        P= P_global
        C = {}
        for e0 in slice:
            if random.random() >= P[0]:
                continue
            Esub, t, index_map = [], [0], {}
            self.add(e0, Esub, t, index_map, self.W[0])

            for j in range(len(self.Bro[e0])):
                if random.random() >= P[1]:
                    continue
                e1 = self.Bro[e0][j]
                self.add(e1, Esub, t, index_map, self.W[1])
                N_e1_older_e0_list = []
                N_e1_older_e0_set = set()
                for e2 in self.Lnb[e1]:
                    if e2 > e0:
                        if e2 not in self.Snb[e0]:
                            N_e1_older_e0_list.append(e2)
                            N_e1_older_e0_set.add(e2)
                            if random.random() >= P[2]:
                                continue
                            self.add(e2, Esub, t, index_map, self.W[2])
                            for e3 in self.Lnb[e2]:
                                if e3 > e0:
                                    if (e3 not in self.Snb[e0]) and (e3 not in self.Snb[e1]):
                                        if random.random() >= P[3]:
                                            continue
                                        self.add(e3, Esub, t, index_map, self.W[3])
                                        self.update(index_map, C)                       #Case 11
                                        self.back(e3, Esub, t, index_map, self.W[3])
                                else:
                                    break
                            for e3 in N_e1_older_e0_list[:-1]:
                                if random.random() >= P[3]:
                                    continue
                                self.add(e3, Esub, t, index_map, self.W[3])
                                self.update(index_map, C)                       #Case 6,7
                                self.back(e3, Esub, t, index_map, self.W[3])

                            self.back(e2, Esub, t, index_map, self.W[2])
                    else:
                        break
                for m in range(j + 1, len(self.Bro[e0])):
                    if random.random() >=P[2]:
                        continue
                    e2 = self.Bro[e0][m]
                    self.add(e2, Esub, t, index_map, self.W[2])
                    for e3 in N_e1_older_e0_list:
                        if random.random() >= P[3]:
                            continue
                        self.add(e3, Esub, t, index_map, self.W[3])
                        self.update(index_map, C)                             #Case 2,4,9,10
                        self.back(e3, Esub, t, index_map, self.W[3])
                    for e3 in self.Lnb[e2]:
                        if e3 > e0:
                            if e3 not in self.Snb[e0]:
                                if e3 not in N_e1_older_e0_set:
                                    if random.random() >= P[3]:
                                        continue
                                    self.add(e3, Esub, t, index_map, self.W[3])
                                    self.update(index_map, C)                              #Case 2,4,9,10
                                    self.back(e3, Esub, t, index_map, self.W[3])
                        else:
                            break

                    for n in range(m + 1, len(self.Bro[e0])):
                        if random.random() >= P[3]:
                            continue
                        e3 = self.Bro[e0][n]
                        self.add(e3, Esub, t, index_map, self.W[3])
                        self.update(index_map, C)                          #Case 1,3,5,8
                        self.back(e3, Esub, t, index_map, self.W[3])

                    self.back(e2, Esub, t, index_map, self.W[2])

                self.back(e1, Esub, t, index_map, self.W[1])

        C = self.Isomorphism_Resolver(C)
        with open('subprocess_result/' + str(taskid), 'wb') as file:
            pickle.dump(C, file)


H_global = None
P_global =None


def census_worker(taskid, slice):
    if H_global.order == 3:
        return H_global.ternary_hmotif_census(taskid, slice)
    elif H_global.order == 4:
        return H_global.quaternary_hmotif_census(taskid, slice)
    else:
        print("No such configuration!")

def multiprocess_run_census(H,P,num_of_processe,num_of_tasks):
    global H_global,P_global
    H_global,P_global=H,P

    tasks_slices = split_into_k_parts(H.E, num_of_tasks)
    print("#taking a census")
    pbar = tqdm(total=num_of_tasks)
    count_start = time.time()
    my_pool=Pool(num_of_processe)
    for taskid in range(num_of_tasks):
        my_pool.apply_async(census_worker, args=(taskid, tasks_slices[taskid]), callback=lambda aaa: pbar.update(1))
    my_pool.close()
    my_pool.join()
    pbar.close()
    count_end = time.time()
    print("#finished, time cost:", count_end - count_start)

    total_C = {}
    total_times = 0
    for taskid in range(num_of_tasks):
        with open('subprocess_result/' + str(taskid), 'rb') as file:
            this_C = pickle.load(file)
        for S in this_C:
            total_times += this_C[S]
            if S in total_C:
                total_C[S] += this_C[S]
            else:
                total_C[S] = this_C[S]
        os.remove('subprocess_result/' + str(taskid))
    print('number of nonisomorphic patterns:', len(total_C))
    print('total counts:', total_times)
    return total_C




num_to_str = {}
for i in range(0, 25):
    if i < 10:
        num_to_str[i] = str(i)
    else:
        num_to_str[i] = chr(ord('A') + i - 10)

str_to_num = {}
for i in range(0, 25):
    if i < 10:
        str_to_num[str(i)] = i
    else:
        str_to_num[chr(ord('A') + i - 10)] = i


def t_to_s(T):    #convert T to string form
    S = ''
    for t in T:
        S += num_to_str[t]
    return S

def s_to_t(S):       #convert S to tuple form
    return tuple(str_to_num[s] for s in S)


def get_new_index_maps(order, W, perms):
    binary_str_map = {}
    for index in range(1, 2 **order):
        binary_str_map[index]= bin(index)[2:].zfill(order)

    initial_weight_np = np.array([W])
    NewIndexMaps = {}
    for i in range(1, len(perms)):
        this_perm = perms[i]
        this_weight = initial_weight_np[:, this_perm][0]  # new weights for this permutation
        this_new_index_map = {}
        for index in range(1, 2 ** order):
            new_index = 0
            for x in range(order):
                new_index += int(binary_str_map[index][x]) * this_weight[x]
            this_new_index_map[index] = new_index
        NewIndexMaps[i] = this_new_index_map
    return NewIndexMaps

def split_into_k_parts(M, k):
    part_size = M // k
    result = []
    start = 0
    for i in range(k):
        end = start + part_size
        if i < M % k:
            end += 1
        result.append(list(range(start, end)))
        start = end
    return result


