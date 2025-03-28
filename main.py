from read_data import *
from hypermotif import *




if __name__ == '__main__':
    #settings
    data_name='iJO1366'  # dataset name

    order=3            #ternary patterns
    P = [1, 1, 1]     # sampling probability

    # order = 4  # quaternary patterns
    # P = [1, 1, 1, 1]  # sampling probability

    num_of_processe = 128  # number of parallel computing processes
    num_of_tasks = num_of_processe * 5  # number of slices


    #census
    N,cooked_edges_list= read_data(data_name)   #read hypergraph
    H = hypergraph(N=N, edges=cooked_edges_list, order=order)
    total_C = multiprocess_run_census(H, P, num_of_processe, num_of_tasks)

    #examine the result
    random_motif = random.choice(list(total_C.keys()))
    print("pattern:",s_to_t(random_motif),"frequency:",total_C[random_motif])

