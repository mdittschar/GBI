#from curses import keyname
import sys
import random  # needed in simulate_previous_generation
import itertools
import math
import numpy as np

class PopGenSimulator:
    def __init__(self, size):
        self.size = size  # population size

    def simulate_previous_generation(self, current_generation, prev_waiting_time, gen_no):
        """
        given a string of letters, each letter representing one individual,
        randomly generate previous generation. For this use use the exponential distribution exp(k(k-1)/2) for the waiting times.
        Record them.
        Then randomly choose a pair which coalesces.
        In previous generation, use lexicographically first letter of all
        children as label of parent, or -, if parent has no child
        :param current_generation: String
        :return: previous_generation
        """
        
        list_cur_gen = list(current_generation)
        

        total_waiting_time = prev_waiting_time
        parent_directory = dict()
        waiting_time =- (math.comb(len(list_cur_gen), 2)**(-1))*np.log(np.random.uniform(0,1))
        if len(list_cur_gen) == 2:
            used_waiting_time = 0
        else: 
            used_waiting_time = waiting_time

        while total_waiting_time + used_waiting_time <= -gen_no + 1 and len(list_cur_gen) > 1:
            gen_combs = [i for i in itertools.combinations(list_cur_gen, 2)]
            print("Gen combs: ", gen_combs)
            
            total_waiting_time = total_waiting_time + waiting_time
            #print("Current waiting time: ", total_waiting_time)
            choice = random.choice(gen_combs)
            parent = choice[0]
            child = choice[1]

            #print("Parent directory is: ", parent_directory)
            if parent not in parent_directory:
                parent_directory[parent] = [child]
            elif isinstance(parent_directory[parent], list):
                
                parent_directory[parent] = parent_directory[parent] + [child]
            else:
                parent_directory[parent] = [parent_directory[parent]] + [child]
            waiting_time =- (math.comb(len(list_cur_gen), 2)**(-1))*np.log(np.random.uniform(0,1))
            list_cur_gen = [i for i in list_cur_gen if i != child and i != parent]
            print("Total waiting time: ", total_waiting_time)
            print("Current waiting time: ", used_waiting_time)
                
            
        #print("Parent directory: ", parent_directory)
        cur_string = "".join([k if k in parent_directory else "-" for k in list(current_generation) ])
        #print(cur_string)
        parent_directory = dict(sorted(parent_directory.items()))
        print("Parent dict:", parent_directory)
        output_nextgen = "".join([k for k in parent_directory])

        #print("Output next gen: ", output_nextgen)
        return total_waiting_time, output_nextgen

    def is_MRCA_of_all_found(self, current_generation):
        """
        have we found the MRCA of all individuals of the initial population?
        :param current_generation: int
        :return:
        """
        if len([i for i in current_generation if i != "-"]) <= 1:
            bool = True
        else:
            bool = False
        return bool

    def make_initial_generation(self, size):
        """
        makes the initial generation, one letter per individual
        :param size: int
        :return:
        """
        genes = "abcdefghijklm"
        genes = genes[:size]
        return genes

    def string_output(self, current_generation, original_generation, gen_no):
        cur_string = "".join([k if k in current_generation else "-" for k in list(original_generation) ])
        if gen_no == 0:
            gen_no = " 0"
        print(f"{gen_no}  {cur_string}")



def main():
    # try:
    size = int(sys.argv[1])
    pop_gen_simulator = PopGenSimulator(size)
    """
        PLEASE IMPLEMENT
            1. generate the initial population

            2. then run the simulator, print out each generation in this format:

                0 abcdefgh
            -1 ab-d-fg-
            -2 a-f---gd
            -3 g-f--a--
            -4 --f---a-
            -5 ---a----

            3. stop when you reach the MRCA of all existing individuals
        """
    # except:
    #     print('run name1_name2_PopGenSimulator.py 4')

    genes = pop_gen_simulator.make_initial_generation(size)
    next_gen = genes
    gen_no = 0
    pop_gen_simulator.string_output(next_gen, genes, gen_no)
    waiting_time = 0
    while not pop_gen_simulator.is_MRCA_of_all_found(next_gen):
        waiting_time, next_gen = pop_gen_simulator.simulate_previous_generation(next_gen, waiting_time, gen_no)
        print("Next gen: ", next_gen)
        gen_no = gen_no -1
        pop_gen_simulator.string_output(next_gen, genes, gen_no)


    


if __name__ == "__main__":
    main()
