#from curses import keyname
import sys
import random  # needed in simulate_previous_generation
import itertools
import math
import numpy as np
import matplotlib.pyplot as plt

class PopGenSimulator:
    def __init__(self, size):
        self.size = size  # population size

    def simulate_previous_generation(self, current_generation):
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
        # set random seed to make the results reproducalbe       
        random.seed(1)
        np.random.seed(50)
        # make a list out of the current generation for workability
        list_cur_gen = list(current_generation)

        total_waiting_time = 0
        # while the waiting time is shorter than a generation, add up
        while total_waiting_time <= 1/(2*len(list_cur_gen)) and len(list_cur_gen) > 1:
            # compute waiting time according to formula
            waiting_time = - (math.comb(len(list_cur_gen), 2)**(-1))*np.log(np.random.uniform(0,1))
            gen_combs = [i for i in itertools.combinations(list_cur_gen, 2)]
            #choose a random pair that coalesces
            random.seed(0)
            choice = random.choice(gen_combs)
           # print('Choice', choice)
            parent = choice[0]
            child = choice[1]
            # if parent not in parent_directory:
            #     parent_directory[parent] = [child]
            # elif isinstance(parent_directory[parent], list):                
            #     parent_directory[parent] = parent_directory[parent] + [child]
            # else:
            #     parent_directory[parent] = [parent_directory[parent]] + [child]
            # waiting_time =- (math.comb(len(list_cur_gen), 2)**(-1))*np.log(np.random.uniform(0,1))
            list_cur_gen = [i for i in list_cur_gen if i != child]
            total_waiting_time = total_waiting_time + waiting_time
                
            
        #print("Parent directory: ", parent_directory)
        # cur_string = "".join([k if k in parent_directory else "-" for k in list(current_generation) ])
        # #print(cur_string)
        # parent_directory = dict(sorted(parent_directory.items()))
        # print("Parent dict:", parent_directory)
        # output_nextgen = "".join([k for k in parent_directory])

        #print("Output next gen: ", output_nextgen)
        return total_waiting_time, list_cur_gen

    def is_MRCA_of_all_found(self, current_generation):
        """
        have we found the MRCA of all individuals of the initial population?
        :param current_generation: int
        :return:
        """
        # find out the length of current genneration excluding dashes
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
        genes = "abcdefghijklmnopqrstuvwxy"
        # take first size characters of genes
        genes = genes[:size]
        return genes

    def string_output(self, current_generation, original_generation, gen_no):
        """For a given generation, print out the string 

        Args:
            current_generation (list): list of genes currently in the gene pool
            original_generation (string): original genes in the first considered generation
            gen_no (int): number of generations back
        """
        # replace every letter of the original generation that isn't in the current generation with a dash
        cur_string = "".join([k if k in current_generation else "-" for k in list(original_generation) ])
        # makes the string output look more uniform by placing a space before the 0
        if gen_no == 0:
            gen_no = " 0"
        print(f"{gen_no}  {cur_string}")
    
    def plot_treesizes(self, sizes, median_arr):
        """Given an array of tree sizes, make a scatter plot of those vs. predicted viszes

        Args:
            sizes (np.array of ints): sizes to consider
            median_arr (np.array): array of medians
        """
        plt.scatter(sizes, median_arr, label="Random trees")
        plt.scatter(sizes, [2*(1-1/n) for n in sizes], label="Expected trees")
        plt.title("Randomly generated vs. Expected tree sizes - Median of 3 runs")
        plt.xlabel("Number of genes")
        plt.ylabel("Tree height")
        plt.ylim(0,None)
        plt.legend()
        plt.savefig("auckenthaler_dittschar_median_vs_expected_trees.png")
        plt.show()
        


def main():
    # try:
    sizes = sys.argv[1:]
    sizes = list(map(int, sizes))
    runs = np.arange(3)
    median_arr = []
    for size in sizes:
        for run in runs:
            # size = int(size)
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
            wait_arr = []
            total_waiting_time = 0
            while not pop_gen_simulator.is_MRCA_of_all_found(next_gen):
                waiting_time, next_gen = pop_gen_simulator.simulate_previous_generation(next_gen)
                gen_no = gen_no -1
                total_waiting_time = total_waiting_time + waiting_time
                pop_gen_simulator.string_output(next_gen, genes, gen_no)
            print()
            wait_arr = np.append(wait_arr,total_waiting_time)
        median_arr = np.append(median_arr, np.median(wait_arr))
    pop_gen_simulator.plot_treesizes(sizes= sizes, median_arr = median_arr)
  


if __name__ == "__main__":
    main()
