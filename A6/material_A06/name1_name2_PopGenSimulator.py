#from curses import keyname
import sys
import random  # needed in simulate_previous_generation
import itertools
import math
import numpy as np

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
        parent_directory = dict()
        print("Parent dict: ", parent_directory)
        for i in list(current_generation):
            subset_generation = [s for s in list(current_generation) if s != i]
            parent = random.choice(subset_generation)
            print("Parent is: ", parent) 
            
            if parent not in parent_directory:
                parent_directory[parent] = i
            elif type(parent_directory[parent]) == "<class 'list'>":
                parent_directory[parent] = parent_directory[parent] +[i]
            else:
                parent_directory[parent] = [parent_directory[parent]] +[i]
            print("Type: ", type(parent_directory[parent]))
            print("Parent_directory: ", parent_directory)
        new_gen = list()
        for i in list(current_generation):
            if i in parent_directory:
                new_gen = new_gen+ [i]
            else:
                new_gen = new_gen + ["-"]
        print("New gen: ",new_gen)
        waiting_time =- math.comb(len(current_generation), 2)**(-1)*np.log(np.random.uniform(0,1))

        return waiting_time

    def is_MRCA_of_all_found(self, current_generation):
        """
        have we found the MRCA of all individuals of the initial population?
        :param current_generation: int
        :return:
        """
        # PLEASE IMPLEMENT
        pass

    def make_initial_generation(self, size):
        """
        makes the initial generation, one letter per individual
        :param size: int
        :return:
        """
        # PLEASE IMPLEMENT
        pass


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
    genes = "abcdefghijklm"
    genes = genes[:size]
    print(genes)
    waiting_time = pop_gen_simulator.simulate_previous_generation(current_generation = genes)
    genes = ["a","b","c","d","e","f","g", "h", "i", "j", "k"]
    genes = genes[:size]
    print("Waiting time: ", waiting_time)
    


if __name__ == "__main__":
    main()
