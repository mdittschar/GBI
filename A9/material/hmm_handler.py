import numpy as np
import getopt
import sys
from io import StringIO
import pandas as pd
from Bio import SeqIO
from collections import OrderedDict

np.seterr(divide='ignore')


def log_data(n):
    """Optional function: Handles the case that a log receives a 0. It returns minus infinity

    Args:
        n (float)

    Returns:
        float or minus inf: log of n or minus infinity if 0 is given
    """
    return -np.inf if n == 0 else np.log(n)


class HMMhandler():
    """A class for handling Hidden Markov Models
    """

    def __init__(self, state_names=[], symbol_names=[], transm_matrix=np.matrix([], []), emission_matrix=np.matrix([], [])):
        """Constructor for the HMM handler class. 
        We recommend saving the parameters of the HMM model as properties of this clas. 
        It is not a must.
        """
        self.state_number = len(state_names)
        self.symbols_number = len(symbol_names)
        self.state_names = state_names
        self.symbol_names = symbol_names
        self.transm_matrix = transm_matrix
        self.emission_matrix = emission_matrix

    def read_hmm(self):
        """Parses an HMM from the fixed format. See script P. 163.
        Relies on the file having the same order and structure from the shown structure.
        """
        argv = sys.argv[1:]
        # get a tuple of the arguments to use
        opts, _ = getopt.getopt(argv, "a:s:", ['hmm', 'sequence'])


        for opt, arg in opts:
            if opt in ["-a", "--hmm"]:
                with open(arg) as f:
                    lines = f.readlines()
                # get the hidden states
                states = lines[3]
                # get the visible states
                vis_states = lines[7]
                vis_states = list(vis_states)[::2]
                vis_states_no = len(vis_states)
                # every second element in string is a state
                states = list(states)[::2]
                states_no = len(states)
                # get the lines that have the emission matrix
                emissions = lines[-states_no -1:-1]
                
                e_mat = np.zeros((states_no, vis_states_no))
                # convert from a text matrix to a numpy array
                for i, e in enumerate(emissions): 
                    #print("E: ",e)
                    e_mat[i, :] = np.squeeze(np.array(pd.read_csv(StringIO(e[:-1]),header=None,  sep=" "))).astype(float)
                #print("E_mat: ", e_mat)
                # get the lines that have the transition matrix
                transitions = lines[-2*states_no -2 :-states_no-2]
                t_mat = np.zeros((states_no , states_no))
                # convert from a text matrix to a numpy array
                for i, t in enumerate(transitions):                     
                    t_mat[i, :] = np.squeeze(np.array(pd.read_csv(StringIO(t[:-2]),header=None,  sep=" "))).astype(float)
                
        self.state_number = states_no
        self.symbols_number = vis_states_no
        self.state_names = states
        self.symbol_names = vis_states
        self.transm_matrix = t_mat
        self.emission_matrix = e_mat         
        

    def get_transition_probability(self, start_state, to_state):
        """Optinal function: Given states from the current HMM, returns the transition probability.

        Args:
            start_state (str): starting state
            to_state (str): ending state

        Returns:
            float: probability of going from start_state to to_state
        """
        pass

    def get_emission_probability(self, state, symbol):
        """Optional function: Given a state and a symbol from the current HMM, returns the emission probability of the symbol under the state.

        Args:
            state (str): current state
            symbol (str): symbol that is emitted

        Returns:
            float: probability of emitting symbol under state
        """
        pass

    def runViterbi(self, sequence):
        """ For the given sequence of symbols, return the decoded sequence of symbols using the Viterbi algorithm.

        Args:
            sequence (str): a string of symbols to decode using the Viterbi algorithm.

        Returns:
            str: decoded sequence of states for the given sequence.
        """

        decoded_states = ""
        # TODO: implement the following functions.
        # You may need to adapt the input parameters for each function.

        # viterbi_init: initializaion of the Viterbi algorithm
        v_mat = HMMhandler.viterbi_init(self, sequence)
        # viterbi_matrix_filling: fills up the matrices
        v_mat = HMMhandler.viterbi_matrix_filling(self, v_mat,  sequence)
        # viterbi_terminate: computes the final step of the recursion.
        # viterbi_traceback: returns the decoded sequence of states for the given sequence of symbols

        return decoded_states

    def viterbi_init(self, sequence):        
        v_mat = np.zeros((self.state_number,len(sequence) + 2))
        v_mat[0, 0] = log_data(1)
        return v_mat

    def viterbi_matrix_filling(self, v_mat, sequence):
        # get states into an ordered dict
        states_dict = OrderedDict()
        vis_states_dict = OrderedDict()
        for i, val in enumerate(self.state_names): 
            states_dict[val] = i 
        for i, val in enumerate(self.symbol_names):   
            vis_states_dict[val] = i
        for i, s in enumerate(sequence):
            l = states_dict[s.lower()]
            h = states_dict[s]
            if i == 0:
                v_mat[l, i + 1] = max(v_mat[0, i]+log_data(self.transm_matrix[0, l]),v_mat[0, i]+log_data(self.transm_matrix[h, l]))
                v_mat[h, i + 1] = max(v_mat[0, i]+ log_data(self.transm_matrix[0, h]),v_mat[0, i]+log_data(self.transm_matrix[l, h]))
            else:
                v_mat[l, i + 1] = max(v_mat[prev_l, i]+ log_data(self.transm_matrix[prev_l, l]),v_mat[prev_h, i]+ log_data(self.transm_matrix[prev_h, l]))
                v_mat[h, i + 1] = max(v_mat[prev_h, i]+ log_data(self.transm_matrix[prev_h, h]),v_mat[prev_l, i]+ log_data(self.transm_matrix[prev_l, h]))
                if i == len(sequence) -1:
                    v_mat[l, i + 2] = max(v_mat[l, i+1]+ log_data(self.transm_matrix[l, 0]),v_mat[h, i+1]+ log_data(self.transm_matrix[h, 0]))
                    v_mat[h, i + 2] = max(v_mat[h, i+1]+ log_data(self.transm_matrix[h, 0]),v_mat[l, i+1]+ log_data(self.transm_matrix[l, 0]))

            prev_l = l
            prev_h = h
            # v_mat[:,i + 1] = self.emission_matrix[:, vis_states_dict[s]]

        for line in v_mat.T[-10:]:
            print(line)
        #print("V_mat: ", v_mat)
        pass

    def viterbi_terminate(self):
        pass

    def viterbi_get_path(self):
        pass


def prettyPrinting(input, decoded):
    """ Formats the given sequences to match the desired output format.

    Args:
        input (str): the original input sequences
        decoded (str): the sequence with the decoded states
    """
    pass


def main():
    argv = sys.argv[1:]
    # get a tuple of the arguments to use
    opts, _ = getopt.getopt(argv, "a:s:", ['hmm','sequence'])
    sequences= [None]*4
    i = 0
    for opt, arg in opts:
        if opt in ["-s", "--sequence"]:        
            for record in SeqIO.parse(arg, "fasta"):
                    se = record.seq
                    sequences[i] = list(str(se))
                    i = i + 1
    hmm_object = HMMhandler()
    hmm_object.read_hmm()#"cpg.hmm")
    hmm_object.runViterbi(sequences[0])
    # TODO Parse fasta file for sequences
    # TODO For each sequence in the fasta file run the viterbi algorithm.
    # TODO Once decoded, print the original and the decoded sequences with the desired output format.


if __name__ == "__main__":
    main()
