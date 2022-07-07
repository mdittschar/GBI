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
        # implemented in viterbi_matrix_filling()
        pass

    def get_emission_probability(self, state, symbol):
        """Optional function: Given a state and a symbol from the current HMM, returns the emission probability of the symbol under the state.

        Args:
            state (str): current state
            symbol (str): symbol that is emitted

        Returns:
            float: probability of emitting symbol under state
        """
        # implemented in viterbi_matrix_filling()
        pass

    def runViterbi(self, sequence):
        """ For the given sequence of symbols, return the decoded sequence of symbols using the Viterbi algorithm.

        Args:
            sequence (list): a list of symbols to decode using the Viterbi algorithm.
            id: 

        Returns:
            path (str): decoded sequence of states for the given sequence.
        """

        # TODO: implement the following functions.
        # You may need to adapt the input parameters for each function.

        # viterbi_init: initializaion of the Viterbi algorithm
        v_mat = HMMhandler.viterbi_init(self, sequence)
        # viterbi_matrix_filling: fills up the matrices            
        # viterbi_terminate: computes the final step of the recursion.
        v_mat, ptr, states_dict = HMMhandler.viterbi_matrix_filling(self, v_mat,  sequence)
        # viterbi_traceback: returns the decoded sequence of states for the given sequence of symbols
        path = HMMhandler.viterbi_get_path(self, v_mat, ptr, states_dict)

        return path

    def viterbi_init(self, sequence):   
        """Generate the initial Viterbi Matrix

        Args:
            sequence (list): sequence to compute Viterbi matrix for

        Returns:
            vmat: (state_nos, sequence_chars + 2): Viterbi matrix with -inf for all entries except (0,0)
        """
        v_mat = np.zeros((self.state_number,len(sequence) + 2))        
        v_mat[:,:] = log_data(0)
        v_mat[0,0] = log_data(1)
        return v_mat

    def viterbi_matrix_filling(self, v_mat, sequence):
        """Fill the Viterbi matrix

        Args:
            vmat: (np.array) Shape: (state_nos, sequence_chars + 2): Viterbi matrix with -inf for all entries except (0,0)
            sequence (list): sequence to compute Viterbi matrix for

        Returns:
            v_mat: (np.array) Shape: (state_nos, sequence_chars + 2): Filled Viterbi matrix
            states_dict: dictionary: Dictionary with states as keys and index as values
            ptr: (np.array) shape: (state_nos, sequence_chars + 2): Traceback matrix with 1 if the previous state was positive,
            0 if it was negative
            states_dict: dict : Dictionary with states and corresponding indices
        """
        # get states into an ordered dict
        states_dict = OrderedDict()
        vis_states_dict = OrderedDict()
        for i, val in enumerate(self.state_names): 
            states_dict[val] = i 
        for i, val in enumerate(self.symbol_names):   
            vis_states_dict[val] = i
        ptr = np.full_like(a=v_mat, fill_value=-np.inf)
        for i, s in enumerate(sequence):
            l = states_dict[s.lower()]
            h = states_dict[s]
            # in the beginning, go from entry (0,0) to the respective states
            if i == 0:
                v_mat[l, i + 1] = max(v_mat[0, i]+ log_data(self.transm_matrix[0, l]),v_mat[0, i]+log_data(self.transm_matrix[h, l]))
                v_mat[h, i + 1] = max(v_mat[0, i]+ log_data(self.transm_matrix[0, h]),v_mat[0, i]+log_data(self.transm_matrix[l, h]))
                # update pointer
                ptr[l, i+1] = 0
            # after that, insert computed probabilities to the current letter
            else:
                l_to_l = v_mat[prev_l, i]+ log_data(self.emission_matrix[l,h-1]*self.transm_matrix[prev_l, l])
                h_to_l = v_mat[prev_h, i]+ log_data(self.emission_matrix[l,h-1]*self.transm_matrix[prev_h, l])
                h_to_h = v_mat[prev_h, i]+ log_data(self.emission_matrix[h,h-1]*self.transm_matrix[prev_h, h])
                l_to_h = v_mat[prev_l, i]+ log_data(self.emission_matrix[l,h-1]*self.transm_matrix[prev_l, h])
                
                v_mat[l, i + 1] = max(l_to_l, h_to_l)
                v_mat[h, i + 1] = max(h_to_h, l_to_h)
                # also update pointer on the right positions
                ptr[l, i + 1] = np.argmax([l_to_l, h_to_l])
                ptr[h, i + 1] = np.argmax([l_to_h, h_to_h])
                
                # terminate (outside function)
                if i == len(sequence) - 1:
                    l_to_0 = v_mat[l, i+1]+ log_data(self.transm_matrix[l, 0])
                    h_to_0= v_mat[h, i+1]+ log_data(self.transm_matrix[h, 0])
                    v_mat[0, i + 2] = max(l_to_0,h_to_0)
                    # update pointer
                    ptr[0, i + 2] = np.argmax([l_to_0, h_to_0])

            prev_l = l
            prev_h = h
        return v_mat, ptr,  states_dict

    def viterbi_terminate(self):
        pass

    def viterbi_get_path(self, v_mat, ptr, states_dict):
        """Get path string from Viterbi matrix

        Args:
            v_mat: (np.array) Shape: (state_nos, sequence_chars + 2): Filled Viterbi matrix
            states_dict: dictionary: Dictionary with states as keys and index as values

        Returns:
            path (str): Traceback path string
        """
        states_arr = []
        # for every column, get letter that corresponds to highest entry in upper or lower half of column (dependent on pointer)
        idx = 0
        for c in np.arange(1,v_mat.shape[1]):
            if ptr[idx,-c] == 1:
                # get key by value
                # Adapted from: https://stackoverflow.com/questions/8023306/get-key-by-value-in-dictionary
                val = list(states_dict.keys())[list(states_dict.values()).index(np.argmax(v_mat[:5,-c-1]))]               
                idx = np.argmax(ptr[:5,-c-1]).astype(int)
            elif ptr[idx, -c] == 0:
                val = list(states_dict.keys())[list(states_dict.values()).index(5 + np.argmax(v_mat[5:,-c-1]))]
                idx = 5 + np.argmax(ptr[5:,-c-1]).astype(int)

            states_arr = np.append(states_arr, val)
        path = "".join(states_arr[::-1])[1:]
        return path


def prettyPrinting(input, id, decoded, filename):
    """ Formats the given sequences to match the desired output format and saves it.

    Args:
        input (str): the original input sequences
        decoded (str): the sequence with the decoded states
    """
    print("Symbols: ")
    print("Viterbi: ")
    #Code from our assignment 2
    pair_no = 60
    input = "".join(input)
    seq_lens = len(input)
    i = 0
    with open(filename, 'w') as file_out:
        file_out.write(f"ID of Sequence: {id}\nSymbols: \nViterbi:\n")
        while i*pair_no < seq_lens:
            # assign strings to be visualised
            # pair_no is number of pairs to visualise in one row
            display_seq0 = input[i*pair_no:(i+1)*pair_no]
            display_seq1 = decoded[i*pair_no:(i+1)*pair_no]
            if i*pair_no < seq_lens - pair_no:
                empty = (pair_no - len(str(seq_lens - pair_no*i)))* " "
                # numbers at the ends of the lines to show number of pairs
                alignment_nos = "".join([str(i*pair_no + 1), empty, str((i+1)*pair_no)])
            else:
                # visualisation if there's an incomplete line left
                empty = " "*(pair_no)
                alignment_nos = "".join([str(i*pair_no + 1), empty, str(seq_lens)])
            # print the combined information
            print(f"\n{alignment_nos}\n{display_seq0}\n{display_seq1}")
            file_out.write(f"\n{alignment_nos}\n{display_seq0}\n{display_seq1}")
            i = i+1

    print()
    
    pass

def get_seqs_and_ids():
    """From a fasta file, get sequences and ids

    Returns:
        sequences (list): sequences of fasta file
        ids (np.array): ids of sequences
    """
    argv = sys.argv[1:]
    # get a tuple of the arguments to use
    opts, _ = getopt.getopt(argv, "a:s:", ['hmm','sequence'])
    sequences= [None]*4
    ids = []
    i = 0
    for opt, arg in opts:
        if opt in ["-s", "--sequence"]:        
            for record in SeqIO.parse(arg, "fasta"):
                    id = record.id
                    se = record.seq
                    sequences[i] = list(str(se))
                    ids = np.append(ids, id)
                    i = i + 1
    return sequences, ids



def main():
    sequences, ids = get_seqs_and_ids()
    # sequences = [list("CGCG")]
    # ids = "1"
    hmm_object = HMMhandler()
    hmm_object.read_hmm()#"cpg.hmm")
    i = 1
    for sequence, id in zip(sequences, ids): 
        path  = hmm_object.runViterbi(sequence)
        prettyPrinting(sequence, id, path, f"auckenthaler_dittschar_viterbi_{id}.txt")
        i = i + 1
    # TODO Parse fasta file for sequences
    # TODO For each sequence in the fasta file run the viterbi algorithm.
    # TODO Once decoded, print the original and the decoded sequences with the desired output format.


if __name__ == "__main__":
    main()
