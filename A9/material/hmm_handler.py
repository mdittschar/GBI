import numpy as np


np.seterr(divide='ignore')


def log_data(n):
    """Optional function: Handles the case that a log recieves a 0. It returns minus infinity

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
        pass

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
        # viterbi_matrix_filling: fills up the matrices
        # viterbi_terminate: computes the final step of the recursion.
        # viterbi_traceback: returns the decoded sequence of states for the given sequence of symbols

        return decoded_states

    def viterbi_init(self):
        pass

    def viterbi_matrix_filling(self):
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
    hmm_object = HMMhandler()
    hmm_object.read_hmm("cpg.hmm")
    # TODO Parse fasta file for sequences
    # TODO For each sequence in the fasta file run the viterbi algorithm.
    # TODO Once decoded, print the original and the decoded sequences with the desired output format.


if __name__ == "__main__":
    main()
