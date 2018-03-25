#!/usr/bin/python3

from tree import Tree
import json
import sys
import os

def GetQ(counts_file_name = None):
    """Reads the counts file and returns the parameters of underlying CFG
    
    Reads the counts file given as argument and return maximul likelihood estimates for the 
    parameters. The counts file contains the rare words replaced by _RARE_ keyword. So, before
    using these parameters, any test data must be preprocessed to replace rare words by _RARE_

    """
    # Indexed by a tuple: (X, Y1, Y2) where X -> Y1 Y2 is the expansion
    q_binary_rules = dict()

    # Indexed by a tuple: (X, W) where X -> W is the expansion
    q_unary_rules = dict()

    # Indexed by a string: X where X is a non-terminal
    q_non_terminal = dict()

    # In the initial iteration, q_binary_rules, q_unary_rules store the counts
    # In the next iteration, these counts are divided by counts of the non-terminals
    ################ FIRST ITERATION #####################
    with open(counts_file_name, "r") as f_counts:
        for line in f_counts:
            tokens = line.strip().split()
            if(tokens[1] == "NONTERMINAL"):
                non_terminal = tokens[2]
                count_non_terminal = int(tokens[0])
                q_non_terminal[non_terminal] = count_non_terminal
             
            elif(tokens[1] == "BINARYRULE"):
                count_binary_rule = int(tokens[0])
                X, Y1, Y2 = tokens[2], tokens[3], tokens[4]
                q_binary_rules[(X, Y1, Y2)] = count_binary_rule
             
            elif(tokens[1] == "UNARYRULE"):
                count_unary_rule = int(tokens[0])
                X, W = tokens[2], tokens[3]
                q_unary_rules[(X, W)] = count_unary_rule
        
    #################### SECOND ITERATION #################
    # Divide the counts of binary rules by the counts of the respective non-terminals
    for X, Y1, Y2 in q_binary_rules:
        q_binary_rules[(X, Y1, Y2)] = q_binary_rules[(X, Y1, Y2)] / q_non_terminal[X]

    # Divide the counts of unary rules by the counts of the respective non-terminals
    for X, W in q_unary_rules:
        q_unary_rules[(X, W)] = q_unary_rules[(X, W)] / q_non_terminal[X]

    return q_binary_rules, q_unary_rules, list(q_non_terminal.keys())


def GetAllWords(counts_file_name = None):
    """Returns the list of all the words in the training data"""

    # Total number of words = 10024
    # Number of rare words = 8615
    # Total number of words after replacement = 10024 - 8615 + 1 (for _RARE_) = 1410
    all_words = set()
    with open(counts_file_name, "r") as f:
        for line in f:
            tokens = line.strip().split()
            if(tokens[1] == "UNARYRULE"):
                all_words.add(tokens[3])
    return all_words

def PreprocessRareWords(words = None, all_words = None):
    """Replace rare words with _RARE_"""
    for i in range(len(words)):
        if(words[i] not in all_words):
            words[i] = "_RARE_"

def getBinaryRulesFor(q_binary_rules, X):
    """Gets the binary rules of form X -> *"""
    binary_rules_with_X = [binary_rule for binary_rule in q_binary_rules if binary_rule[0] == X]
    return binary_rules_with_X    

def toJSONArray(bp = None, root_val = None, n = None):
    """Computes JSON representation of the underlying parse tree from bp dictionary"""

    # represents empty tree as of now
    json_array = []
    queue = list()
    queue.append([json_array, root_val, 0, n - 1])

    while(len(queue) != 0):

        # This means that a subtree rooted at `root_val` needs to be constructed
        subtree, root_val, l, r = queue[0]

        # Get the expansion rule and the split point
        expansion_rule, s = bp[(l, r, root_val)]

        # A binary rule is used for expansion
        if(len(expansion_rule) == 3):
            this_root_left_child = []
            this_root_left_child_root = expansion_rule[1]

            this_root_right_child = []
            this_root_right_child_root = expansion_rule[2]

            # Modify the subtree in place to include this root
            subtree.append(root_val)
            subtree.append(this_root_left_child)
            subtree.append(this_root_right_child)
            
            # queue the operations to construct left and right sub tree
            queue.append([this_root_left_child, this_root_left_child_root, l, s])
            queue.append([this_root_right_child, this_root_right_child_root, s + 1, r])

        # A unary rule is used for expansion
        elif(len(expansion_rule) == 2):
            word = expansion_rule[1]

            # Make the subtree X -> W
            subtree.append(root_val)
            subtree.append(word)

        else:
            raise Exception("toJSONArray: Invaluid expansion rule")
        
        # delete the current operation entry
        del queue[0]

    return json_array

def CKY(words, q_binary_rules, q_unary_rules, N):
    """Runs the dynamic programming based CKY on the given sentence
    The `words` has been preprocessed already to replace rare words with keyword rare.
    """
    # Indexed by start and end point and the non-terminal spanning the range
    pi = dict()
    bp = dict()

    #################### INITIALIZATION ##########################
    for i in range(len(words)):
        for X in N:
            rule = (X, words[i])
            if(rule in q_unary_rules):
                pi[(i, i, X)] = q_unary_rules[(X, words[i])]
                bp[(i, i, X)] = ((X, words[i]), -1)
            else:
                pi[(i, i , X)] = 0
    
    ############## MAIN LOOP OF THE ALGORITHM ##########
    n = len(words)

    # List of valid root non-terminals for the whole sentence
    valid_root_vals = []
    for l in range(2, n + 1):
        for i in range(0, n - l + 1):
            j = i + l - 1

            for X in N:
                pi[(i, j, X)] = 0

                # Stores the binary rule that gives the max probability
                max_binary_rule = None
                # Stores the split point that gives the max probability   
                max_s = None
                
                flag = False 
                
                for binary_rule in getBinaryRulesFor(q_binary_rules, X):
                    for s in range(i, j):
                        X_, Y, Z = binary_rule

                        # Should always evalulate to true
                        assert(X_ == X)

                        this_prob = q_binary_rules[binary_rule] * pi[(i, s, Y)] * pi[(s + 1, j, Z)]
                        assert(this_prob >= 0)
                        if(this_prob > pi[(i, j, X)]):
                                pi[(i, j, X)] = this_prob
                                max_binary_rule = binary_rule
                                max_s = s
                                flag = True

                # Store the back pointer information
                bp[(i, j, X)] = (max_binary_rule, max_s)
                
                # Append this non-terminal only if a max_binary_rule and a valid split point was
                # found
                if(i == 0 and j == n - 1 and flag == True):
                    valid_root_vals.append(X)
    
    # Handling the case where the sentence is a fragment
    root_val = None
    if(pi[(0, n - 1, 'S')] != 0):
        root_val = 'S'
    else:
        max_prob = 0
        for X in valid_root_vals:
            if(pi[(0, n - 1, X)] > max_prob):
                max_prob = pi[(0, n - 1, X)]
                root_val = X
    ##################### BUILD THE PARSE TREES OUT OF BACKPOINTERS ####################
    assert(root_val is not None)                
    parse_tree_as_array = toJSONArray(bp = bp, root_val = root_val, n = n)
    parse_tree_as_json = json.dumps(parse_tree_as_array)
    return parse_tree_as_json
         
def ParseTestData(test_data_file_name = None, counts_file_name = None, 
        test_predictions_file_name = None):
    """Computes the parse trees for the test data

    Reads the test data file line by line. Each line contains a single sentence. The sentence is
    preprocessed to replace rare words by _RARE_ (the parameters use the same keyword).

    """
    # Compute the name of the outut key file
    # For parse_test.dat, the output file name is parse_dev.key
    test_data_key_file_name = test_predictions_file_name 

    # Get the list of all words from the counts file and give it to PreprocessRareWords
    all_words = GetAllWords(counts_file_name = counts_file_name)
    
    # Calculate the parameters of the model
    q_binary_rules, q_unary_rules, N = GetQ(counts_file_name = counts_file_name)

#    # Sanity checks on probabilty
#    for binary_rule in q_binary_rules:
#        assert(q_binary_rules[binary_rule] > 0)
#    for unary_rule in q_unary_rules:
#        assert(q_unary_rules[unary_rule] > 0)

    with open(test_data_file_name, "r") as f_test_data_input, open(test_data_key_file_name, "w+") as f_test_data_output:
         for line in f_test_data_input:
            line = line.strip()
            words = line.split()
            
            # Replace rare words with _RARE_
            PreprocessRareWords(words = words, all_words = all_words)
            
            # Run the CKY on this sentence
            parse_tree_as_json = CKY(words, q_binary_rules, q_unary_rules, N)
            
            # Write the JSON to the prediction file
            f_test_data_output.write(parse_tree_as_json + "\n")

            

if __name__ == "__main__":
    
    # Parse the command line arguments    
    train_file_name = sys.argv[1]
    test_file_name = sys.argv[2]
    test_predictions_file_name = sys.argv[3]
    counts_file_name = "cfg_q5.counts"

    # Generate the counts file from the new train file: parse_train.RARE.dat
    cmd_counts_file_generation = "./count_cfg_freq.py %s > %s" % (
                                        train_file_name, counts_file_name)
    os.system(cmd_counts_file_generation)
        
    # Calculate parse trees for the test data
    ParseTestData(test_data_file_name = test_file_name, 
                    counts_file_name = counts_file_name,
                  test_predictions_file_name = test_predictions_file_name) 

    # Delete the counts file
    cmd_delete_counts_file = "rm -rf %s" % (counts_file_name)
    os.system(cmd_delete_counts_file)

    # Generate the evaluation results
    cmd_generate_evaluation_results = "python eval_parser.py parse_dev.key %s > q5_eval.txt" % (
                                    test_predictions_file_name)
    os.system(cmd_generate_evaluation_results)

