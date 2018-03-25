#!/usr/bin/python3

from tree import Tree
import json
import sys
import os

def GetParameters(counts_file_name = None):
    """Reads the counts file and returns the parameters of underlying CFG: q_binary_rules, q_unary_rules, N, all_words
    
    Reads the counts file given as argument and return maximul likelihood estimates for the 
    parameters. The counts file contains the rare words replaced by _RARE_ keyword. So, before
    using these parameters, any test data must be preprocessed to replace rare words by _RARE_

    """
    # Indexed by the non-terminal X: All the rules X -> Y Z are stored in q_binary_rules[X]
    # q_binary_rules[X] = (Y, Z, count of the binary rule)
    q_binary_rules = dict()

    # Indexed by the non-termianl X: All the rules X -> W are stored in q_unary_rules[X]
    q_unary_rules = dict()

    # Indexed by a string: X where X is a non-terminal. This is used when computing ML estimates
    q_non_terminal = dict()

    # list of all the words. This list includes rare words classified into the category: _RARE_
    all_words = set()

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

                # If no rule `X -> Y Z` has been seen, intialize q_binary_rules[X] to an empty dictionary
                if(X not in q_binary_rules):
                    q_binary_rules[X] = []
                    
                q_binary_rules[X].append((Y1, Y2, count_binary_rule))
 
            elif(tokens[1] == "UNARYRULE"):
                count_unary_rule = int(tokens[0])
                X, W = tokens[2], tokens[3]

                if(X not in q_unary_rules):
                    q_unary_rules[X] = dict()

                q_unary_rules[X][W] = count_unary_rule
                all_words.add(W)
        
    #################### SECOND ITERATION #################
    # Divide the counts of binary rules by the counts of the respective non-terminals
    for X in q_binary_rules:
        count_of_X = q_non_terminal[X]
        for index, (Y1, Y2, count) in enumerate(q_binary_rules[X]):
            q_binary_rules[X][index] = (Y1, Y2, count / count_of_X)

    # Divide the counts of unary rules by the counts of the respective non-terminals
    for X in q_unary_rules:
        count_of_X = q_non_terminal[X]
        for W in q_unary_rules[X]:
            q_unary_rules[X][W] = q_unary_rules[X][W] / count_of_X

    return q_binary_rules, q_unary_rules, list(q_non_terminal.keys()), all_words

def PreprocessRareWords(words = None, all_words = None):
    """Replace rare words with _RARE_"""
    for i in range(len(words)):
        if(words[i] not in all_words):
            words[i] = "_RARE_"

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
        # the X non-terminal must expand to a word. Those words are present in q_unary_rules
        for X in q_unary_rules:
            if(words[i] in q_unary_rules[X]):
                pi[(i, i, X)] = q_unary_rules[X][words[i]]

                # The -1 is appended because the bp below stores (rule, split_point)
                # The -1 here is useless but helps code be compatabile 
                bp[(i, i, X)] = ((X, words[i]), -1)
    
    ############## MAIN LOOP OF THE ALGORITHM ##########
    n = len(words)
    
    # Get the rules with binary expansions. These are used in one of the loops of CKY
    NT_with_binary_expansions = list(q_binary_rules.keys())

    # List of valid root non-terminals for the whole sentence
    valid_root_vals = []
    for l in range(2, n + 1):
        for i in range(0, n - l + 1):
            j = i + l - 1
            
            # The recursive rule expands this X into 2 non-terminals, thus, this X must have a 
            # binary expansion => it must belong to the q_binary_rules dictionary / NT_with_binary_expansions   
            for X in NT_with_binary_expansions:
                pi[(i, j, X)] = 0

                # Stores the binary rule that gives the max probability
                max_binary_rule = None
                
                # Stores the split point that gives the max probability   
                max_s = None
                
                flag = False 
                
                for Y, Z, count_of_binary_rule in q_binary_rules[X]:
                    for s in range(i, j):                      
                        # if these tuples are not present => pi = 0  for them => this_prob = 0
                        # and hence, this non-terminal can be skipped
                        if((i, s, Y) in pi and (s + 1, j, Z) in pi):
                            this_prob = count_of_binary_rule * pi[(i, s, Y)] * pi[(s + 1, j, Z)]
                            
                            if(this_prob > pi[(i, j, X)]):
                                    pi[(i, j, X)] = this_prob
                                    max_binary_rule = (X, Y, Z)
                                    max_s = s
                                    flag = True

                # Store the back pointer information
                bp[(i, j, X)] = (max_binary_rule, max_s)
                
                # Append this non-terminal only if a max_binary_rule and a valid split point was
                # found
                # if(i == 0 and j == n - 1 and flag == True):
                #   valid_root_vals.append(X)
    
    # Handling the case where the sentence is a fragment
    root_val = None
    if(pi[(0, n - 1, 'S')] != 0):
        root_val = 'S'
    else:
        max_prob = 0
        for X in NT_with_binary_expansions:
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
 
    # Calculate the parameters of the model: q_binary_rules, q_unary_rules and Non-terminals and 
    # the word list to be used by PreprocessRareWords
    ###################### VERIFIED THE `all-words` ##########################
    q_binary_rules, q_unary_rules, N, all_words = GetParameters(counts_file_name = counts_file_name)
    
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
    counts_file_name = "cfg_q6.counts"

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
    cmd_generate_evaluation_results = "python eval_parser.py parse_dev.key %s > q6_eval.txt" % (
                                    test_predictions_file_name)
    os.system(cmd_generate_evaluation_results)

