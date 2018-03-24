#!/usr/bin/python3

import json
import sys
import os

def getRareWords(file_name = None):
    """Return a list of rare words given the file name

    Returns the list of rare words after reading the given `counts` file 
    A word x is a rare word if Count(x) < 5 
    
    """
    f = open(file_name, "r")

    word_counts = dict()
    for line in f:
        
        # The lines with token UNARYRULE are the ones that contain termianl symbols
        tokens = line.strip().split(" ")
        if(tokens[1] == "UNARYRULE"):
            word =  tokens[3]

            # If the word has not been seen before, initialize its count to zero
            if(word not in word_counts):
                word_counts[word] = 0

            # Increment the count of the word by the current count
            this_count = int(tokens[0])
            word_counts[word] += this_count
     
    # A word is a rare word if its count < 5
    rare_words = [word for word in word_counts if word_counts[word] < 5]
    return rare_words


def findWordsInTree(tree = None):
    """Returns the words at the fringes of the `tree`

    Each tree is a list. If a tree is a Unary rule, it has a length == 2. If a list has a length
    == 3, it is a binary rule. If a tree is rooted at a binary rule, its first element is a Non-
    terminal and hence, not a word. So, the function calls itself recursively on the two childre
    n, unions the set of words returned.

    If the tree is a unary rule and has length == 2, the function returns the second element of t
    he list

    """
    # Unary Rule => second element is the word, return it as a set
    if(len(tree) == 2):
        return {tree[1]}
    elif(len(tree) == 3):
        words_in_left_subtree = findWordsInTree(tree[1])
        words_in_right_subtree = findWordsInTree(tree[2])
        return words_in_left_subtree | words_in_right_subtree
    else:
        raise Exception("findWordsInTree: Tree's length is not valid")

def replaceRareWordsInTree(tree = None, rare_words = None, rare_keyword = '_RARE_'):
    """Replaces the rare words in the given tree with the word specified by `rare_keyword`
    
    This function returns the tree after replacing the rare keywords in place.
    The words are always present in the Unary rules. For binary rules, the function calls itself.
    """

    # Unary rule => second element is the word, replace it in place
    if(len(tree) == 2):
        word = tree[1]
        if(word in rare_words):
            tree[1] = rare_keyword
    # Binary rule => call it recursively on the left and right sub-trees
    elif(len(tree) == 3):
        replaceRareWordsInTree(tree = tree[1], rare_words = rare_words, 
                rare_keyword = rare_keyword)         
        replaceRareWordsInTree(tree = tree[2], rare_words = rare_words, 
                rare_keyword = rare_keyword)                


def ReplaceRareWords(input_file_name = None, output_file_name = None, rare_words = None):
    """Read the `input_file_name`, replace the rare_words, save them into the `output_file_name`
    
    Reads the training file specified by `input_file_name`. Checks if a word is rare by checking 
    if it is present in rare_words (list), replaces it with a reserved keyword - _RARE_. Writes  
    the new data into `output_file_name`

    The old training and the new training files are written in JSON format.

    """
    # The output file name is opened for writing. Replaced if exists. Created if doesn't exist
    with open(input_file_name, "r") as f_input, open(output_file_name, "w+") as f_output:
        
        for line_number, line in enumerate(f_input):
            tree = json.loads(line)

            ######################################
            # Each `tree` is a list
            # Length of each tree == 3 (Asserted)
            ######################################
            
            # Replace the rare words in place
            replaceRareWordsInTree(tree = tree, rare_words = rare_words, rare_keyword = '_RARE_')
            assert(len(tree) == 3) 
            # Convert the new tree in json array
            tree_as_json = json.dumps(tree)
            
            # Write the json array to the new train file
            f_output.write(tree_as_json + "\n")

if __name__ == "__main__":
    
    # Parse the command line arguments to get the original train file name and the new 
    # train file name
    original_train_file = sys.argv[1]
    new_train_file = sys.argv[2]
    counts_file_name = "cfg_q4.counts"

    # Generate the cfg.counts file to get the list of the rare words
    # Generated everytime the code is run
    cmd_generate_counts_file = "./count_cfg_freq.py %s > %s" % (
                                original_train_file, counts_file_name)
    os.system(cmd_generate_counts_file)

    # Number of rare words found = 8615
    # Total number of words = 10024
    rare_words = getRareWords(file_name = counts_file_name)
    
    ReplaceRareWords(input_file_name = original_train_file, output_file_name = new_train_file
            , rare_words = rare_words)
    
    # Delete the counts file
    cmd_delete_counts_file = "rm -rf %s" % (counts_file_name)
    os.system(cmd_delete_counts_file)
