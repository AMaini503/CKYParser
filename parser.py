#!/usr/bin/python

import sys
import os

def start():
    """Entry point for all the programs"""
    question_number = sys.argv[1]
    assert(question_number == "q4" or
           question_number == "q5" or
           question_number == "q6")

    if(question_number == "q4"):
        # SAMPLE USAGE: python parser.py q4 parse_train.dat parse_train.RARE.dat    

        # Original Train file
        original_train_file = sys.argv[2]

        # Train file with rare words replaced with _RARE_
        new_train_file = sys.argv[3]

        # Run the q4.py to generate the new training file
        cmd = "./q4.py %s %s" % (original_train_file, new_train_file)
        os.system(cmd)

    elif(question_number == "q5"):
        # SAMPLE USAGE: python parser.py q5 parse_train.RARE.dat parse_dev.dat q5_prediction_file
        
        train_file_name = sys.argv[2]
        test_file_name = sys.argv[3]
        test_predictions_file_name = sys.argv[4]

        # Run the q5 to produce predictions and evaluation results
        cmd =  "./q5.py %s %s %s" % (train_file_name, test_file_name, test_predictions_file_name)
        os.system(cmd)
    elif(question_number == "q6"):
        # SAMPLE USAGE: 
        # python parser.py q6 parse_train_vert.RARE.dat parse_dev.dat q6_prediction_file
        
        train_file_name = sys.argv[2]
        test_file_name = sys.argv[3]
        test_predictions_file_name = sys.argv[4]
        
        # Run the q6 to produce predictions and evaluation results
        cmd = "./q6.py %s %s %s" % (train_file_name, test_file_name, test_predictions_file_name)
        os.system(cmd)     
if __name__ == "__main__":
    start()


