import os
import itertools
import collections
import pprint
import sys
import os.path
import glob
import re

mypath = os.path.abspath(os.getcwd())  # get path of current dir

print("Directory path detected : ", mypath)

logfile = glob.glob("*log*.txt")  # getting the log filename
print(logfile)
file_path = os.path.join(mypath, logfile[0])


print("\nFile path detected : ", file_path)

with open(logfile[0], "r") as f:
    i = 0
    for line in f:
        if "-+" in line:
            nextline = next(f)
            i = i + 1

            nextlinearray = (
                nextline.split()
            )  # splitting the first row in different values
            bind_aff = nextlinearray[1]  # getting the binding affinity of first pose

            with open("output.txt", "a") as myfile:
                print(i, " : ", bind_aff, end="\n", file=myfile)

    print("Done! The result is provided in the output.txt file.")
