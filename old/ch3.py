import re

from sys import argv

filename = argv 

data = open("/Users/lauraporras/projects/ch3.txt") 

print (re.findall("[A-Za-z]", data))