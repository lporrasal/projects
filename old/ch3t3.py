import re

with open("Ch4.txt", "r") as fobj:
    text = fobj.read()

  print .join(re.findall("[^A-Z]+[A-Z]{3}([a-z])[A-Z]{3}[^A-Z]+", text))