filetoclean = raw_input ("file:")

s = open(filetoclean,'r').read()

chars = ('$','%','^','*') # etc
for c in chars:
  s = ''.join( s.split(c) )

out_file = open('nospchar.csv','w')
out_file.write(s)
out_file.close()

s = open('nospchar.csv','r').read()

chars = (' ') # etc
for c in chars:
  s = '_'.join( s.split(c) )

out_file = open('nospcs.csv','w')
out_file.write(s)
out_file.close()