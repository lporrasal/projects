s = open('myfile.csv','r').read()

chars = ('$','%','^','*') # etc
for c in chars:
  s = ''.join( s.split(c) )

out_file = open('myfile_nospchar.csv','w')
out_file.write(s)
out_file.close()

s = open('myfile_nospchar.csv','r').read()

chars = (' ') # etc
for c in chars:
  s = '_'.join( s.split(c) )

out_file = open('myfile_nospcs.csv','w')
out_file.write(s)
out_file.close()