import re

with open("Ch4.txt", "r") as fobj:
    text = fobj.read()

# Use word boundary anchors (\b) so only five-digit numbers are matched.
# Otherwise, 123456 would also be matched (and the match result would be 12345)!

output = re.findall(r'[A-Za-z]', text)
# Join the matches together

out_str = ",".join(output)
# Write them to a file, again using "with" so the file will be closed.

with open("output.txt", "w") as outp:
    outp.write(out_str)

print out_str 

#print "".join(re.findall("[^A-Z]+[A-Z]{3}([a-z])[A-Z]{3}[^A-Z]+", out_str))

#with findall(out_str):
#	findall
	