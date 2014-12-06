def word_count(path):


	counts = {
		'lines' : 0,
		'words' : 0,
		'chars' : 0
	}

	with open(path, 'r') as f:
		for line in f:
			counts['lines'] += 1
			counts['chars'] += len(line)
			counts['words'] += len(line.split())

	return counts

if __name__ == '__main__':     
	print(word_count('ex15_sample.txt'))
