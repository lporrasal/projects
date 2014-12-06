import csv

with open("nba_salary.csv", 'rU') as data:
        reader = csv.DictReader(data)
        for row in reader:
            print (row ["SALARY"])
            #print (row["PER"])
