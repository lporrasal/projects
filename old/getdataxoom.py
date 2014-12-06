import requests
import csv

website = raw_input ("sitio:")

response = requests.get(website)
if response.status_code == 200:
   
   with open (response.content) as csvfile:
   	csvfile.seek ("USD")
   	print csvfile.seek 