import requests

site = raw_input("Direccion web:"), 

response = requests.get (site)
if response.status_code == 200:
   print response.content