from websource import websource
import requests

#def get_country():
response = requests.get (websource)
if response.status_code == 200: 
	print response.content
