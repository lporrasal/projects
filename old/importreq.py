import requests
response = requests.get("http://www.pythonchallenge.com/pc/def/equality.html")

if response.status_code == 200:
	print response.content