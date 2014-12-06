import re
import requests
import difflib

response = requests.get("https://www.xoom.com/mexico/fees-fx")

if response.status_code == 200:
	#response.content 
		data = response.content 
		#data.get = ('USD =')
		#def exchange_rate
		#exchange = (re.findall('USD = %%%%%%' , data))
		exchange = (get_field('USD = %%%%%%' , data))
		print exchange 
		#exchange_rate = 'USD = %', % exchange_rate
		#print 'exchange is' '%'
		#exchange2 =+ 1
		#rate = 
		#rate = (USD %%%%%%)
		#print "USD : %" % data.get('USD = %')



#print 'Rate:', exchange_rate


#print data 

#print exchange2
	#print response.content