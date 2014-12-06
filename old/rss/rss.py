

#requirements.txt

#feedparser==5.1.3
#requests==2.4.1
#wsgiref==0.1.2

#rss.py
#!/usr/bin/env python
 
##########################################################################
## Imports
##########################################################################
 
import os
import json
import requests
import feedparser
from datetime import datetime
 
##########################################################################
## Module Constants
##########################################################################
 
BASEDIR = os.path.abspath(os.path.dirname(__file__))
FIXTURES = os.path.join(BASEDIR, 'fixtures')
FEEDS = os.path.join(FIXTURES, 'feeds.json')
DOWNLOADS = os.path.join(FIXTURES, 'downloads')
 
##########################################################################
## Functions
##########################################################################
 
def feeds(path=FEEDS):
#"""
#Opens up our JSON list of RSS feeds and returns it as a Python list
#"""

with open(path, 'r') as data:
return json.load(data)
 
def ingest(path=FEEDS, downloads=DOWNLOADS):
"""
Ingest RSS feeds from a list of feeds to a downloads directory.
"""

for feed in feeds(path):
download_dir = os.path.join(downloads, feed["id"])
if not os.path.exists(download_dir):
os.mkdir(download_dir)
 
response = requests.get(feed['rss_url'])
if response.status_code == 200:
now = datetime.now().strftime('%Y%m%d')
name = "%s-%s.rss" % (feed['id'], now)
name = os.path.join(download_dir, name)
with open(name, 'w') as out:
for chunk in response.iter_content(4096):
out.write(chunk)
 
f = feedparser.parse(name)
for entry in f.entries:
yield feed['name'], entry.title
 
 
if __name__ == '__main__':
from collections import defaultdict
 
results = defaultdict(list)
for f,e in ingest():
results[f].append(e)
 
for name, entries in results.items():
print "Downloaded %i Entries from %s" % (len(entries), name)

