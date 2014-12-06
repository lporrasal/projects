import sys

filename = sys.argv[1]

lines = open(filename).readlines()

print '<?xml version="1.0"?>'
print '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1 plus MathML 2.0 plus SVG 1.1//EN"'
print '              "http://www.w3.org/Math/DTD/mathml2/xhtml-math11-f.dtd">'
print '<?xml-stylesheet type="text/xsl" href="mathml.xsl"?>'
print '<html xmlns="http://www.w3.org/1999/xhtml">'


for line in lines[1:]:
    line = line.strip()
    if line.startswith("<meta"): line = line + "</meta>"
    if line.startswith("<link"): line = line + "</link>"
    if line.startswith("<li>"): line = line + "</li>"
    if line.startswith("<img"):
       i = line.index(">")
       line = line[:i] + " />" + line[i+1:]
    line = line.replace("</li></ul>","</ul>")
    line = line.replace("<br></td>","</td>")
    line = line.replace("&ldquo;","&#8220;")
    line = line.replace("&rdquo;","&#8221;")
    line = line.replace("&iuml;","&#239;")
    line = line.replace("<br>","<br />")
    line = line.replace("<dd>","</dt><dd>")
    if "<p>" in line:
        if not "</p>" in line:
            line = line.replace("<p>","")
    if "<img" in line:
        i = line.index("<img")
        i = line.index(">", i)
        if line[i-1]!="/":
           line = line[:i] + "/" + line[i:]
    line = line.replace('<p class="noindent">',"")
    line = line.replace("index.html","index.xhtml")
    if "<dd>" in line: line = line + "</dd>"
    print line
