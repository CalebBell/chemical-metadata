from codecs import open
import glob
import json
import os

try:
    from HTMLParser import HTMLParser
except ImportError:
    from html.parser import HTMLParser
    
'''Runs on unix. Requires htmltopdf.
Not fully python 3 compliant.

'''
def checkCAS(CASRN):
    try:
        check = CASRN[-1]
        CASRN = CASRN[::-1][1:]
        productsum = 0
        i = 1
        for num in CASRN:
            if num == '-':
                pass
            else:
                productsum += i*int(num)
                i += 1
        return (productsum % 10 == int(check))
    except:
        return False

h = HTMLParser()
folder = os.path.join(os.path.dirname(__file__))

def parse_scifinder(data):
    d = {}
    possibilities = data.split('CAS Registry Number')[1].split('<b>')[1].split('</b>')[1].split('<br/>')#[0]
    if len(possibilities) > 2 and 'Manual Registration' in possibilities[2]:
        name = possibilities[1]
    else:
        name = possibilities[0]
    
    d['Name'] = name
    if 'Alternate CAS Registry Numbers' in data:
        alt = data.split('Alternate CAS Registry Numbers')[1].split('<b>')[0].split('</b>')[1]
        alt = alt.split('<')[0] # remove the end tag
        alt = alt.split('; ') # If there are multiple hits - or split on ',' and include br
        d['Alternate CAS'] = alt
    if 'Deleted CAS Registry Numbers' in data:
        deleted = data.split('Deleted CAS Registry Numbers')[1].split('\n')[0].split('</b>')[1].replace('<br/>', '').replace(' ', '').split(',')
        d['Deleted CAS'] = deleted
    if 'Other Names' in data:
        hits = data.split('Other Names')[1]#.split('\n')[0].split('</b>')[1].replace('<br/>', '')
        if 'Deleted CAS Registry Numbers' in hits:
            hits = hits.split('Deleted CAS Registry Numbers')[0]
        if 'Alternate CAS Registry Numbers' in hits:
            hits = hits.split('Alternate CAS Registry Numbers')[0]
        else:
            hits = hits.split('References')[0]
        hits = hits.split('\n')[0].split('</b>')[1].replace('<br/>', '')
        
        for bad_text in ['Alternate CAS Registry Numbers', 'Incompletely Defined Substance', 'Coordination Compound', 'Manual Registration', 'Tabular Inorganic Substance']:
            hits = hits.split(bad_text)[0]
        hits = hits.replace('<b>', '')
        names = [i.strip() for i in hits.split(';')]
        d['Other Names'] = names
    return d



if len(os.listdir('pdf')) != len(os.listdir('html')):
    os.system('cd %s'%os.path.join(folder, 'pdf') + '; for i in *.pdf ; do pdftohtml -i -enc UTF-8 -noframes  "$i" "../html/${i%.*}.html" ; done' )


parsed_data = {}
for html in glob.glob(os.path.join(os.path.join(folder,'html'), '*.html')):
    CAS = html.split('.html')[0].split('/')[1]
    text = h.unescape(open(os.path.join(folder, html)).read().decode('utf-8'))
    dat = parse_scifinder(text)
    parsed_data[CAS] = dat
    
#print(parsed_data)

#print(len([i for i in parsed_data.values() if ('Alternate CAS' in i and i['Alternate CAS'])]))
#assert 13 == len([i for i in parsed_data.values() if ('Alternate CAS' in i and i['Alternate CAS'])])
#print(len([i for i in parsed_data.values() if ('Deleted CAS' in i and i['Deleted CAS'])]))
#assert 226 == len([i for i in parsed_data.values() if ('Deleted CAS' in i and i['Deleted CAS'])])
#for CAS, i in parsed_data.items():
#    if 'Other Names' in i:
#        print(i['Other Names']) 

f = open('Parsed scifinder metadata2.json', 'w')
#print(parsed_data)
json.dump(parsed_data, f, indent=2, separators=(',', ': '), encoding="utf-8", sort_keys=True)
f.close()
