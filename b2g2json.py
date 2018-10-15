#!/usr/bin/env python

import csv
import sys
import re
import json

b2g_file = sys.argv[1]
obo_file = sys.argv[2]
jsn_file = sys.argv[3]

'''
b2g_file = 'jobs_data/KTWVSQKAHA1D7WLJV31K/STEP_3/ANNOT/OUT.annot'
obo_file = 'bin/blast2go/go.obo'
jsn_file = sys.argv[1]
'''

def get_def(id):
	if 'GO:' in id:
		global obo_file
		current_name = ''
		current_idlist = []
		for linha in open(obo_file):
			if linha == '[Term]\n':
				if id in current_idlist:
					return current_name
				current_idlist = []
				current_name = ''
			if re.match("id: .+",linha):
				current_idlist.append(linha.split(': ',1)[1].replace('\n',''))
			if re.match("name: .+",linha):
				current_name = linha.split('name: ',1)[1].replace('\n','')
			if re.match("alt_id: .+",linha):
				current_idlist.append(linha.split(': ',1)[1].replace('\n',''))
		else:
			return current_name
	else:
		return 'Enzyme Code'
proteins = []

protein_data = {}

for linha in csv.reader(open(b2g_file),delimiter='\t'):
	if 'id' in protein_data:
		if linha[0] == protein_data['id']:
			pass
		else:
			proteins.append(protein_data)
			protein_data = {}
			protein_data['id'] = linha[0]
			protein_data['GO'] = []
	else:
		protein_data['id'] = linha[0]
		protein_data['GO'] = []
	print linha
	protein_data['GO'].append(
                                    (linha[1],get_def(linha[1]))
                                 )
proteins.append(protein_data)

resultado = open(jsn_file,'w')
resultado.write(json.dumps({'DATA':b2g_file,'PROT':proteins}))
print jsn_file
resultado.close()
