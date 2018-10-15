#!/usr/bin/env python

import os
import sys
import subprocess
import shutil
import shlex
import csv
import json
import argparse
import re
from Bio import SeqIO
from Bio import Seq
import warnings
warnings.filterwarnings("ignore")
import MySQLdb
import multiprocessing

database_connection = MySQLdb.connect('', 'ms6', '', '') # edit this line and put your accession information
database_cursor = database_connection.cursor()

#######################################
## arguments ##########################
#######################################

ArgumentParser = argparse.ArgumentParser()
ArgumentParser.add_argument('-m','--mgf',required=True)
ArgumentParser.add_argument('-p','--proteome')
ArgumentParser.add_argument('-g','--genome')
ArgumentParser.add_argument('-r','--replicates',type=str,default=1)
ArgumentParser.add_argument('-o','--outputdir',type=str,default='./')
ArgumentParser.add_argument('-j','--job_id',type=str)
ArgumentParser.add_argument('-ps','--peptideshaker_parameters')
ArgumentParser.add_argument('-c','--contaminants')
Arguments = ArgumentParser.parse_args()

#######################################
## constants ##########################
#######################################

XtandemDir = ""   # edit this line and put your path to the xtandem dir
jobDataDir = ""   # edit this line and put your path to the "/job_data" dir, inside the MS6 dir
crapProtDir = ""  # edit this line and put your path to the "crap_protein" fasta file

CWD = os.getcwd()
BIN_PATHS = {
	'xtandem':'/home/cdtec/Frederico/ms6/bin/xtandem/bin/./tandem.exe',
	'searchGUI':'/home/cdtec/Frederico/ms6/bin/searchGUI/SearchGUI-2.6.4.jar',
	'peptideshaker':'/home/cdtec/Frederico/ms6/bin/peptideshaker/PeptideShaker-1.7.4.jar',
	'interproscan':'/home/cdtec/Frederico/ms6/bin',
	'blast2go':'/home/cdtec/Frederico/ms6/bin/blast2go',
	'swissprot':'/home/cdtec/Frederico/ms6/bin/swissprot',
	'msgfplus':'/home/cdtec/Frederico/re_analise_sibele/programas/MSGFPlus.jar'
}

XtandemOutPath = "{0}/STEP_1".format(Arguments.outputdir)
PeptideShakerOutPath = "{0}/STEP_2".format(Arguments.outputdir)
BLAST2GOOutPath = "{0}/STEP_3".format(Arguments.outputdir)
ProteomeFastaPath = os.path.abspath('{0}/proteome.fasta'.format(Arguments.outputdir))
ProteomeFastaDecoyPath = os.path.abspath('{0}/proteome_concatenated_target_decoy.fasta'.format(Arguments.outputdir))

IdentifiedProteinList = []

#######################################
## functions ##########################
#######################################

def SetStatus(status):
	global Arguments
	global database_connection
	global database_cursor
	database_cursor.execute("UPDATE jobs SET status = '{0}' WHERE job_id = '{1}'".format(
		status,Arguments.job_id))
	database_connection.commit()

def SplitFilePathString(RawFilePathString):
	try:
		SplittedFilePaths = shlex.split(RawFilePathString)
		return SplittedFilePaths
	except:
		sys.exit()

def prepare_dir():
	global Arguments
	global BIN_PATHS
	global CWD
	global XtandemOutPath
	global PeptideShakerOutPath
	global BLAST2GOOutPath
	if not os.path.isdir(Arguments.outputdir):
		os.mkdir(Arguments.outputdir)
	if not os.path.isdir(XtandemOutPath):
		os.mkdir(XtandemOutPath)
	if not os.path.isdir(PeptideShakerOutPath):
		os.mkdir(PeptideShakerOutPath)
		os.mkdir('{0}/MGF'.format(PeptideShakerOutPath))
		os.mkdir('{0}/XML'.format(PeptideShakerOutPath))
		os.mkdir('{0}/REPORTS'.format(PeptideShakerOutPath))
		os.mkdir('{0}/LOGS'.format(PeptideShakerOutPath))
	if not os.path.isdir(BLAST2GOOutPath):
		os.mkdir('{0}'.format(BLAST2GOOutPath))
		os.mkdir('{0}/BLAST_XML'.format(BLAST2GOOutPath))
		os.mkdir('{0}/CDD_BLAST'.format(BLAST2GOOutPath))
		#os.mkdir('{0}/INTERPROSCAN'.format(BLAST2GOOutPath))
		os.mkdir('{0}/ANNOT'.format(BLAST2GOOutPath))
	shutil.copy('/home/cdtec/Frederico/ms6/bin/xtandem/bin/input.xml',
	'{0}/input.xml'.format(XtandemOutPath))
	shutil.copy('/home/cdtec/Frederico/ms6/bin/xtandem/bin/taxonomy.xml',
	'{0}/taxonomy.xml'.format(XtandemOutPath))
	shutil.copy('/home/cdtec/Frederico/ms6/bin/xtandem/bin/default_input.xml',
	'{0}/default_input.xml'.format(XtandemOutPath))
	shutil.copy('/home/cdtec/Frederico/ms6/bin/xtandem/bin/tandem-input-style.xsl',
	'{0}/tandem-input-style.xml'.format(XtandemOutPath))
	shutil.copy('/home/cdtec/Frederico/ms6/bin/xtandem/bin/tandem-style.css',
	'{0}/tandem-style.css'.format(XtandemOutPath))
	shutil.copy('/home/cdtec/Frederico/ms6/bin/xtandem/bin/tandem-style.xsl',
	'{0}/tandem-style.xsl'.format(XtandemOutPath))
	return True

def extract_proteome():
	print 'EXTRAINDO PROTEOMA ...'
	global Arguments
	global BIN_PATHS
	global CWD
	global ProteomeFastaPath
	global ProteomeFastaDecoyPath
	global record
	ProteomeFastaContent = ''
	ProteomeFastaHandle = open(ProteomeFastaPath,'w')
	ProteomeCDSIndex = 0
	for FilePath in SplitFilePathString(Arguments.genome):
		FileHandle = open(FilePath)
		for Scaffold in SeqIO.parse(FileHandle,'genbank'):
			print Scaffold.id
			for Feature in Scaffold.features:
				if Feature.type == 'CDS':
					ProteomeCDSIndex += 1
					CDSSeq = Feature.location.extract(Scaffold)
					CDSProtSeq = Seq.translate(CDSSeq)
					if 'locus_tag' in Feature.qualifiers.keys():
						LocusTag = Feature.qualifiers['locus_tag'][0]
					else:
						LocusTag = 'MISSING_LOCUS_TAG'
					if 'product' in Feature.qualifiers.keys():
						Product = Feature.qualifiers['product'][0]
					else:
						Product = 'MISSING_PRODUCT'
					ProteomeFastaContent += '>{0}|{0} {1} {2} {3}\n{4}\n'.format(
	                						str(ProteomeCDSIndex), Scaffold.id,
	                                        LocusTag,Product.replace("'",""),
	                                        str(CDSProtSeq))

	# add crap contaminant proteins

	crapProteinsHandler = open('/home/cdtec/Frederico/ms6/bin/crap.fasta')
	crapProteinsParser = SeqIO.parse(crapProteinsHandler,'fasta')
	for crapProteinIndex,crapProtein in enumerate(crapProteinsParser):
		ProteomeCDSIndex += 1
		ProteomeFastaContent += '>CONTAMINANT_CRAP_{0}|{0} {1} {2}\n{3}\n'.format(
	                			str(ProteomeCDSIndex), crapProteinIndex,
	                            crapProtein.description.replace("'",""),
	                            str(crapProtein.seq))

	# add custom contaminant proteins

	if record[8]:
		customContaminantProteinsHandler = open('/home/cdtec/Frederico/ms6/jobs_data/%s/contaminants.fasta'%Arguments.job_id)
		customContaminantProteinsParser = SeqIO.parse(customContaminantProteinsHandler,'fasta')
		for customContaminantIndex,customContaminant in enumerate(customContaminantProteinsParser):
			ProteomeCDSIndex += 1
			ProteomeFastaContent += '>CONTAMINANT_CUSTOM_{0}|{0} {1} {2}\n{3}\n'.format(
	                			str(ProteomeCDSIndex), customContaminantIndex,
	                            customContaminant.description.replace("'",""),
	                            str(customContaminant.seq))

	ProteomeFastaHandle.write(ProteomeFastaContent)
	ProteomeFastaHandle.close()
	return True
def prepare_proteome():
	print 'PREPARANDO PROTEOMA'
	global Arguments
	global BIN_PATHS
	global CWD
	global ProteomeFastaPath
	global ProteomeFastaDecoyPath
	global record
	ProteomeFastaContent = ''
	try:
		handle = open(Arguments.proteome)
		parser = SeqIO.parse(handle,'fasta')
		record_index = 0
		for record_index, recordSeq in enumerate(parser):
			ProteomeFastaContent += ">{0}|{1} {2}\n{3}\n".format(
				            'MS6', record_index,
				            recordSeq.description.replace("'",""),
				            str(recordSeq.seq))
		if not record_index:
			sys.exit()

		# add crap contaminant proteins

		crapProteinsHandler = open('/home/cdtec/Frederico/ms6/bin/crap.fasta')
		crapProteinsParser = SeqIO.parse(crapProteinsHandler,'fasta')
		for crapProteinIndex,crapProtein in enumerate(crapProteinsParser):
			record_index += 1
			ProteomeFastaContent += '>CONTAMINANT_CRAP_{0}|{0} {1} {2}\n{3}\n'.format(
	                				str(record_index), crapProteinIndex,
	                            	crapProtein.description.replace("'",""),
	                            	str(crapProtein.seq))

		# add custom contaminant proteins

		if record[8]:
			print record
			customContaminantProteinsHandler = open('/home/cdtec/Frederico/ms6/jobs_data/%s/contaminants.fasta'%Arguments.job_id)
			customContaminantProteinsParser = SeqIO.parse(customContaminantProteinsHandler,'fasta')
			for customContaminantIndex,customContaminant in enumerate(customContaminantProteinsParser):
				record_index += 1
				ProteomeFastaContent += '>CONTAMINANT_CUSTOM_{0}|{0} {1} {2}\n{3}\n'.format(
	                			str(record_index), customContaminantIndex,
	                            customContaminant.description.replace("'",""),
	                            str(customContaminant.seq))
		ProteomeFastaHandle = open(ProteomeFastaPath,'w')
		ProteomeFastaHandle.write(ProteomeFastaContent)
		ProteomeFastaHandle.close()
	except Exception, err:
		print 'ERROR PROTEOMA'
		sys.exit()


def gerar_decoy_database():
	print 'GERANDO BANCO DE DADOS ...'
	global Arguments
	global BIN_PATHS
	global CWD
	global XtandemOutPath
	global PeptideShakerOutPath
	os.chdir(Arguments.outputdir)
	CMD_GENERATEDECOYDATABASE = ('sudo java -cp :{0}: eu.isas.searchgui.cmd.FastaCLI '
		                         '-in {1}/proteome.fasta -decoy').format(
		                         BIN_PATHS['searchGUI'],Arguments.outputdir)
	os.chdir(CWD)
	print CMD_GENERATEDECOYDATABASE
	subprocess.call(CMD_GENERATEDECOYDATABASE,shell=True#,
					#stdout=open(os.devnull,'w'),
					#stderr=open(os.devnull,'w')
        )

def rodar_xtandem():
	print 'RODANDO X!TANDEM ...'
	global Arguments
	global BIN_PATHS
	global CWD
	global XtandemOutPath
	global PeptideShakerOutPath
	MGFCount = 0
	for SpectraFile in SplitFilePathString(Arguments.mgf):
		MGFCount += 1
		shutil.copy(SpectraFile,'{0}/{1}.mgf'.format(XtandemOutPath,str(MGFCount)))
		os.chdir(XtandemOutPath)
		process_inputxml = open('/home/cdtec/Frederico/ms6/bin/xtandem/bin/input.xml'.format(
			               BIN_PATHS['xtandem'])).read().replace("input.mgf","{0}.mgf".format(
						   str(MGFCount)))
		process_outxml = open('input.xml','w')
		process_outxml.write(process_inputxml)
		process_outxml.close()
		process_taxonomyin = open('taxonomy.xml').read().replace("$database$",
		os.path.abspath('../proteome_concatenated_target_decoy.fasta'))
		process_taxonomyout = open('taxonomy.xml','w')
		process_taxonomyout.write(process_taxonomyin)
		process_taxonomyout.close()
		subprocess.call('sudo {0} input.xml'.format(BIN_PATHS['xtandem']),
					shell=True,
					stdout=open(os.devnull,'w'),
					stderr=open(os.devnull,'w'))
		
		shutil.copy('{0}.mgf'.format(str(MGFCount)),'{1}/MGF/{2}.mgf'.format(
			        CWD,PeptideShakerOutPath,str(MGFCount)))
		subprocess.call('java -Xmx5000M -jar {0} -s {2}/MGF/{1}.mgf -d ../proteome_concatenated_target_decoy.fasta -o {2}/XML/{1}.mzid'.format(
			             BIN_PATHS['msgfplus'],str(MGFCount),os.getcwd().replace("STEP_1","STEP_2")),shell=True)
		for FileName in os.listdir('.'):
			if re.match('output\..+\.t.xml',FileName):
				shutil.copy(FileName,'{1}/XML/{2}.t.xml'.format(
				CWD,PeptideShakerOutPath,str(MGFCount)))
				subprocess.call('sudo rm output.*',shell=True,
				stdout=open(os.devnull,'w'),stderr=open(os.devnull,'w'))
				break
		os.chdir(CWD)
	os.chdir(CWD)

def rodar_peptideshaker():
	global Arguments
	global BIN_PATHS
	global CWD
	global XtandemOutPath
	global PeptideShakerOutPath
	os.chdir(PeptideShakerOutPath)
	MGFCount = len(SplitFilePathString(Arguments.mgf))
	CMD_GENERATEREPORT = ('sudo java -cp :{0}: eu.isas.peptideshaker.cmd.ReportCLI '
	                     ' -in {1}/result.cpsx -out_reports {1}/REPORTS '
						 ' -reports 0,1,2,3,4,5,6,7').format(
                         BIN_PATHS['peptideshaker'],os.getcwd())

	CMD_GENERATEPARAMETER = ('sudo java -cp :{0}: eu.isas.peptideshaker.cmd.IdentificationParametersCLI'
	                        ' -out parameters.par -db {1}').format(BIN_PATHS['peptideshaker'],
                             ProteomeFastaDecoyPath)
	CMD_CALLPEPTIDESHAKER = 'sudo java -cp :{0}: eu.isas.peptideshaker.cmd.PeptideShakerCLI'.format(BIN_PATHS['peptideshaker'])

	PeptideShakerCommandline = ('{0} -experiment MS6 -sample MS6_DATA -spectrum_files {1} '
                            '-identification_files {2}  -replicate {3} -out {4}/result.cpsx -id_params '
                            '{4}/parameters.par').format(CMD_CALLPEPTIDESHAKER,os.getcwd()+'/MGF' ,
                            os.getcwd()+'/XML',Arguments.replicates,os.getcwd())
	if not record[7]:
		subprocess.call(CMD_GENERATEPARAMETER,shell=True,stdout=open(os.devnull,'w')
					,stderr=open(os.devnull,'w'))
	else:
		shutil.copy('../parameters.par','%s/parameters.par'%(os.getcwd()))
	subprocess.call(PeptideShakerCommandline,shell=True,stdout=open(os.devnull,'w')
					,stderr=open(os.devnull,'w'))
	subprocess.call(CMD_GENERATEREPORT,shell=True,stdout=open(os.devnull,'w')
					,stderr=open(os.devnull,'w'))
	os.chdir(CWD)

def select_proteins():
	global Arguments
	global BIN_PATHS
	global CWD
	global XtandemOutPath
	global PeptideShakerOutPath
	global BLAST2GOOutPath
	global IdentifiedProteinList
	protein_csv = "{0}/REPORTS/MS6_MS6_DATA_{1}_Default_Protein_Report.txt".format(
		PeptideShakerOutPath,Arguments.replicates)
	for row_i, row in enumerate(csv.reader(open(protein_csv),delimiter='\t')):
		if row_i > 0:
			IdentifiedProteinList.append(row[1].split(' ')[0])
	print IdentifiedProteinList

filteredProteomePATH = '{0}/filtered_proteome.fasta'.format(BLAST2GOOutPath)

def rodar_blast2go():
	global Arguments
	global BIN_PATHS
	global CWD
	global XtandemOutPath
	global PeptideShakerOutPath
	global BLAST2GOOutPath
	global IdentifiedProteinList
	global filteredProteomeHANDLE
	global filteredProteomePATH
	proteomeParser = SeqIO.parse(open("{0}/../proteome_concatenated_target_decoy.fasta".format(PeptideShakerOutPath)),'fasta')
	filteredProteomeFASTA = ''
	for protein in proteomeParser:	
		if (protein.id in IdentifiedProteinList or protein.id.split("|")[0] in IdentifiedProteinList or protein.id.split("|")[1] in IdentifiedProteinList) and ('_REVERSED' not in protein.description):
			filteredProteomeFASTA += protein.format('fasta') + '\n'
	filteredProteomePATH = '{0}/filtered_proteome.fasta'.format(BLAST2GOOutPath)
	filteredProteomeHANDLE = open(filteredProteomePATH,'w')
	filteredProteomeHANDLE.write(filteredProteomeFASTA)
	filteredProteomeHANDLE.close()
	
	# filtered proteome
	
	BLAST_SWISSPROT_COMMANDLINE = 'blastp -db {0}/swissprot -query {1}/filtered_proteome.fasta -evalue 1e-250 -out {1}/BLAST_XML/OUT.xml -num_threads 8 -show_gis -outfmt 5'.format(
		                           BIN_PATHS['swissprot'],BLAST2GOOutPath)

	os.system(BLAST_SWISSPROT_COMMANDLINE)
	os.chdir(BIN_PATHS['blast2go'])
	BLAST2GO_COMMANDLINE = ('java -cp *:ext/*:jdom-1.1.3.jar:mysql-connector-java-5.1.37-bin.jar:biojava.jar:'
	                      ' es.blast2go.prog.B2GAnnotPipe -in {1}/BLAST_XML/OUT.xml -out {1}/ANNOT/OUT -prop b2gPipe.properties -v -annot -dat').format(
	                      BIN_PATHS['blast2go'],BLAST2GOOutPath)
	subprocess.call(BLAST2GO_COMMANDLINE,shell=True,stdout=open(os.devnull,'w'),stderr=open(os.devnull,'w'))
	os.chdir(CWD)
	subprocess.call(('python b2g2json.py {0}/ANNOT/OUT.annot bin/blast2go/go.obo '
					 '{0}/ANNOT/OUT.json').format(BLAST2GOOutPath),stdout=open(os.devnull,'w'),
					 stderr=open(os.devnull,'w'),shell=True)
		
rpsbproc_dict = {}

def rodar_CDDSCcan():
	global Arguments
	global filteredProteomePATH
	global BLAST2GOOutPath
	global rpsbproc_dict
	subprocess.call(('rpsblast+ -query {0} -db bin/cdd-blast/database/Cdd -out '
			  '{1}/CDD_BLAST/cdd.xml -outfmt 5 -num_threads 2').format(
		      filteredProteomePATH,BLAST2GOOutPath),shell=True)
	subprocess.call(('bin/cdd-blast/./rpsbproc -c /home/cdtec/Frederico/ms6/'
					 'bin/cdd-blast/rpsbproc.ini '
			         '-i {0}/CDD_BLAST/cdd.xml -o '
			         '{0}/CDD_BLAST/cdd.rpsbproc').format(BLAST2GOOutPath),
					 shell=True)
	in_domains = False
	for linha in open('{0}/CDD_BLAST/cdd.rpsbproc'.format(BLAST2GOOutPath)):
		if linha[0] == '#':
			continue
		if linha[0:6] == 'QUERY\t':
			protein_id = linha.split('\t')[4].split(" ")[0]
			rpsbproc_dict[protein_id] = []
		if linha == 'DOMAINS\n':
			in_domains = True
		elif linha == 'ENDDOMAINS\n':
			in_domains = False
		if in_domains:
			if linha != 'DOMAINS\n':
				accession = linha.split('\t')[8]
			 	short_name = linha.split('\t')[9]
				domain_start = linha.split('\t')[5]
				domain_end = linha.split('\t')[6]
				identification = linha.split('\t')[3]
				print identification, domain_start,domain_end
				rpsbproc_dict[protein_id].append((accession, short_name, identification))
	#print rpsbproc_dict
	json_rpsproc = open('{0}/CDD_BLAST/cdd.json'.format(BLAST2GOOutPath),'w')
	json_rpsproc.write(json.dumps(rpsbproc_dict))
	json_rpsproc.close()

def gerar_csv_final():
	global Arguments
	global filteredProteomePATH
	global BLAST2GOOutPath
	global PeptideShakerOutPath
	global rpsbproc_dict
	filteredProteomeHANDLE = open(filteredProteomePATH)
	filteredProteomePARSER = SeqIO.parse(filteredProteomeHANDLE,'fasta')
	b2g_json = json.loads(open('{0}/ANNOT/OUT.json'.format(BLAST2GOOutPath)).read())
	protein_csv = csv.reader(open("{0}/REPORTS/MS6_MS6_DATA_{1}_Default_Protein_Report.txt".format(
		PeptideShakerOutPath,Arguments.replicates)))
	print b2g_json
	print protein_csv
	print rpsbproc_dict
	return 0
#######################################
## main ###############################
#######################################

database_cursor.execute('SELECT * FROM jobs WHERE job_id="%s"'%Arguments.job_id)

record = database_cursor.fetchone()

try:
	SetStatus('PROCESSING_DIRS')
	prepare_dir()
except:
	SetStatus('ERROR_DIRS')
	sys.exit()	
try:
	SetStatus('PROCESSING_REFERENCE')
	if Arguments.proteome:
		prepare_proteome()
	elif Arguments.genome:
		extract_proteome()
except:
	SetStatus('ERROR_GENOME')
	sys.exit()

try:
	SetStatus('PROCESSING_DATABASE')
	gerar_decoy_database()
except:
	SetStatus('ERROR_DECOY')
	sys.exit()	

try:
	SetStatus('PROCESSING_IDENTIFICATION')
	rodar_xtandem()
except:
	SetStatus('ERROR_XTANDEM')
	sys.exit()
try:
	rodar_peptideshaker()
except:
	SetStatus('ERROR_PEPTIDESHAKER')
	sys.exit()
try:
	select_proteins()
except:
	SetStatus('ERROR_PROTEIN_IDS')
	sys.exit()

try:
	SetStatus('PROCESSING_FUNCTIONAL_ANNOTATION')
	rodar_blast2go()
except:
	SetStatus('ERROR_BLAST2GO')
	sys.exit()

try:
	SetStatus('PROCESSING_CDD_DOMAIN_ANNOTATION')
	rodar_CDDSCcan()
except:
	SetStatus('ERROR_CDD')
	sys.exit()

try:
	SetStatus('FINAL_PROCESSING')
	gerar_csv_final()
except:
	SetStatus('ERROR_FINAL_PROCESSING')
	sys.exit()
os.system('rm -r {0}/spectra'.format(Arguments.outputdir))
#os.system('rm -r {0}/proteome*'.format(Arguments.outputdir))
#os.system('rm {0}/STEP_1/*.mgf'.format(Arguments.outputdir))
os.system('rm {0}/*raw_ref'.format(Arguments.outputdir))
os.chdir('jobs_data')
os.system('tar -zcvf {0}.tgz {0}'.format(Arguments.outputdir))
os.chdir('..')
os.system("sudo chmod -R 777 jobs_data/*")
SetStatus('FINISHED')

#except:
#	pass
