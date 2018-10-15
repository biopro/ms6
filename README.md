# ms6

## About

MS6 provides an automated pipeline to identify proteins from MS/MS data by integrating 
the programs SearchCLI ([Vaudel et al., 2011](http://www.ncbi.nlm.nih.gov/pubmed/12622365))
, X!Tandem ([Fenyö and Beavis, 2003](http://www.ncbi.nlm.nih.gov/pubmed/12622365)) 
and PeptideShaker ([Vaudel et al., 2015](http://www.ncbi.nlm.nih.gov/pubmed/25574629)). 
Additionally, our pipeline also generates an structural and 
functional report of the identified proteins using the CDD database 
([Marchler-Bauler et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383992/)) 
and BLAST2GO ([Conesa et al., 2005](http://www.ncbi.nlm.nih.gov/pubmed/16081474)), 
as well provides all intermediate files. MS6 only requires the spectra files (“.mgf”) 
and a reference genome in Genbank format (“.gb”) or a proteome FASTA (“.fasta”) 
to be executed. The front-end of MS6 was built on top of the 
[Bootstrap](http://getbootstrap.com/) and [DataTables](https://datatables.net/) 
frameworks.

The MS6 webserver is available at http://labbioinfo.ufpel.edu.br/ms6/eng/.

## Deploy

1. Install the `BioPython` package.
1. Run the `ms6_database.sql` script on your system.
2. Edit the `ms6_analysis.py` script and change the paths to programs and directories.
3. Run the program.

## Developers
MSc. Frederico Schmitt Kremer  
*Master of Science (Biotechnology).*
**contact:** fred.s.kremer@gmail.com.

Martina Bianca Fuhrmann 
*Master degree student (Genetics).*
**contact:** martinabfuhrmann@gmail.com.

Dr. Luciano da Silva Pinto 
*Professor at Universidade Federal de Pelotas.*
**contact:** ls_pinto@hotmail.com.

### BioPro - Laboratório de Bioinformática e Proteômica

MS6 was developed by the BioPro Laboratory team at Universidade Federal Federal de Pelotas. Our group develops researches on many areas of biotechnology, bioinformatics, genomics and proteomics, including: production of recombinant proteins, isolation and caracterization of plant lectins, microbial genomics and software development for bioinformatics.

### Universidade Federal de Pelotas
The Federal University of Pelotas (Portuguese: Universidade Federal de Pelotas, UFPEL) is a higher-learning facility with campuses in Pelotas and Capão do Leão in Rio Grande do Sul, Brazil. The university offers more than a hundred different graduation courses.

## Citing MS6

If you have used MS6 in your research please cite:

Kremer, F.S. et al.(2017)

### The methodology which have inspired the development of MS6:

Esperotto, R.L. _et al_.(2017) Proteomic analysis of _Toxocara canis_ excretory and secretory (TES) proteins. **Molecular and Biochemical Parasitology**, v. 211, p39-47. DOI: [10.1016/j.molbiopara.2016.09.002](http://dx.doi.org/10.1016/j.molbiopara.2016.09.002). 
PubMed: [27638150](https://www.ncbi.nlm.nih.gov/pubmed/27638150). 
