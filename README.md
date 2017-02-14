# ms6

## About

MS6 provides an automated pipeline to identify proteins from MS/MS data by integrating 
the programs SearchCLI (Vaudel et al., 2011), X!Tandem (Fenyö and Beavis, 2003) and PeptideShaker 
(Vaudel et al., 2015). Additionally, our pipeline also generates an structural and 
functional report of the identified proteins using the CDD database 
(Marchler-Bauler et al., 2015) and BLAST2GO (Conesa et al., 2005), 
as well provides all intermediate files. MS6 only requires the spectra files 
(“.mgf”) and a reference genome in Genbank format (“.gb”) or a proteome FASTA (“.fasta”) 
to be executed. The front-end of MS6 was built on top of the Bootstrap and DataTables frameworks.
