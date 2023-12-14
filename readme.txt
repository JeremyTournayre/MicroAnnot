>Requirements
-tRNAscan-SE must be installed and available in the path
-perl
-Bioperl with Bio::Seq, Bio::SeqIO and Bio::SearchIO;
-perl modules: Storable qw(dclone) and Clone qw(clone);
-ncbi-blast-2.13.0+
-glimmer3.02
-Bioperl modification: add sort on row 1264 on genbank.pm (usr/share/perl5/Bio/SeqIO)
    # foreach my $tag (sort keys %{ $fth->field } ) { 
	# instead of     
	# foreach my $tag (keys %{ $fth->field } ) {

# add sort on row 956 on embl.pm  (usr/share/perl5/Bio/SeqIO)
# foreach my $tag (sort keys %{$fth->field} ) {
# instead of     
# foreach my $tag (keys %{ $fth->field } ) {
	
# 'eq' instead of '==' on row 350 on Simple.pm (usr/share/perl5/Bio/Location)


>Use of the microannot.pl script:
perl microannot.pl input_file learn_Glim_orf_file_icm min_orf_size bool_interpro evalue_alignment evalue_small evalue_TE glimmer_size debug list_dat_name list_dat_id dir_blast dir_db dir_glimmer

input_file: Input sequence(s) (fasta file)
learn_Glim_orf_file_icm: Training dataset (used only if less than 50 CDS are identified by homology)
min_orf_size: Minimum size for ORF finding (min: 200 nt)
bool_interpro: Use of InterproScan: 1 for true or 0 for false
evalue_alignment: Evalue for the blastp alignment for CDS prediction by comparative approach (CDS >= 80 AA)
evalue_small: Evalue for the blastx alignment for small CDS prediction by comparative approach (CDS < 80 AA)
evalue_TE: Evalue for the tblastx alignment for transposable elements detection
glimmer_size: Minimum size for CDS prediction (min: 300)
debug: Use of the debug printing: 1 for true or 0 for false
list_dat_name: List of databases to perform blastp alignment for the CDS prediction by comparative approach. Choose all databases 'A algerae,E bieneusi,E cuniculi,N ceranae' or only some databases 'E bieneusi,E cuniculi'
list_dat_id: Active or not a print of the list_dat_name in the log file: 1 for true or 0 for false
dir_blast: Path of the ncbi-blast-2.13.0+/bin directory
dir_db: Path of the databases directory
dir_glimmer: Path of the /glimmer3.02/bin directory


>Example of an execution of this script with the sample: input/sample.fsa_nt
perl microannot.pl input/sample.fsa_nt A_algerae_complete_CDS.icm 240 1 -15 -5 -10 300 1 'A algerae,E bieneusi,E cuniculi,N ceranae' 0 ./tools/ncbi-blast-2.13.0+/bin ./db ./tools/glimmer3.02/bin

>Contents of the results directory:

**input_file**: the input file

#Alignment/small CDS < 80 aa
small.bls: blastx results against the small CDS database
#tRNA
trna.txt: tRNAscan-SE results
#rRNA
rRNA.bls: blastn results against the 16_rRNA database

#Alignment/CDS >= 80 aa
orf.fa: List of ORFs found in the input file
orf.bls: blastp results against the orf.fa file
#glimmer
learn_Glim_orf.fa: List of CDS used for train glimmer 
res_glimmer.predict: glimmer Results 
res_glimmer.detail: glimmer Results

#Transposable element
CDS_gene_nt.fa: List of CDS used for the blast of transposable element
TE.bls: tblastx results against the transposable element database

#InterproScan
CDS_aa.fa: List of CDS used for InterproScan
Res_Interproscan: Results of InterproScan

#Intermediate results
Res_gb: Intermediate MicroAnnot results in gb format without anotations of interproscan and transposable elements  
Res_embl: Intermediate results of MicroAnnot in embl format without interproscan and transposable elements anotations

#MicroAnnot results:
Res_gff_annot: in gff format
Res_gb_annot: in gb format
Res_embl_annot: in embl format
Res_warnings: Warnings only

#MicroAnnot results in tar.gz :
Res_gb_annot.tar.gz
Res_embl_annot.tar.gz
Res_gff_annot.tar.gz
Res_warnings.tar.gz