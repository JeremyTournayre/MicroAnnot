use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Storable qw(dclone);
use Clone qw(clone);


my $input_file=$ARGV[0];
my $learn_Glim_orf_file_icm=$ARGV[1];
my $min_orf_size=$ARGV[2];
my $bool_interpro=$ARGV[3];
my $evalue_alignment="1e".$ARGV[4];
my $evalue_small="1e".$ARGV[5];
my $evalue_TE="1e".$ARGV[6];
my $glimmer_size=$ARGV[7];
my $debug=$ARGV[8];
my $list_dat_name=$ARGV[9];
my $list_dat_id=$ARGV[10];

my $dir_blast=$ARGV[11]; #with the /bin
my $dir_db=$ARGV[12];
my $dir_glimmer=$ARGV[13]; #with the /bin

if ($min_orf_size eq "" ){
	$min_orf_size=240;
}
if ($min_orf_size <= 200){
	$min_orf_size=200;
}
if ($glimmer_size eq "" ){
	$glimmer_size=300;
}
if ($glimmer_size <= 300){
	$glimmer_size=300;
}
if ($evalue_alignment eq "" ){
	$evalue_alignment=1e-15;
}
if ($evalue_small eq "" ){
	$evalue_small=1e-5;
}
if ($evalue_TE eq "" ){
	$evalue_TE=1e-10;
}

$learn_Glim_orf_file_icm=$dir_db."/db_ref_Glimmer/".$learn_Glim_orf_file_icm;

if (!-e $learn_Glim_orf_file_icm){
	$learn_Glim_orf_file_icm=~m/(.*)\.icm$/;
	my $learn_Glim_orf_file=$1;
	#Please ask to jeremy.tournayre@inrae.fr to have the databases
	`$dir_glimmer/build-icm -r $learn_Glim_orf_file_icm < $learn_Glim_orf_file`;
}

my $mess_interpro="Activate";
if ($bool_interpro==0){
	$mess_interpro="Disabled";
}

print "Options:\n";
print 'Sequence(s) file:'.'"'.$input_file.'"'."\n";
print 'Minimum ORF size (min: 200):'.'"'.$min_orf_size.'"'."\n";
print 'Glimmer size (min: 300):'.'"'.$glimmer_size .'"'."\n";
print 'Training dataset for Glimmer (used if there is less than 50 CDS):'.'"'.$learn_Glim_orf_file_icm.'"'."\n";
print 'Interproscan:'.'"'.$mess_interpro.'"'."\n";
print 'Evalue for the alignment blast:'.'"'.$evalue_alignment.'"'."\n";
print 'Evalue for the small cds blast:'.'"'.$evalue_small.'"'."\n";
print 'Evalue for the transposable element blast:'.'"'.$evalue_TE.'"'."\n";
if ($debug==1){print 'Debug:'.'"'.$debug.'"'."\n";}
if ($list_dat_id ne 0){print 'list_dat_name:'.'"'.$list_dat_name.'"'."\n";}

print "\n";

$input_file=~/.*\/(.*)$/;
my $file_name=$1;
my $new_dir_o="output/".$file_name;
`mkdir -p $new_dir_o`;
my $new_input_file=$new_dir_o."/".$file_name;

open (F,$input_file);
open (O,">".$new_input_file);
my %id_2_name;
my %number_name;
my $id=0;
while(<F>){
	chomp($_);
	$_=~s/ //g;
	$_=~s/\t//g;
	$_=~s/\r//g;
	$_=~s/\n//g;
	$_=uc($_);
	if ($_=~m/>/){
		$id++;
		$id_2_name{$id}=$_;		
		print O '>'.$id."\n";
	}elsif($id!=0){
		print O $_."\n";
	}
}
close(F);
close(O);

`dos2unix $new_input_file`;

my $in = Bio::SeqIO ->new (-file => $new_input_file,-format => 'Fasta');
my $orf_file=$new_dir_o."/".'orf.fa';

############ORF
open (O_ORF,'>'.$orf_file);
my %annotation;

while (my $seq=$in->next_seq()) {
	$annotation{$seq->display_id}{"bioseq"}=$seq;
	#6 frames
	my $revcom=$seq->revcom;
	for(my $i_cadre=0;$i_cadre<6;$i_cadre++){
		my $n_cadre=$i_cadre;
		my $i_seq=$seq;
		my $cadre="+".($i_cadre+1);
		if ($i_cadre>2){
			$n_cadre=$i_cadre-3;
			$i_seq=$revcom;
			$cadre="-".($n_cadre+1);
		}
								
		my $translate=$i_seq->translate("","",$n_cadre)->seq;
		my @split_orf=split('\*',$translate);
		my $bool_extrem_5prim=1;
		my $i=0;
		my $start=1+$n_cadre;
		foreach (@split_orf){
			my $end=$start+(length($_)*3)-1;
			if (($end-$start+1)>=$min_orf_size){
				my $extrem="";
				if ($bool_extrem_5prim==1){
					$extrem="_extr_5prim";
				}
				elsif ($i == $#split_orf){
					$extrem="_extr_3prim";
				}
				print O_ORF ">".$seq->display_id."_orf_".($cadre)."_".($start)."-".$end.$extrem."\n".$_."\n";
				$annotation{$seq->display_id}{$seq->display_id."_orf_".($cadre)."_".($start)."-".$end.$extrem}{"seqAA"}=$_;
				
				my $start_nt=$start-20;
				my $add_nt=20;
				if ($start_nt <=0){
					$add_nt=$start;
					$start_nt=0;
				}
				$annotation{$seq->display_id}{$seq->display_id."_orf_".($cadre)."_".($start)."-".$end.$extrem}{"seqNT-add-n"}=$add_nt;
				$annotation{$seq->display_id}{$seq->display_id."_orf_".($cadre)."_".($start)."-".$end.$extrem}{"seqNT-add"}=substr($i_seq->seq,$start_nt,$end);
			}
			$i++;
			$bool_extrem_5prim=0;
			$start+=(length($_)*3)+3;
		}	
	}
}
close(O_ORF);

#####
my $result;
my $display_id;
my $start;
my $amont_M_20nt;
my $note_signal;
my $bool_signal;

sub signal {
	my $tmp_amont_M_20nt=$amont_M_20nt;
	my $AT= $tmp_amont_M_20nt=~s/A|T/X/g;
	my $perc_AT=$AT/length($amont_M_20nt);
	$note_signal="";
	$bool_signal=0;
	my $bool_strong=0;
	if ($perc_AT >=0.8){
		$note_signal.="(percent of AT >= 80% into the 20nt before the M)";
		$bool_signal=1;
		$bool_strong=1;
	}
	
	if ($amont_M_20nt =~m/CCC/){
		$note_signal.="(CCC into the 20nt before the M)";
		$bool_signal=1;
		$bool_strong=1;
	}
	if ($amont_M_20nt =~m/GGG/){
		$note_signal.="(GGG into the 20nt before the M)";
		$bool_signal=1;
		$bool_strong=1;
	}	
	if ($bool_strong==1){
		$note_signal="-> strong signal ".$note_signal."";
	}
	if ($bool_strong==0){
		if ($amont_M_20nt =~m/CC/){
			$note_signal.="(CC into the 20nt before the M)";
			$bool_signal=1;
		}
		if ($amont_M_20nt =~m/GG/){
			$note_signal.="(GG into the 20nt before the M)";
			$bool_signal=1;
		}
		if ($bool_signal==1){
			$note_signal="-> weak signal ".$note_signal."";
		}
	}
}


#####Function annotation alignment
my $hit;
my $hsp;
my %hit_2_frameshift;
my %small_cds;

sub annot_alignment {
	my $bool_hit=$_[0];
	$amont_M_20nt="";
	my $amont10="";
	my @split_orf;
	my @split2_orf;
	if ($debug){print $result->query_name."\n";}
	@split_orf=split("_",$result->query_name);
	@split2_orf=split("-",$split_orf[3]);
	$display_id=$split_orf[0];
	if ($split_orf[2] =~m/^-/){
		$hit_2_frameshift{$display_id}{"-"}{$split2_orf[1]}{"orf"}=$result->query_name;
		$hit_2_frameshift{$display_id}{"-"}{$split2_orf[1]}{"hit"}{$hit->name()}=1;
	}else{
		$hit_2_frameshift{$display_id}{"+"}{$split2_orf[0]}{"orf"}=$result->query_name;
		$hit_2_frameshift{$display_id}{"+"}{$split2_orf[0]}{"hit"}{$hit->name()}=1;
	}
	
	if ($bool_hit==1){
		return;
	}
	my $anno_start="";
	my $anno_end="";
	my $note="";
	my $warning=0;
	if ($hsp->hit_string =~m/^M/ && $hsp->query_string =~m/^M/ && $hsp->start("hit")== 1){
		$anno_start=$hsp->start("query");
		$anno_start--;	
		$note="Alignment starting with M";
		$warning=0;

	}elsif($result->query_name=~m/extr/ && 
		(($split_orf[2] =~m/^-/ && $result->query_name =~m/_5prim/) || 
		($split_orf[2] !~m/^-/ && $result->query_name =~m/_5prim/))){
		$note="Start undetermined";
		$anno_start=0;
	}else{
		my $substr_a;
		my $substr_b=10;
		
		$substr_a=(($hsp->start("query")-$hsp->start("hit"))-5);
		if ($substr_a<0){
			$substr_b=$substr_b+$substr_a;
			$substr_a=0;
		}		
		if (($substr_a+$substr_b)>=($hsp->start("query")+5)){
			$substr_b=($hsp->start("query")-1);
		}	
		my $amont10="";
		if ($substr_b >0){
			$amont10=substr($annotation{$display_id}{$result->query_name}{"seqAA"},$substr_a,$substr_b);
		}

		if ($amont10=~m/(.*)M/){
			$anno_start=$substr_a+(length($1));	
			$note="Alignment with a M into the 5AA before or after the start of the hit";
			$warning=0;				
		}else{
			my $start=-3;
			my $amont;
			$amont=substr($annotation{$display_id}{$result->query_name}{"seqAA"},0,($hsp->start("query")-1));
			my $bool=0;
			my $bool_get_firstM="";
			$bool_signal=0;
			while ($amont=~m/(.*?)M/g){
				my $amont_M=$1;
				$start+=3;
				my $amont_M_nt;
				if ($split_orf[2] =~m/^-/){
					$start+=(length($amont_M)*3);	
				}else{
					$start+=(length($amont_M)*3);
				}
				$amont_M_nt=substr($annotation{$display_id}{$result->query_name}{"seqNT-add"},0,(($start-1)+$annotation{$display_id}{$result->query_name}{"seqNT-add-n"}));
				$amont_M_20nt=substr($amont_M_nt,-20);
				signal();
				if ($bool_signal==1){
					$anno_start=$start/3;
					$note="Alignment without M but in upstream there is a M with a signal ".$note_signal;
					$warning=0;				
					last;
				}
				if ($bool==0){
					$bool_get_firstM=$start/3;
				}
				$bool=1;
			}
			if ($bool_signal==0 && $bool==1){
				$anno_start=$bool_get_firstM;
				$note="Alignment without M but in upstream there is a M without a signal.";
				$warning=1;
			}elsif ($bool_signal==0){
				$anno_start=$hsp->start("query");
				$anno_start--;
				$note="Alignment without M and there is no M at all upstream of the alignment.";
				$warning=1;
			}
		}
	}
	$anno_start=$anno_start*3;
	if ($split_orf[2] =~m/^-/){
		$anno_end=length($annotation{$display_id}{"bioseq"}->seq)-$split2_orf[0]-$anno_start+1;
		$anno_start=length($annotation{$display_id}{"bioseq"}->seq)-$split2_orf[1]+1;
	}else{
		$anno_start+=$split2_orf[0];
		$anno_end=$split2_orf[1];
	}
	
	$annotation{$display_id}{$result->query_name}{"CDS"}{"start"}=$anno_start;
	$annotation{$display_id}{$result->query_name}{"CDS"}{"end"}=$anno_end;
	$annotation{$display_id}{$result->query_name}{"CDS"}{"note"}=$note;
	if ($debug==1){
		if ($amont10 ne ""){
			$annotation{$display_id}{$result->query_name}{"CDS"}{"note"}.="    The 5AA +/- start(query - hit):".$amont10.".";
		}			
		if ($amont_M_20nt ne ""){
			$annotation{$display_id}{$result->query_name}{"CDS"}{"note"}.="    The 20nt upstream:".$amont_M_20nt.".";
		}
	}
	$annotation{$display_id}{$result->query_name}{"CDS"}{"warning"}=$warning;
	$annotation{$display_id}{$result->query_name}{"CDS"}{"hit"}=$hit->name();
	$annotation{$display_id}{$result->query_name}{"CDS"}{"evalue"}=$hsp->evalue();		
}

my $orf_stop="";
my $orf_start="";
my $anno_start="";
my $anno_end="";
sub find_orf{
	my $strand=$_[0];
	my $frame=$_[1];
	my $start_query=$_[2];
	my $end_query=$_[3];
	if ($strand ==-1){
		my $a=(($annotation{$display_id}{"bioseq"}->length-$frame)%3)+1;
		my $b=$start_query-1;
		if ($debug==1){print 'translate_a'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
		if ($start_query > 3 && $annotation{$display_id}{"bioseq"}->trunc($a,$b)->revcom->translate->seq =~ m/(.*?)\*/){
			$anno_start=$start_query -(length($1)*3+3);	
		}else{
			$anno_start=(($annotation{$display_id}{"bioseq"}->length-($frame))%3)+1; 
		}
		$orf_stop=$anno_start;
		$a=$end_query+1;
		$b=$annotation{$display_id}{"bioseq"}->length-$frame;
		if ($debug==1){print 'translate_b'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
		if ($a < $b && $annotation{$display_id}{"bioseq"}->trunc($a,$b)->revcom->translate->seq =~ m/.*\*(.*)/){
			$orf_start=(length($1)*3)+$end_query;  
		}else{
			$orf_start=$annotation{$display_id}{"bioseq"}->length-($frame);
		}
	}else{
		my $a=$end_query+1;
		my $b=($annotation{$display_id}{"bioseq"}->length-($annotation{$display_id}{"bioseq"}->length-$frame)%3);
		if ($debug==1){print 'translate_c'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
		if ($end_query+1 < ($annotation{$display_id}{"bioseq"}->length-3) && $annotation{$display_id}{"bioseq"}->trunc($a,$b)->translate->seq =~ m/(.*?)\*/){
			
			$anno_end=$end_query+(length($1)*3+3);	
		}else{		
			$anno_end=$annotation{$display_id}{"bioseq"}->length-(($annotation{$display_id}{"bioseq"}->length-($frame))%3); 
		}
		$orf_stop=$anno_end;
		$a=$frame+1;
		$b=$start_query-1;
		if ($debug==1){print 'translate_d'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
		if ($frame+1 < $start_query-1 && $annotation{$display_id}{"bioseq"}->trunc($a,$b)->translate->seq =~ m/(.*)\*/){
			$orf_start=(length($1)*3+3)+$frame+1;
		}else{
			$orf_start=$frame+1;
		}
	}
}


my %orf_hit_small;
sub annot_alignment_small {
	$display_id=$result->query_name();
	$anno_start="";
	$anno_end="";
	my $note="";
	my $warning=0;
	my $pos_small=0;
	my $bool_small=0;
	my $amont10="";
	$amont_M_20nt="";
	$orf_stop="";
	$orf_start="";
	find_orf($hsp->strand,$hsp->frame("query"),$hsp->start("query"),$hsp->end("query"));
	if (!defined($orf_hit_small{$orf_stop."-".$hsp->strand})){
		if ($hsp->hit_string =~m/^M/ && $hsp->query_string =~m/^M/ && $hsp->start("hit")== 1){
			if ( ($hsp->strand == 1)){		
				$anno_start=$hsp->start("query");
			}else{
				$anno_end=$hsp->end("query");
			}	
			$bool_small=1;
			$note="Alignment starting with M.";
			$warning=0;
		}else{
			my $substr_a;
			my $substr_b=10*3;
			if ($hsp->strand == 1){
				$substr_a=((($hsp->start("query"))-($hsp->start("hit")*3))+3-(5*3));
				if ($substr_a<=0){
					$substr_b=$substr_b+$substr_a-3-$hsp->frame("query")-1;
					$substr_a=$hsp->frame("query")+1;
				}		
				if ($substr_b>0){
					my $a=$substr_a;
					my $b=($substr_a+$substr_b)-1;
					if ($debug==1){print 'translate_e'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
					$amont10=$annotation{$display_id}{"bioseq"}->trunc($a,$b)->translate->seq;
				}
			}else{
				$substr_a=((($hsp->end("query"))+($hsp->start("hit")*3))-3+(5*3));
				if ($substr_a>length($annotation{$display_id}{"bioseq"}->seq)){
					$substr_b=$substr_b-($substr_a+3-$annotation{$display_id}{"bioseq"}->length+$hsp->frame("query"));
					$substr_a=$annotation{$display_id}{"bioseq"}->length-$hsp->frame("query");
				}		
				if ($substr_b>0){	
					my $a=$substr_a-$substr_b+1;
					my $b=($substr_a);
					if ($debug==1){print 'translate_f'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
					$amont10=$annotation{$display_id}{"bioseq"}->trunc($a,$b)->revcom->translate->seq;	
				}
			}
			if ($amont10 ne "" && $amont10=~m/(.*\*)(.*)/){
				$amont10=$2;
				if ($hsp->strand == 1){
					$substr_a=$substr_a+(length($1)*3);
				}else{
					$substr_a=$substr_a-(length($1)*3);
				}
			}
			if ($amont10 ne "" && (($hsp->strand == 1 && $orf_start <= $substr_a)	|| ($hsp->strand == -1 && $orf_start >= $substr_a) )
					&& $amont10=~m/(.*)M/){
				
				if ($hsp->strand == 1){
					$anno_start=$substr_a+(length($1)*3);	
				}else{
					$anno_end=$substr_a-(length($1)*3);	
				}
				$note="Alignment with a M into the 5AA before or after the start of the hit.";
				$warning=0;				
			}else{
				my $start=0;
				my $amont;
				if ( $hsp->strand == 1){	
					my $a=$hsp->frame("query")+1;
					my $b=$hsp->start("query")-1;
					if ($debug==1){print 'translate_g'.((($b)-($a)+1)%3)."\n";	if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
					$amont="";
					if ($a < $b){					
						$amont=$annotation{$display_id}{"bioseq"}->trunc($a,$b)->translate->seq;
					}
					if ($amont=~m/(.*\*)(.*)/){
						$amont=$2;				
						$pos_small=(length($1)*3)+$hsp->frame("query")+1;
					}else{
						$pos_small=$hsp->frame("query")+1;
					}	
				}else{
					my $a=$hsp->end("query")+1;
					my $b=$annotation{$display_id}{"bioseq"}->length-$hsp->frame("query"); 
					if ($debug==1){print 'translate_h'.((($b)-($a)+1)%3)."\n";	if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
					$amont="";
					if ($a < $b){
						$amont=$annotation{$display_id}{"bioseq"}->trunc($a,$b)->revcom->translate->seq	;	
					}
					if ($amont=~m/(.*\*)(.*)/){
						$amont=$2;				
						$pos_small=length($annotation{$display_id}{"bioseq"}->seq)-(length($1)*3)-$hsp->frame("query");
					}else{
						$pos_small=length($annotation{$display_id}{"bioseq"}->seq)-$hsp->frame("query");
					}					
				}		
				my $bool=0;
				my $bool_get_firstM="";
				$bool_signal=0;
				while ($amont=~m/(.*?)M/g){
					my $amont_M=$1;
					if ($amont_M ne ""){
						$amont_M_20nt="";
						$start+=(length($amont_M)*3);
						my $add=20;
						if ( ($hsp->strand == 1)){	
							if (($pos_small+$start)-$add < 0){
								$add=$pos_small+$start-1;
							}	
							$amont_M_20nt=$annotation{$display_id}{"bioseq"}->trunc(($pos_small+$start)-$add,($pos_small+$start)-1)->seq;
						}else{
							if ((($pos_small-$start)+$add) > $annotation{$display_id}{"bioseq"}->length){
								$add=$annotation{$display_id}{"bioseq"}->length-($pos_small-$start);
							}
							$amont_M_20nt=$annotation{$display_id}{"bioseq"}->trunc(($pos_small-$start+1),(($pos_small-$start)+$add))->revcom->seq	;
						}
						signal();
						if ($bool_signal==1){
							if ( ($hsp->strand == 1)){	
								$anno_start=($pos_small+$start);
							}else{
								$anno_end=($pos_small-$start);
							}
							$note="Alignment without M but in upstream there is a M with a signal. ".$note_signal;
							$warning=0;				
							last;
						}					
					}
					if ($bool==0){
						if ( ($hsp->strand == 1)){		
							$bool_get_firstM=$pos_small+$start;
						}else{
							$bool_get_firstM=($pos_small-$start);
						}						
					}					
					$bool=1;
				}
				if ($bool_signal==0 && $bool==1){
					if ( ($hsp->strand == 1)){		
						$anno_start=$bool_get_firstM;
					}else{
						$anno_end=$bool_get_firstM;
					}
					$note="Alignment without M but in upstream there is a M without a signal.";
					$warning=1;
				}elsif ($bool_signal==0){					
					if ( ($hsp->strand == 1)){		
						$anno_start=$hsp->start("query");
					}else{
						$anno_end=$hsp->end("query");
					}					
					$note="Alignment without M and there is no M at all upstream of the alignment.";
					$warning=1;
				}
			}
		}
		if ($hsp->strand() == -1){
			$anno_start=$anno_start+3;
		}else{
			$anno_end=$anno_end-3;
		}
		#Verifcation if the CDS is at least out of one extremity
		my $end=$anno_end;	
		my $tmp;
		my $seq_obj=$annotation{$display_id}{"bioseq"};
		my $a=$anno_start;
		my $b=$anno_end;
		if ($hsp->strand()==-1){
			$a-=3;
		}else{
			$b+=3;
		}			
		if ($debug==1){print 'translate_i'.((($b)-($a)+1)%3)."\n";		if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}	}
		if ($hsp->strand() == -1){			
			$tmp=$seq_obj->trunc($a,$b)->revcom->translate("","",0)->seq;
		}else{
			$tmp=$seq_obj->trunc($a,$b)->translate("","",0)->seq;
		}
		my $bool_start=0;
		my $bool_stop=0;
		if ($tmp =~m/\*$/){chop($tmp);$bool_stop=1;}
		if ($tmp=~m/^M/ && $warning==0){
			$bool_start=1;
		}
		my $bool_check_extremity=1;
		if ($hsp->strand()==-1){
			if ($bool_stop==0 || ($orf_start>=($seq_obj->length-3) && $bool_start==0)){
				$bool_check_extremity=0;
			}		
		}else{
			if ($bool_stop==0 || ($orf_start<=3 && $bool_start==0)){	
				$bool_check_extremity=0;
			}
		}			
		if (($bool_check_extremity==0 && ($anno_end - $anno_start+1)/3 <= 120*$hit->length/100) || (($anno_end - $anno_start+1)/3 >= 80*$hit->length/100 && ($anno_end - $anno_start+1)/3 <= 120*$hit->length/100)){
			$small_cds{$display_id}{$anno_start}{"end"}=$anno_end;
			$small_cds{$display_id}{$anno_start}{"note"}=$note;
			if ($debug==1){
				$small_cds{$display_id}{$anno_start}{"note"}.="      length : hit80%<CDS>hit80%)".':'.$hit->length.'  '.(80*$hit->length/100).'<='.(($anno_end - $anno_start+1)/3).'<='.(120*$hit->length/100).'.';
				if ($amont10 ne ""){
					$small_cds{$display_id}{$anno_start}{"note"}.="    The 5AA +/- start(query - hit):".$amont10.".";
				}			
				if ($amont_M_20nt ne ""){
					$small_cds{$display_id}{$anno_start}{"note"}.="    The 20nt upstream:".$amont_M_20nt.".";
				}
				$small_cds{$display_id}{$anno_start}{"note"}.="    orf_start:".$orf_start.".";
				$small_cds{$display_id}{$anno_start}{"note"}.="    orf_stop:".$orf_stop.".";
			}
			$small_cds{$display_id}{$anno_start}{"warning"}=$warning;
			$small_cds{$display_id}{$anno_start}{"hit"}=$hit->name();
			$small_cds{$display_id}{$anno_start}{"evalue"}=$hsp->evalue();		
			$small_cds{$display_id}{$anno_start}{"strand"}=$hsp->strand();
			$small_cds{$display_id}{$anno_start}{"orf"}=$orf_start."-".$orf_stop;
			$orf_hit_small{$orf_stop."-".$hsp->strand}=1;
		}
	}
}

############CDS alignment
$list_dat_name=~s/ /\./g;
#Please ask to jeremy.tournayre@inrae.fr to have the databases
my $db_blast_CDS="Proteomes_".$list_dat_name.".txt";
if (!-e "$dir_db/db_blast/".$db_blast_CDS.".pdb"){
	my @split=split(",",$list_dat_name);
	foreach(@split){
		$_=~s/ /\./g;
		$_="Proteomes_".$_.".txt";
		`echo " " >>$dir_db/db_blast/$db_blast_CDS`;
		`cat $dir_db/db_blast/$_ >>$dir_db/db_blast/$db_blast_CDS`;
		`$dir_blast/makeblastdb -in $dir_db/db_blast/$db_blast_CDS -out $dir_db/db_blast/$db_blast_CDS  -dbtype prot`;
	}
}
my $orf_file_bls=$new_dir_o."/orf.bls";
`$dir_blast/blastp -word_size 3  -num_alignments  20 -matrix BLOSUM45 -evalue $evalue_alignment -db $dir_db/db_blast/$db_blast_CDS -query $orf_file -out $orf_file_bls `;

my $Ech = new Bio::SearchIO(-format => 'blast', -file => $orf_file_bls); 
while ($result = $Ech->next_result){
	my $bool_hit=0;
	while ($hit = $result->next_hit){
		$hsp = $hit->next_hsp;
		annot_alignment($bool_hit);
		$bool_hit=1;
	}
}
############Small CDS
my $file_small_bls=$new_dir_o."/small.bls";
#Please ask to jeremy.tournayre@inrae.fr to have the databases
if (!-e "$dir_db/small_db/small_db.fa.pdb"){
	`$dir_blast/makeblastdb -in $dir_db/small_db/small_db.fa -out $dir_db/small_db/small_db.fa  -dbtype prot`;
}
`$dir_blast/blastx -word_size 3 -num_alignments 500 -seg no -evalue $evalue_small -matrix BLOSUM45 -db  $dir_db/small_db/small_db.fa -query $new_input_file -out $file_small_bls`;
$Ech = new Bio::SearchIO(-format => 'blast', -file => $file_small_bls); 
while ($result = $Ech->next_result){
	while ($hit = $result->next_hit){
		while ($hsp = $hit->next_hsp){
			if ($hsp->query_string =~m/\*/){next;}
			annot_alignment_small();
		}
	}	
}

############Frameshift
foreach (sort keys %hit_2_frameshift){
	my $display_id=$_;
	my %already_done;
	foreach (sort keys %{$hit_2_frameshift{$display_id}}){
		my $strand=$_;
		foreach (sort {$a <=> $b} keys %{$hit_2_frameshift{$display_id}{$strand}}){
			my $order=$_;
			my $orf=$hit_2_frameshift{$display_id}{$strand}{$order}{'orf'};
			if (!defined($already_done{$orf}) && defined($annotation{$display_id}{$orf})){
				$already_done{$orf}=1;
				my %list_hit;
				foreach (keys %{$hit_2_frameshift{$display_id}{$strand}{$order}{'hit'}}){
					$list_hit{$_}=1;
				}
				my %get;

				my @split=split("_",$orf);my @split2=split("-",$split[3]);
				my $s=$split2[0];my $e=$split2[1];

				$get{$orf}=1;	
				while (1){
					my $bool_orf=0;
					foreach (sort {$a <=> $b} keys %{$hit_2_frameshift{$display_id}{$strand}}){
						my $order_b=$_;
						my $orf_b=$hit_2_frameshift{$display_id}{$strand}{$order_b}{'orf'};
						if (!defined($already_done{$orf_b})){
							if ($annotation{$display_id}{$orf_b}{"CDS"}{"note"} =~ m/Alignment starting with M/ || $annotation{$display_id}{$orf_b}{"CDS"}{"note"} =~ m/Alignment with a M into the 5AA before or after the start of the hit/){
								next;
							}			
							if (!defined($get{$orf_b})){
								my $bool=0;
								my @split=split("_",$orf_b);my @split2=split("-",$split[3]);
								my $for_s=$split2[0];my $for_e=$split2[1];	
								foreach (keys %{$hit_2_frameshift{$display_id}{$strand}{$order_b}{'hit'}}){
									if (defined($list_hit{$_}) && $for_s <= ($e+50) ){
										$bool=1;
										if ($s>$for_s){
											$s=$for_s;
											$bool_orf=1;
										}
										if ($e<$for_e){
											$e=$for_e;
											$bool_orf=1;
										}
										$get{$orf_b}=1;
									}
									if ($bool==1){
											last;
									}
								}
								if ($bool==1){
									foreach (keys %{$hit_2_frameshift{$display_id}{$strand}{$order_b}{'hit'}}){
										if (!defined($list_hit{$_})){
											$list_hit{$_}=1;
											$bool_orf=1;
										}
									}		
								}
							}
						}
					}
					if ($bool_orf==0){
						last;
					}
				}
				
				my $i=0;
				my @list_orf;
				foreach (sort keys %get){
					push(@list_orf,$_);
					$i++;
					$already_done{$_}=1;
				}
				if ($i>1){
					my $bool_strand=0;
					my $anno_start=0;
					my $anno_end=0;
					my $note="";
					my $hits="";
					my $anno_orf="";
					my @split;
					my @split2;
					my $evalues="";
					my $bool_stop="";
					foreach (@list_orf){
						my $orf=$_;
						$anno_orf.=$orf.",";
						@split=split("_",$orf);
						@split2=split("-",$split[3]);
						if ($anno_start == 0 || $anno_start>$annotation{$display_id}{$orf}{"CDS"}{"start"}){
							$anno_start=$annotation{$display_id}{$orf}{"CDS"}{"start"};
							if ($split[2] !~m/^-/){	
								$note="Frameshift or intron. ".$annotation{$display_id}{$orf}{"CDS"}{"note"};
							}else{
								$a=$annotation{$display_id}{$orf}{"CDS"}{"start"}-3;
								$b=$annotation{$display_id}{$orf}{"CDS"}{"end"};
								my $seq_obj=$annotation{$display_id}{"bioseq"};
								if ($a <=0){$a=$a+3;}
								my $tmp=$seq_obj->trunc($a,$b)->revcom->translate("","",0)->seq;
								if ($tmp =~m/\*$/){$bool_stop=1;}else{$bool_stop=0;}
							}								
						}
						if ($anno_end == 0 || $anno_end<$annotation{$display_id}{$orf}{"CDS"}{"end"}){
							$anno_end=$annotation{$display_id}{$orf}{"CDS"}{"end"};
							if ($split[2] =~m/^-/){	
								$note="Frameshift or intron. ".$annotation{$display_id}{$orf}{"CDS"}{"note"};
							}else{
								$a=$annotation{$display_id}{$orf}{"CDS"}{"start"};
								$b=$annotation{$display_id}{$orf}{"CDS"}{"end"}+3;
								my $seq_obj=$annotation{$display_id}{"bioseq"};
								if ($b >$seq_obj->length){$b=$b-3;}
								my $tmp=$seq_obj->trunc($a,$b)->translate("","",0)->seq;
								if ($tmp =~m/\*$/){$bool_stop=1;}else{$bool_stop=0;}
							}
						}
						$evalues.=$annotation{$display_id}{$orf}{"CDS"}{"evalue"}." ";
						$hits.=$annotation{$display_id}{$orf}{"CDS"}{"hit"}." ";
						delete $annotation{$display_id}{$orf};
					}
					chop($evalues);
					chop($anno_orf);
					chop($hits);
					$annotation{$display_id}{$anno_orf}{"CDS"}{"frame_or_intron"}=$bool_stop;	
					$annotation{$display_id}{$anno_orf}{"CDS"}{"start"}=$anno_start;
					$annotation{$display_id}{$anno_orf}{"CDS"}{"end"}=$anno_end;
					$annotation{$display_id}{$anno_orf}{"CDS"}{"note"}=$note;
					$annotation{$display_id}{$anno_orf}{"CDS"}{"hit"}="multiple: ".$hits;
					$annotation{$display_id}{$anno_orf}{"CDS"}{"evalue"}="multiple: ".$evalues;
					$annotation{$display_id}{$anno_orf}{"CDS"}{"warning"}=1;
				}
			}
		}
	}
}

############ORF without alignment
my $learn_Glim_orf_file;
$learn_Glim_orf_file=$new_dir_o."/".'learn_Glim_orf.fa';
open (O_CDS_glimmer,'>'.$learn_Glim_orf_file);
my $nb_CDS=0;
my $CDS_orf_file=$new_dir_o."/".'CDS_homo_orf.fa';

foreach (keys %annotation){
	my $display_id=$_;
	foreach (keys %{$annotation{$display_id}}){
		my $orf=$_;
		if ($orf eq "bioseq"){
			next;
		}
		if (defined($annotation{$display_id}{$orf}{"CDS"})){
			my @split=split("_",$orf);
			my @split2=split("-",$split[3]);
			my $seq="";
			if ($split[2] =~m/^-/){
				$seq=$annotation{$display_id}{"bioseq"}->trunc($annotation{$display_id}{$orf}{"CDS"}{"start"},$annotation{$display_id}{$orf}{"CDS"}{"end"})->revcom->seq;	
			}else{
				$seq=$annotation{$display_id}{"bioseq"}->trunc($annotation{$display_id}{$orf}{"CDS"}{"start"},$annotation{$display_id}{$orf}{"CDS"}{"end"})->seq;
			}
			if ($annotation{$display_id}{$orf}{"CDS"}{"warning"}==0){
				print O_CDS_glimmer '>'.$orf."-".$annotation{$display_id}{$orf}{"CDS"}{"hit"}."\n";
				print O_CDS_glimmer $seq."\n";
				$nb_CDS++;
			}		
		}
	}	
}
close (O_CDS_glimmer);
############TrnaSCAN-se
my $trna_file=$new_dir_o."/".'trna.txt';
if (-e $trna_file){
	`rm $trna_file`;
}
`tRNAscan-SE -q --score 50 -o $trna_file $new_input_file`;
open (F, $trna_file);
my $skip=<F>;$skip=<F>;$skip=<F>;
my %trna;
while (<F>){
	chomp($_);
	my @split=split('\s+',$_);
	$trna{$split[0]}{$split[2]}{$split[3]}=1;
}
close (F);

############Ribosomal RNA
my $rRNA_bls=$new_dir_o."/rRNA.bls";
#Please ask to jeremy.tournayre@inrae.fr to have the databases
if (!-e "$dir_db/db_blast/16S_rRNA.txt.nin"){
	`$dir_blast/makeblastdb -in $dir_db/db_blast/16S_rRNA.txt -out $dir_db/db_blast/16S_rRNA.txt  -dbtype nucl`;
}
`$dir_blast/blastn -num_alignments  500 -word_size 7  -gapopen 4  -gapextend 2 -penalty -1 -reward 1 -evalue 1e-50 -db $dir_db/db_blast/16S_rRNA.txt -query $new_input_file -out $rRNA_bls `;
$Ech = new Bio::SearchIO(-format => 'blast', -file => $rRNA_bls); 
my %rrna;
while ($result = $Ech->next_result){
	$display_id=$result->query_name();
	while ($hit = $result->next_hit){
		my @tab_pos_neg=(-1,1);
		my %start;
		my %end;
		my %temp_rrna;
		while ($hsp = $hit->next_hsp){
			my $strand=$hsp->strand("subject");
			my $start=$hsp->start("query");
			my $end=$hsp->end("query");
			if ($strand == -1){
				$start= $start-2500;
			}else{
				$end= $end+2500;
			}
			if ($start < 1){
				$start=1;
			}
			if ($end > $annotation{$display_id}{"bioseq"}->length){
				$end=$annotation{$display_id}{"bioseq"}->length;
			}					
			$rrna{$display_id}{$start}{$end}{"c"}=$strand;
		}
		last;
	}
}

############Glimmer
if ($nb_CDS >= 50){
	$learn_Glim_orf_file_icm=$new_dir_o."/learn_Glim_orf.icm";
	`$dir_glimmer/build-icm -r $learn_Glim_orf_file_icm < $learn_Glim_orf_file`;
}
my $res_glimmer_file=$new_dir_o."/".'res_glimmer';

`$dir_glimmer/glimmer3 -X -A atg -l -o10 -g $glimmer_size -t 30 $new_input_file $learn_Glim_orf_file_icm $res_glimmer_file`;

my %glim;


open (F,$res_glimmer_file.".predict");
$display_id="";
while (<F>){
	chomp($_);
	if ($_ =~m/^>(.*)/){
		$display_id=$1;

	}else{
		my @split=split('\s+',$_);
		$bool_signal=0;	
		my $complement=1;
		my $size_check_AA=30;
		my $i=0;
		################Problem of size
		if ($split[1] > $annotation{$display_id}{"bioseq"}->length()){
			$split[1]=$split[1]-3;
		}
		if ($split[2] > $annotation{$display_id}{"bioseq"}->length()){
			$split[2]=$split[2]-3;
		}				
		if ($split[1] <= 0){
			$split[1]=$split[1]+3;
		}	
		if ($split[2] <= 0){
			$split[2]=$split[2]+3;
		}
		################END OF Problem of size
		my $strand=1;
		if ($split[3]=~m/-/){
			$strand=-1;
		}
		$split[3]=~s/-//;
		$split[3]--;
		
		$orf_stop="";
		$orf_start="";
		if ($strand == 1){
			find_orf($strand,$split[3],$split[1],$split[2]);
		}else{
			$split[3]=($annotation{$display_id}{"bioseq"}->length()-$split[2]+1)%3;
			find_orf($strand,$split[3],$split[2],$split[1]);
		}
	
		$note_signal="";
		if ($split[2]<$split[1]){
			$complement=-1;
			my $a=$split[2];
			my $b=$orf_start;
			if ($debug==1){print 'translate_j'.((($b)-($a)+1)%3)."\n";		if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
			my $prot=$annotation{$display_id}{"bioseq"}->trunc($a,$b)->revcom->translate->seq;
			my $check_AA=$prot;
			if (length($prot)>=(($b-$split[1])/3)+$size_check_AA+1){
				$check_AA=substr($prot,0,(($b-$split[1])/3)+$size_check_AA+1);
			}				
			while ($check_AA=~m/(.*?)M/g){
				$i+=length($1)*3;

				my $a=$split[2];
				my $b=$orf_start;
				$b=$b-$i;
				my $end=$b+21;
				if ($end>$annotation{$display_id}{"bioseq"}->length()+1){
					$end=$annotation{$display_id}{"bioseq"}->length()+1;
				}
				if (($b+1) >= $end-1){
					$amont_M_20nt="";
				}else{
					$amont_M_20nt=$annotation{$display_id}{"bioseq"}->trunc(($b+1),(($end-1)))->revcom->seq;
				}
				$note_signal="";
				if ($amont_M_20nt ne ""){
					signal();
					if ($i < $orf_start-$split[1] && $note_signal=~m/weak/){
						$bool_signal=0;
					}
				}
				if ($bool_signal==1){
					$note_signal="M with a signal ".$note_signal;
					last;
				}else{
					$note_signal="M without a signal";
				}	
				$i+=3;				
			}

			if ($bool_signal==1){
				# $i=$i+3;
				if((($i-($orf_start-$split[1]))/3)!=0){
					if ($i < $orf_start-$split[1]){
						$note_signal.=". The selected M is at ".(-($i-($orf_start-$split[1]))/3)." AA before the M selected by Glimmer because it had a signal";
					}else{
						$note_signal.=". The selected M is at ".(($i-($orf_start-$split[1]))/3)." AA after the M selected by Glimmer because it had a signal";
					}
					$split[1]=$orf_start-$i;	
				}		
			}
			my $tmp=$split[1];
			$split[1]=$split[2];
			$split[2]=$tmp;	
			$split[1]=$split[1]+3;							
		}else{
			my $a=$orf_start;
			my $b=$split[2];
			if ($debug==1){print 'translate_k'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}						
			my $prot=$annotation{$display_id}{"bioseq"}->trunc($a,$b)->translate->seq;
			my $check_AA=$prot;
			if (length($prot)>=(($split[1]-$a)/3)+$size_check_AA+1){
				$check_AA=substr($prot,0,(($split[1]-$a)/3)+$size_check_AA+1);
			}
			while ($check_AA=~m/(.*?)M/g){
				$i+=length($1)*3;
				my $a=$orf_start;
				my $b=$split[2];
				$a=$a+$i;
				my $start=($a-21);
				if ($start<0){
					$start=0;
				}
				if ($start+1 >= ($a-1)){
					$amont_M_20nt="";
				}else{
					$amont_M_20nt=$annotation{$display_id}{"bioseq"}->trunc(($start+1),$a-1)->seq;
				}
				$note_signal="";
				if ($amont_M_20nt ne ""){
					signal();
					if ($i < $split[1]-$orf_start && $note_signal=~m/weak/){
						$bool_signal=0;
					}					
				}
				if ($bool_signal==1){
					$note_signal="M with a signal ".$note_signal;
					last;
				}else{
					$note_signal="M without a signal";
				}	
				$i+=3;
			}
			if ($bool_signal==1){
				if((($i-($split[1]-$orf_start))/3)!=0){
					if ($i < $split[1]-$orf_start){
						$note_signal.=". The selected M is at ".(-($i-($split[1]-$orf_start))/3)." AA before the M selected by Glimmer because it had a signal";
					}else{
						$note_signal.=". The selected M is at ".(($i-($split[1]-$orf_start))/3)." AA after the M selected by Glimmer because it had a signal";
					}
				}
				$split[1]=$orf_start+$i;		
			}
			$split[2]=$split[2]-3;		
		}		
		if ($note_signal eq ""){
			$note_signal="Prediction without M";
		}
		$glim{$display_id}{$split[1]}{$split[2]}{"c"}=$complement;
		$glim{$display_id}{$split[1]}{$split[2]}{"s"}=$note_signal;
		$glim{$display_id}{$split[1]}{$split[2]}{"o"}=$orf_start."-".$orf_stop;
	}
}
close(F);

############ OUTFILES
#######Create CDS file for TE and Interproscan
my $file_CDS_gene_nt=$new_dir_o."/CDS_gene_nt.fa";
open (FILE_CDS_GENE_NT,'>'.$file_CDS_gene_nt);
my %list_file_CDS_aa;
my $feat;
my $cds_start_end_end;
my $bool_display_id_infile_warnings;
my $display_id_infile;
my $tag;
sub start_end_extremity{
	my $tmp=$_[0];
	my $strand=$_[1];
	my $start=$_[2];	
	my $end=$_[3];
	my $length=$_[4];
	my $translation=$_[5];
	my $start_ORF=$_[6];
	my $end_ORF=$_[7];
	my $note="";
	if (defined($_[8])){
		$note=$_[8];
	}
	my $bool_start=0;
	my $bool_stop=0;
	my $return="";
	$cds_start_end_end=$start."-".$end."\n";
	if ($translation ne "no"){
		if ($translation eq "frameshift_stop"){$bool_stop=1;}
		elsif($translation eq "frameshift_nostop"){$bool_stop=0;}
		elsif ($tmp =~m/\*$/){chop($tmp);$bool_stop=1;}
		if ($tmp=~m/^M/){$bool_start=1;}	
		if ($translation eq "all"){
			$feat->add_tag_value('translation',$tmp);
		}
		if ($debug==1 && $note =~m/there is a M/ && $bool_start==0){
			print "verif: sequence dont begin by M\n";
		}		
		if ($strand==-1){
			my $bool_extremity=0;
			if ($bool_stop==0){
				$bool_extremity=1;
				$feat->start("<1");
				$feat->add_tag_value('note',"stop undetermined");
				$cds_start_end_end="<1-".$end."\n";
			}else{
				$feat->start($start);
			}		
			if ($end>=($length-3) && $bool_start==0){
				$feat->end(">".$end);
				$cds_start_end_end=$start."-".">".$end."\n";
			}else{
				$feat->end($end);
			}	
			if ($bool_extremity || ($start_ORF>($length-3))){
				$feat->add_tag_value('note',"extremity");
			}
			if (($start_ORF-$end) > 300){
				my $warning_start_far="Warning: the start is at ".($start_ORF-$end)." nucleotides from the start of the ORF";
				$feat->add_tag_value('note',$warning_start_far);
				if ($bool_display_id_infile_warnings==0){
					print FILE_WARNINGS '>'.$display_id_infile."\n";
				}
				$bool_display_id_infile_warnings=1;
				print FILE_WARNINGS $tag."-".$cds_start_end_end.$warning_start_far."\n";
			}			
		}else{
			my $bool_extremity=0;
			if ($start<=3 && $bool_start==0){
				$feat->start("<1");
				$cds_start_end_end="<1"."-".$end."\n";
			}else{
				$feat->start($start);
			}
			if ($bool_stop==0){
				$bool_extremity=1;
				$feat->end(">".$length);
				$feat->add_tag_value('note',"stop undetermined");
				$cds_start_end_end=$start."-".">".$length."\n";
			}else{
				$feat->end($end);
			}
			if ($bool_extremity || ($start_ORF<=3)){
				$feat->add_tag_value('note',"extremity");
			}
			if (($start-$start_ORF) > 300){
				my $warning_start_far="Warning: the start is at ".($start-$start_ORF)." nucleotides from the start of the ORF";
				$feat->add_tag_value('note',$warning_start_far);
				if ($bool_display_id_infile_warnings==0){
					print FILE_WARNINGS '>'.$display_id_infile."\n";
				}
				$bool_display_id_infile_warnings=1;
				print FILE_WARNINGS $tag."-".$cds_start_end_end.$warning_start_far."\n";
			}			
		}
	}else{
		$feat->start($start);
		$feat->end($end);
	}
}



my $dir_gb=$new_dir_o."/Res_gb";
if (-e $dir_gb){`rm -R $dir_gb`;}
`mkdir $dir_gb`;
my $dir_embl=$new_dir_o."/Res_embl";
if (-e $dir_embl){`rm -R $dir_embl`;}
`mkdir $dir_embl`;
my $dir_warnings=$new_dir_o."/Res_warnings";
if (-e $dir_warnings){`rm -R $dir_warnings`;}
`mkdir $dir_warnings`;

#modification Bioperl: add sort on row 1264 on genbank.pm (usr/share/perl5/Bio/SeqIO)
    # foreach my $tag (sort keys %{ $fth->field } ) { 
	# instead of     
	# foreach my $tag (keys %{ $fth->field } ) {

# add sort on row 956 on embl.pm  (usr/share/perl5/Bio/SeqIO)
# foreach my $tag (sort keys %{$fth->field} ) {
# instead of     
# foreach my $tag (keys %{ $fth->field } ) {
	
# 'eq' instead of '==' on row 350 on Simple.pm (usr/share/perl5/Bio/Location)

my %list_name;
foreach (sort {$a<=>$b} keys %annotation){
	my $display_id=$_;
	my $seq_obj = $annotation{$display_id}{"bioseq"};
	my $tmp=$id_2_name{$display_id};
	$tmp=~s/ /_/g;
	$tmp=~s/>//g;
	$tmp=substr($tmp,0,20);
	my $tmp_file=$tmp;
	my $i_file=2;
	while (defined($list_name{$tmp_file})){
		$tmp_file=$tmp."_".$i_file;
		$i_file++;
	}
	$list_name{$tmp_file}=1;
	my $out = Bio::SeqIO ->new (-file => '>'.$dir_gb.'/'.$tmp_file.'.gb',-format => 'genbank'); 
	my $out2 = Bio::SeqIO ->new (-file => '>'.$dir_embl.'/'.$tmp_file.'.embl',-format => 'embl');
	my $file_warnings=$dir_warnings.'/'.$tmp_file.'.txt';
	
	open (FILE_WARNINGS,'>'.$file_warnings);
	$display_id_infile=$tmp;
	$bool_display_id_infile_warnings=0;
	$seq_obj->display_id($display_id_infile);
	#Order of annotation
	my %order;
	my %cds_start_end;
	foreach (sort keys %{$annotation{$display_id}}){
		my $orf=$_;
		if ($orf ne "bioseq"){
			if (defined($annotation{$display_id}{$orf}{"CDS"})){
				if ($debug==1 && $annotation{$display_id}{$orf}{"CDS"}{"start"} > $annotation{$display_id}{$orf}{"CDS"}{"end"}){
					print '1_verif_cds_start_end'."\n"; 
				}
				$order{$annotation{$display_id}{$orf}{"CDS"}{"start"}}=$orf;
				$cds_start_end{$annotation{$display_id}{$orf}{"CDS"}{"start"}}{$annotation{$display_id}{$orf}{"CDS"}{"end"}}=1;
			}
		}
	}

	#small_cds : only if there is no overlapping CDS already defined
	my %tmp_cds_start_end = %{ clone \%cds_start_end };
	foreach (sort {$a<=>$b} keys %{$small_cds{$display_id}}){
		my $anno_start=$_;
		my $anno_end=$small_cds{$display_id}{$anno_start}{"end"};
		my $bool_print=1;
		foreach (sort {$a<=>$b} keys %tmp_cds_start_end){	
			my $CDS_start=$_;
			if ($anno_end < $CDS_start){
				last;
			}			
			foreach (sort {$a<=>$b} keys %{$tmp_cds_start_end{$CDS_start}}){	
				my $CDS_end=$_;
				if ($anno_start <= $CDS_end && $anno_end >= $CDS_start){
					$bool_print=0;
					last;
				}
				if ($anno_start > $CDS_end){
					delete($tmp_cds_start_end{$CDS_start}{$CDS_end});
				}									
			}
		}	
		if ($bool_print==1){
			if ($debug==1 && $anno_start > $anno_end){
				print '2_verif_cds_start_end'."\n";
			}			
			$order{$anno_start}="small_cds";	
			$cds_start_end{$anno_start}{$anno_end}=1;
		}
	}

	#tRNA :
	my %others_start_end;
	foreach (sort {$a<=>$b} keys %{$trna{$display_id}}){
		my $anno_start=$_;
		foreach (sort {$a<=>$b} keys %{$trna{$display_id}{$anno_start}}){
			my $anno_end=$_;
			my $order_start=$anno_start;
			my $order_end=$anno_end;
			if ($anno_start > $anno_end){
				$order_start=$anno_end;
				$order_end=$anno_start;
			}
			$order{$order_start}="trna:".$anno_start.'-'.$anno_end;		
			$others_start_end{$order_start}{$order_end}=1;
		}
	}

	#rRNA : 
	foreach (sort {$a<=>$b} keys %{$rrna{$display_id}}){
		my $anno_start=$_;
		foreach (sort {$a<=>$b} keys %{$rrna{$display_id}{$anno_start}}){
			my $anno_end=$_;
			my $order_start=$anno_start;
			my $order_end=$anno_end;
			if ($anno_start > $anno_end){
				$order_start=$anno_end;
				$order_end=$anno_start;
			}
			$order{$order_start}="rrna:".$anno_start.'-'.$anno_end;		
			$others_start_end{$order_start}{$order_end}=1;
		}
	}

	#Glimmer : only if there is no overlapping CDS already defined
	my %get;
	my $low_start;
	my $max_end=0;
	my $i=0;
	my %to_check;
	foreach (sort {$a <=> $b} keys %{$glim{$display_id}}){
		my $start=$_;
		foreach (sort {$a <=> $b} keys %{$glim{$display_id}{$start}}){
			my $end=$_;

			my $bool_print=1;
			foreach (sort {$a<=>$b} keys %others_start_end){	
				my $CDS_start=$_;
				if ($end < $CDS_start){
					last;
				}
				foreach (sort {$a<=>$b} keys %{$others_start_end{$CDS_start}}){	
					my $CDS_end=$_;
					if ($start <= $CDS_end && $end >= $CDS_start){
						$bool_print=0;
						last;
					}
					if ($start > $CDS_end){
						delete($others_start_end{$CDS_start}{$CDS_end})
					}									
				}
			}
			if ($bool_print==1){
				if ($max_end==0){
					$max_end=$end;
					$low_start=$start;
				}
				if ($start<=$max_end && $end >= $low_start){
					if ($end>$max_end ){
						$max_end=$end;
					}
				}else{
					$max_end=$end;
					$low_start=$start;
					$i++;
				}
				$get{$i}{"l_m"}=$low_start."-".$max_end;
				$get{$i}{"list"}{$end-$start}{$start.'-'.$end}=1;
				$to_check{$i}++;
			}
		}
	}
	foreach (sort {$a <=> $b} keys %get){
		my $i=$_;
		my @split=split("-",$get{$i}{"l_m"});
		my $start=$split[0];
		my $end=$split[1];
		foreach (sort {$a <=> $b} keys %cds_start_end){
			my $CDS_start=$_;
			foreach (sort {$a <=> $b} keys %{$cds_start_end{$CDS_start}}){
				my $CDS_end=$_;
				if ($start<=$CDS_end && $end >= $CDS_start){
					my $get_max=0;
					foreach (sort {$b <=> $a} keys %{$get{$i}{"list"}}){
						$get_max=$_;
						last;
					}
					$get{$i}{"list"}{($get_max+1)}{$CDS_start."-".$CDS_end}="alignment";
					$to_check{$i}++;
				}
			}
		}
	}
	foreach (sort {$a <=> $b} keys %get){
		my $i=$_;
		while ($to_check{$i}){
			my $original_alignment=0;
			my $original_start=0;
			my $original_end;	
			my $max_end=0;
			my $low_start;			
			my $first_small_o_big;
			foreach (sort {$b <=> $a} keys %{$get{$i}{"list"}}){
				my $size=$_;
				my $small_o_big="small";
				if ($size>500){
					$small_o_big="big";
				}		
				foreach (sort keys %{$get{$i}{"list"}{$size}}){
					my @split=split("-",$_);
					my $start=$split[0];
					my $end=$split[1];
					if ($max_end==0){
						$max_end=$end;
						$low_start=$start;
						$first_small_o_big=$small_o_big;
						$original_start=$start;
						$original_end=$end;							
						if ($get{$i}{"list"}{($size)}{$start."-".$end} eq 1){
							$order{$start}="Glimmer:".$start."-".$end.":".$glim{$display_id}{$start}{$end}{"s"};
						
						}else{
							$first_small_o_big="big";
							$original_alignment=1;
						}
					}
					if ($get{$i}{"list"}{($size)}{$start."-".$end} ne 1){
						$small_o_big="big";
					}
					if ($start<=$max_end && $end >= $low_start){
						my $bool_alignment=$get{$i}{"list"}{($size)}{$start."-".$end} eq 1;
						delete($get{$i}{"list"}{$size}{$start.'-'.$end});
						$to_check{$i}--;
						if ($small_o_big eq "big" && $first_small_o_big eq "big" && (($end>=$low_start && $end<=$max_end && $end-$low_start+1<=50 ) || ($start>=$low_start && $start<=$max_end && $max_end-$start+1<=50))){
							if ($bool_alignment){
								$order{$start}="Glimmer:".$start."-".$end.":".$glim{$display_id}{$start}{$end}{"s"}." Warning: Overlapping with an other predicted CDS\n";
								if ($original_alignment == 0 && $original_start != 0 && $start != $original_start && $end != $original_end){
									$order{$original_start}="Glimmer:".$original_start."-".$original_end.":".$glim{$display_id}{$original_start}{$original_end}{"s"}." Warning: Overlapping with an other predicted CDS\n";

									$original_start=0;
								}elsif($original_alignment){
								}
							}
							if ($start<$low_start ){
								$low_start=$start;
							}
							if ($end>$max_end ){
								$max_end=$end;
							}
						}
					}
				}
			}
		}
	}
	
	foreach (sort {$a<=>$b} keys %order){
		my $a_order=$_;
		$feat = new Bio::SeqFeature::Generic(); 
		my $nt="";my $tmp="";my $a=0;my $b=0;my $start_orf=0;my $end_orf=0;
		if (!defined($annotation{$display_id}{$order{$a_order}}{"CDS"}) && $order{$a_order} ne "small_cds"){
			my @split=split(":",$order{$a_order});
			my @split2=split("-",$split[1]);
			my $mode_translate="all";
			my $complement="+";
			my $note="";
			$tag="";
			if ($split[0] eq "Glimmer"){
				$tag="CDS";
				$feat->primary_tag('CDS'); 
				$note="CDS found by ".$split[0].". ".$split[2];
				if (defined($split[3])){
					$note.=":".$split[3];
				}
				$feat->add_tag_value('note',$note); 
				my @split3=split("-",$glim{$display_id}{$split2[0]}{$split2[1]}{"o"});
				$start_orf=$split3[0];
				$end_orf=$split3[1];
				if ($glim{$display_id}{$split2[0]}{$split2[1]}{"c"}==-1){
					$complement="-";
				}
			}elsif ($split[0] eq "trna"){
				$feat->primary_tag('tRNA'); 	
				$feat->add_tag_value('note',"tRNA found by tRNAscan-SE"); 
				$mode_translate="no";
				if ($split2[0] > $split2[1]){
					$complement="-";
				}				
			}elsif ($split[0] eq "rrna"){
				$feat->primary_tag('rRNA'); 	
				$feat->add_tag_value('note',"rRNA unit"); 
				$mode_translate="no";
				if ($rrna{$display_id}{$split2[0]}{$split2[1]}{"c"}==-1){
					$complement="-";
				}				
			}					
			if ($complement eq "-"){
				$feat->strand("-");		
				$a=$split2[0]-3;
				if ($a <=0){$a=$a+3;}
				$b=$split2[1];				
				if ($mode_translate eq "all"){	
					if ($debug==1){print 'translate_l'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
					$tmp=$seq_obj->trunc($a,$b)->revcom->translate("","",0)->seq;
					$nt=$seq_obj->trunc($a,$b)->revcom->seq;
				}					
				start_end_extremity($tmp,-1,$a,$b,$seq_obj->length,$mode_translate,$start_orf,$end_orf);
			}else{
				$feat->strand("+");	
				$a=$split2[0];
				$b=$split2[1]+3;	
				if ($b >$seq_obj->length){
					$b=$b-3;
				}						
				if ($mode_translate eq "all"){							
					if ($debug==1){print 'translate_m'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}	
					$tmp=$seq_obj->trunc($a,$b)->translate("","",0)->seq;
					$nt=$seq_obj->trunc($a,$b)->seq;
				}					
				start_end_extremity($tmp,1,$a,$b,$seq_obj->length,$mode_translate,$start_orf,$end_orf);	
			}
			if (defined($split[3])){
				if ($debug==1){
					$feat->add_tag_value('note_warning',1); 
				}		
				if ($bool_display_id_infile_warnings==0){
					print FILE_WARNINGS '>'.$display_id_infile."\n";
				}
				print FILE_WARNINGS $tag."-".$cds_start_end_end.$note."\n";
			}else{
				if ($debug==1){
					$feat->add_tag_value('note_warning',0); 
				}						
			}			
		}	
		elsif( $order{$a_order} eq "small_cds"){
			my $mode_translate="no";
			$tag="";
			if ($small_cds{$display_id}{$a_order}{"note"} =~ m/Alignment without M and there is no M at all upstream of the alignment/){
				if ($debug==1){$feat->primary_tag('small_gene');$tag="small_gene";
				}else{$feat->primary_tag('gene');$tag="gene";}
				$small_cds{$display_id}{$a_order}{"note"}.=" Possible frameshift or intron or N-terminal truncated protein.";
				$mode_translate="extremity";
			}else{
				if ($debug==1){$feat->primary_tag('small_CDS');$tag="small_CDS"; }
				else{$feat->primary_tag('CDS');$tag="CDS";}
				$mode_translate="all";
			}
			
			my $end=$small_cds{$display_id}{$a_order}{"end"};	
			my @split=split("-",$small_cds{$display_id}{$a_order}{"orf"});
			$tmp="";$nt="";$a=0;$b=0;$start_orf=$split[0];$end_orf=$split[1];
			if ($small_cds{$display_id}{$a_order}{"strand"} == -1){
				$a=$a_order-3;
				if ($a <=0){$a=$a+3;}
				$b=$end;				
				if ($mode_translate eq "all" || $mode_translate eq "extremity"){	
					if ($debug==1){print 'translate_n'.((($b)-($a)+1)%3)."\n";	if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
					$tmp=$seq_obj->trunc($a,$b)->revcom->translate("","",0)->seq;
					$nt=$seq_obj->trunc($a,$b)->revcom->seq;
				}					
				
			}else{
				$a=$a_order;
				$b=$end+3;
				if ($b >$seq_obj->length){
					$b=$b-3;
				}	
				if ($mode_translate eq "all"  || $mode_translate eq "extremity"){					
					if ($debug==1){print 'translate_o'.((($b)-($a)+1)%3)."\n";	if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
					$tmp=$seq_obj->trunc($a,$b)->translate("","",0)->seq;
					$nt=$seq_obj->trunc($a,$b)->seq;
				}			
			}			
			start_end_extremity($tmp,$small_cds{$display_id}{$a_order}{"strand"},$a,$b,$seq_obj->length,$mode_translate,$start_orf,$end_orf,$small_cds{$display_id}{$a_order}{"note"});	
			$feat->add_tag_value('hit',$small_cds{$display_id}{$a_order}{"hit"}); 
			$feat->add_tag_value('evalue',$small_cds{$display_id}{$a_order}{"evalue"}); 

			$feat->add_tag_value('note',$small_cds{$display_id}{$a_order}{"note"}); 
			if ($debug==1){
				$feat->add_tag_value('note_warning',$small_cds{$display_id}{$a_order}{"warning"}); 
			}
			if ($small_cds{$display_id}{$a_order}{"warning"} == 1){
				if ($bool_display_id_infile_warnings==0){
					print FILE_WARNINGS '>'.$display_id_infile."\n";
				}
				$bool_display_id_infile_warnings=1;
				print FILE_WARNINGS $tag."-".$cds_start_end_end.$small_cds{$display_id}{$a_order}{"note"}."\n";
			}
			$feat->strand($small_cds{$display_id}{$a_order}{"strand"});

		}else{
			my $orf=$order{$a_order};
			
			my @split=split("_",$orf);
			my @split2=split("-",$split[3]);
			my @split3=split(",",$orf);
			$start_orf=0;
			$end_orf=0;
			my $tmp_end_orf=0;
			
			foreach (@split3){
				my @split4=split("_",$_);
				my @split5=split("-",$split4[3]);
				if ($start_orf==0 || $split5[0] < $start_orf){
					$start_orf=$split5[0];
				}
				if ($end_orf==0 || $split5[1] > $end_orf){
					$end_orf=$split5[1];
				}
			}
			if ($split[2] =~m/^-/){
				$start_orf=$annotation{$display_id}{"bioseq"}->length-$start_orf+1;
				$end_orf=$annotation{$display_id}{"bioseq"}->length-$end_orf+1;
			}			
			if ($tmp_end_orf != 0){
				$end_orf=$tmp_end_orf;
			}
			my $start=$annotation{$display_id}{$orf}{"CDS"}{"start"}; 
			my $end=$annotation{$display_id}{$orf}{"CDS"}{"end"};
			my $mode_translate="";
			$tag="";
			if (defined($annotation{$display_id}{$orf}{"CDS"}{"frame_or_intron"})){
				if ($annotation{$display_id}{$orf}{"CDS"}{"frame_or_intron"} == 1){
					$mode_translate='frameshift_stop';
				}else{
					$mode_translate='frameshift_nostop';
				}
				$feat->primary_tag('gene');
				$tag="gene";
			}elsif ($annotation{$display_id}{$orf}{"CDS"}{"note"} =~ m/Alignment without M and there is no M at all upstream of the alignment/){
				$feat->primary_tag('gene');
				$tag="gene";
				$annotation{$display_id}{$orf}{"CDS"}{"note"}.=" Possible frameshift or intron or N-terminal truncated protein.";
				$mode_translate="extremity";
			}else{
				$mode_translate="all";
				$feat->primary_tag('CDS');		
				$tag="CDS";
			}
			if ($split[2] =~m/^-/){
				$feat->strand("-");
				$tmp="";$nt="";$a=0;$b=0;;
				$a=$annotation{$display_id}{$orf}{"CDS"}{"start"}-3;
				if ($a <=0){$a=$a+3;}
				$b=$annotation{$display_id}{$orf}{"CDS"}{"end"};				
				if ($mode_translate eq "all" || $mode_translate eq "extremity"){
					if ($debug==1){print 'translate_p'.((($b)-($a)+1)%3)."\n";if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}			
					$tmp=$seq_obj->trunc($a,$b)->revcom->translate("","",0)->seq; 
					$nt=$seq_obj->trunc($annotation{$display_id}{$orf}{"CDS"}{"start"},$annotation{$display_id}{$orf}{"CDS"}{"end"})->revcom->seq; 
				}
				start_end_extremity($tmp,-1,$a,$b,$seq_obj->length,$mode_translate,$start_orf,$end_orf,$annotation{$display_id}{$orf}{"CDS"}{"note"});				
			}else{
				$feat->strand("+");
				$tmp="";$nt="";
				$a=$annotation{$display_id}{$orf}{"CDS"}{"start"};
				$b=$annotation{$display_id}{$orf}{"CDS"}{"end"}+3;
				if ($b >$seq_obj->length){
					$b=$b-3;
				}				
				if ($mode_translate eq "all" || $mode_translate eq "extremity"){
					if ($debug==1){print 'translate_q'.((($b)-($a)+1)%3)."\n";	if (((($b)-($a)+1)%3) != 0){print 'warning'."\n";}}
					$tmp=$seq_obj->trunc($a,$b)->translate("","",0)->seq;
					$nt=$seq_obj->trunc($annotation{$display_id}{$orf}{"CDS"}{"start"},$annotation{$display_id}{$orf}{"CDS"}{"end"})->seq;
				}
				start_end_extremity($tmp,1,$a,$b,$seq_obj->length,$mode_translate,$start_orf,$end_orf,$annotation{$display_id}{$orf}{"CDS"}{"note"});					
			}

			if (defined($annotation{$display_id}{$orf}{"CDS"}{"hit"})){
				$feat->add_tag_value('hit',$annotation{$display_id}{$orf}{"CDS"}{"hit"}); 
				$feat->add_tag_value('evalue',$annotation{$display_id}{$orf}{"CDS"}{"evalue"}); 
			}
			$feat->add_tag_value('note',$annotation{$display_id}{$orf}{"CDS"}{"note"});
			if ($debug==1){			
				$feat->add_tag_value('note_warning',$annotation{$display_id}{$orf}{"CDS"}{"warning"}); 
			}
			if ($annotation{$display_id}{$orf}{"CDS"}{"warning"} == 1){
				if ($bool_display_id_infile_warnings==0){
					print FILE_WARNINGS '>'.$display_id_infile."\n";
				}
				$bool_display_id_infile_warnings=1;				
				print FILE_WARNINGS $tag."-".$cds_start_end_end.$annotation{$display_id}{$orf}{"CDS"}{"note"}."\n";
			}
		}
		if ($order{$a_order} eq "small_cds" || defined($annotation{$display_id}{$order{$a_order}}{"CDS"}) || $order{$a_order} =~m/^Glimmer:/){
			my $start=$feat->start;
			$start =~ s/[^0-9]//g;
			my $end=$feat->end;
			$end =~ s/[^0-9]//g;
			my $get_nt=$nt;
			if ($get_nt eq ""){
				if ($feat->strand == 1){
					$get_nt=$seq_obj->trunc($a,$b)->seq;
					$tmp=$seq_obj->trunc($a,$b)->translate("","",0)->seq; 
				}else{
					$get_nt=$seq_obj->trunc($a,$b)->revcom->seq;
					$tmp=$seq_obj->trunc($a,$b)->revcom->translate("","",0)->seq; 					
				}
			}
			print FILE_CDS_GENE_NT '>'.$start.'-'.$end."\n".$get_nt."\n";
			if ($tmp =~m/\*$/){chop($tmp);}
			if ($tmp !~ m/\*/){
				$list_file_CDS_aa{$start.'-'.$end}=$tmp;
			}
		}
		
		$seq_obj->add_SeqFeature($feat); 
	}
	$out->write_seq($seq_obj);
	$out2->write_seq($seq_obj);
}
close(FILE_CDS_GENE_NT);

############Transposable Element TE 
my $file_CDS_aa=$new_dir_o."/CDS_aa.fa";
my $TE_bls=$new_dir_o."/TE.bls";

#Please ask to jeremy.tournayre@inrae.fr to have the databases
if (!-e "$dir_db/db_blast/ConsensusTE.txt.nin"){
	`$dir_blast/makeblastdb -in $dir_db/db_blast/ConsensusTE.txt -out $dir_db/db_blast/ConsensusTE.txt  -dbtype nucl`;
}
`$dir_blast/tblastx -matrix BLOSUM45 -word_size 3 -num_alignments  10 -evalue $evalue_TE -db $dir_db/db_blast/ConsensusTE.txt -query $file_CDS_gene_nt -out $TE_bls`;

$Ech = new Bio::SearchIO(-format => 'blast', -file => $TE_bls); 
my %id_2_te;
while ($result = $Ech->next_result){
	while ($hit = $result->next_hit){
		my $name=$result->query_name;
		my $hit_name=$hit->name();
		if (!defined($id_2_te{$name})){
			$id_2_te{$name}="TE found by alignment: ".$hit_name;
		}else{
			$id_2_te{$name}.=", ".$hit_name;
		}
	}
}
open (FILE_CDS_AA,'>'.$file_CDS_aa);
foreach (sort keys %list_file_CDS_aa){
	print FILE_CDS_AA '>'.$_."\n".$list_file_CDS_aa{$_}."\n";
}
close(FILE_CDS_AA);

############Interproscan
my $dir_interpro=$new_dir_o."/Res_Interproscan";
if (-e $dir_interpro){
	`rm -R $dir_interpro`;
}
`mkdir $dir_interpro`;
if ($bool_interpro){
	`/my_interproscan/interproscan-5.60-92.0/interproscan.sh --goterms --iprlookup --output-dir $dir_interpro -i $file_CDS_aa -f xml`;
}

my %id_2_interproscan;
my %convert;
$convert{"BIOLOGICAL_PROCESS"}="BP";
$convert{"MOLECULAR_FUNCTION"}="MF";
$convert{"CELLULAR_COMPONENT"}="CC";
my %id_list;
if (-e $dir_interpro."/CDS_aa.fa.xml"){
	open (F,$dir_interpro."/CDS_aa.fa.xml");
	while (<F>){
		chomp($_);
		if($_ =~m/<sequence md5=/){
			%id_list=();
		}		
		if($_ =~m/^^\s+\<xref\sid\=\"(\d+-\d+)/){
			$id_list{$1}=1;
		}
		if(/^\s+\<entry\s+ac\=\"(\w+)\"\s+desc\=\"(.+?)\"\s+/){
			foreach (keys %id_list){
				$id_2_interproscan{$_}{"function"}{$2}=1;
				$id_2_interproscan{$_}{"db_xref"}{$1}=1;
			}
		}
		if(/^\s+\<go\-xref\s+category\=\"\w+\"\s+db\=\"GO"\s+id\=\"(\w+\:\w+)\"/){ 
			foreach (keys %id_list){
				$id_2_interproscan{$_}{"db_xref"}{$1}=1;
			}
		}
		if(/^\s+\<go\-xref\scategory\=\"(\w+\_\w+)\"\s+db\=\"GO\"\s+id\=\"\w+\:\w+\"\s+name\=\"(.+?)\"/){ 
			foreach (keys %id_list){
				$id_2_interproscan{$_}{"function"}{$2}=1;
			}
		}
	}
	close(F);
}

my $dir_gb_annot=$new_dir_o."/Res_gb_annot";
if (-e $dir_gb_annot){`rm -R $dir_gb_annot`;}
`mkdir $dir_gb_annot`;
my $dir_embl_annot=$new_dir_o."/Res_embl_annot";
if (-e $dir_embl_annot){`rm -R $dir_embl_annot`;}
`mkdir $dir_embl_annot`;
my $dir_gff_annot=$new_dir_o."/Res_gff_annot";
if (-e $dir_gff_annot){`rm -R $dir_gff_annot`;}
`mkdir $dir_gff_annot`;
my @ls=`ls $dir_gb`;
foreach (@ls){
	chomp($_);
	$_=~/(.*).gb/;
	my $wo_ext=$1;
	$in = Bio::SeqIO ->new (-file => $dir_gb.'/'.$_,-format => 'genbank'); 
	my $out = Bio::SeqIO ->new (-file => '>'.$dir_gb_annot.'/'.$_,-format => 'genbank');
	my $out2 = Bio::SeqIO ->new (-file => '>'.$dir_embl_annot.'/'.$wo_ext.'.embl',-format => 'embl');
	open(O_GFF,'>'.$dir_gff_annot.'/'.$wo_ext.'.gff');
	while (my $seq_obj=$in->next_seq()) {
		my $new_seq_obj =  clone $seq_obj;
		$new_seq_obj->remove_SeqFeatures();
		print O_GFF "#".$seq_obj->display_id."\n";
		for my $feat_object ($seq_obj->get_SeqFeatures) {
			my $a=$feat_object->start.'-'.$feat_object->end;
			if (defined($id_2_te{$a})){
				$feat_object->primary_tag('mobile_element');
				foreach ($feat_object->get_all_tags()){
					$feat_object->remove_tag($_);
				}
				$feat_object->add_tag_value("note","Putative transposable element");
			}elsif (defined($id_2_interproscan{$a})){
				foreach (sort keys %{$id_2_interproscan{$a}}){
					my $b=$_;
					foreach (sort keys %{$id_2_interproscan{$a}{$b}}){
						my $c=$_;
						$feat_object->add_tag_value($b,$c); 
					}
				}
			}
			print O_GFF $feat_object->gff_string."\n";
			$new_seq_obj->add_SeqFeature($feat_object); 
		}
		$out->write_seq($new_seq_obj);
		$out2->write_seq($new_seq_obj);
	}
	close(O_GFF);
}
chdir $new_dir_o;
`tar -czf Res_gb_annot.tar.gz Res_gb_annot`;
`tar -czf Res_embl_annot.tar.gz Res_embl_annot`;
`tar -czf Res_gff_annot.tar.gz Res_gff_annot`;
`tar -czf Res_warnings.tar.gz Res_warnings`;