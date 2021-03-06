#!/usr/bin/perl
use strict;
use warnings;
use threads;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
my $expression;
my $knowledge;
my $threads;
my $model = "HFODE";
my $fitness_type = 1;
my $generation_no = 1000;
my $output_dir;
my $tmp_dir;
my $cc = 0;
my $help;
my $program = abs_path($0);
#======================
my %gene_index;
my $min = 0;



sub Usage {
	print STDERR "Usage: perl $program [Options]
Options:
	-i	[FILE]	time-series profile
	-o	[PATH]	Output Directory
	-k	[FILE]	knowledge of regulatory
	-m	[model] {s-system,HFODE} Default is HFODE
	-f	[Fitness] {0,1,2,3,4} Default is 1
	-c	[cc]	{0,1} Default is 0, consider correlation coefficient
	-g	[Generation] Number of generation, default is 1000
	-t	[No]	Number of threads
	-min	[value] ex: 0.01 [min value] default is 0
	-h	Show the usage
	";
}
sub read_expression {
	my ($profile,$log_file) = @_;
	my ($line,@ele);
	my $profile_tag = 0;
	my $total_rep = 0;
	my $total_point = 0;
	my %profile_expression;
	open PROFILE,"<",$profile;
	while($line=<PROFILE>){
		chomp $line;
		if($line =~ /repeat_number=/){
			@ele = split(/=/,$line);
			$total_rep = $ele[1];
			print STDERR "total replication:".$total_rep."\n";
		}elsif($line =~ /timepoint_number=/){
			@ele = split(/=/,$line);
			$total_point = $ele[1];
			print STDERR "total timepoint:".$total_point."\n";
		}
		elsif($line =~ /^rep/){
			@ele = split(/\t/,$line);
			my $tps = scalar @ele - 2;
			if($tps != $total_point){
				print STDERR $line."\n";die;
			}
			my $rep_no = $ele[0];
			my $gene_name = $ele[1];
			my $profile_info;
			$rep_no =~ s/rep//;
			#print STDERR "rep".$rep_no."\n";
			if($rep_no > $total_rep){
				print STDERR "Error, the replication number is larger total replication\n";die;
			}
			for(my $i=2;$i< scalar @ele; $i++){#profile from column2
				$profile_info .= $ele[$i]."\t";
			}
			$profile_info =~ s/\t$//;#remove the last tab
			$profile_expression{$rep_no}{$gene_name} = $profile_info;
		}
	}
	close PROFILE;
	my @all_gene = keys %{$profile_expression{1}};#get the gene name
	my $gene_no = 0;
	open LOG,">",$log_file;
	print LOG "RUN_ID\tGENE_NAME\n";
	foreach my $i (@all_gene){
		#print STDERR $i."\n";
		$gene_index{$gene_no} = $i;
		print LOG $gene_no."\t".$i."\n";
		$gene_no++;
	}
	return ($gene_no,$total_rep,$total_point,%profile_expression);
}

sub read_knowledge {
	my ($total_gene_no,$know_init_file) = @_;
	my ($line,@ele);
	my $format_tag = 0;
	my %knowledge;
	my %fix_status_generation;
	open KNOW,"<",$knowledge;
	while($line=<KNOW>){
		chomp $line;
		#print STDERR $line."\n";
		if($line eq "TF\tGENE\tREGULATORY"){
			print STDERR "format is correct\n";
			$format_tag = 1;
			next;
		}
		if($format_tag){
			@ele = split(/\t/,$line);
			my ($tf,$gene,$regulation) = ($ele[0],$ele[1],$ele[2]);
			$knowledge{$tf}{$gene} = $regulation;
			#print STDERR "TF:".$tf."\tGene:".$gene."\treg:".$regulation."\n";
		}
	}
	close KNOW;
	open INIT,">",$know_init_file;
	for(my $i=0;$i<$total_gene_no;$i++){#each target gene
		my $gene_name = $gene_index{$i};
		my $gene_name_generation = $gene_name."_0";#0 is the initial generation
		my $regulation_knowledge;
		#print STDERR $gene_name.": ";
		for(my $j=0;$j<$total_gene_no;$j++){#each TF
			my $tf_name = $gene_index{$j};
			#print STDERR "TF:".$tf_name."\n";
			if($knowledge{$tf_name}{$gene_name} eq '0'){
				print INIT "0 ";
				$regulation_knowledge .= "0 ";
			}elsif(!$knowledge{$tf_name}{$gene_name}){
				print INIT "? ";
				$regulation_knowledge .= "? ";
			}else{
				print INIT $knowledge{$tf_name}{$gene_name}." ";
				$regulation_knowledge .= $knowledge{$tf_name}{$gene_name}." ";
			}
		}
		$regulation_knowledge =~ s/ $//;
		$fix_status_generation{$gene_name_generation} = $regulation_knowledge;
		#print STDERR "fix_status_generation{".$gene_name_generation."}:".$fix_status_generation{$gene_name_generation}."\n";
		print INIT "\n";
	}
	close INIT;
	return %fix_status_generation;
}
sub generate_config {
	my ($configure,$total_gene_no,$total_rep,$total_point,%profile_expression) = @_;
	open CONF,">",$configure;
	print CONF "num_var ".$total_gene_no."\n";
	print CONF "num_tf ".$total_gene_no."\n";
	print CONF "num_run ".$total_rep."\n";
	print CONF "num_trial ".$total_point."\n";
	print CONF "nochange 1000\n";
	print CONF "delta_t 1\n";
	print CONF "b_up 1\n";
	print CONF "b_down 0\n";
	print CONF "transMax_up 8\n";
	print CONF "transMax_down 0\n";
	print CONF "n_up 4\n";
	print CONF "n_down -4\n";
	print CONF "k_up 4\n";
	print CONF "k_down 0\n";
	print CONF "degrade_up 1\n";
	print CONF "degrade_down 0\n";
	print CONF "list_data\n";
	for(my $i=1;$i<=$total_rep;$i++){
		print CONF $i."\n";
		for(my $j=0;$j<$total_gene_no;$j++){
			my $gene_name = $gene_index{$j};
			if(!$profile_expression{$i}{$gene_name}){
			    print STDERR "i:".$i."\tgene_name:".$gene_name."\n";die;
			}
			print CONF $profile_expression{$i}{$gene_name}."\n";
		}
		#print CONF "\n";
	}

	close CONF;
}
sub run_quantification {
	my ($gen,$use_know,$use_config,$total_gene_no,%confidence_level) = @_;
	my $n_start;
	my $n_end;
	my @threads;
	my $thr_count = 0;
	my $gene_name;
	if($model eq "HFODE"){
		$n_start = 2;
		$n_end = 2+$total_gene_no;
	}
	my $final_results = $output_dir."/final_results.txt";
	open FINAL,">",$final_results;
	print FINAL "TF\tGENE\tREGULATORY\tCONFIDENCE_LEVEL\n";
	for(my $i=0; $i<$total_gene_no; $i++){
		$threads[$thr_count] = threads->new(\&run_iga, $i,$gen,$use_know,$use_config);
		sleep(1) while(scalar threads->list(threads::running) >= $threads);
		$thr_count++;
	}
	foreach my $i (@threads){
		$i->join();
	}
	for(my $i=0; $i<$total_gene_no; $i++){
		$gene_name = $gene_index{$i};
		my $iga_results = $output_dir."/".$gene_name.".txt";
		if(!-e $iga_results){
			print STDERR $gene_name."'s iga result does not exist\n";die;
		}
		my @select_iga_results = select_results(1,$iga_results);#select top 1
		print STDERR "Gene:".$gene_name." ".$select_iga_results[0]."\n";
		my @iga = split(/\t/,$select_iga_results[0]);
		for(my $r=$n_start; $r<$n_end; $r++){
			my $tf_no = $r-2;
			my $tf_name = $gene_index{$tf_no};
			my $role;
			if($iga[$r]>0){
				$role = "+";
			}elsif($iga[$r]<0){
				$role = "-";
			}else{
				$role = 0;
			}
			my $relationship = $tf_name."_".$gene_name;
			print FINAL $tf_name."\t".$gene_name."\t".$role."\t".$confidence_level{$relationship}."\n";
		}
	}
	close FINAL;
}
sub run_EMA {
	my ($gen,$know_filename,$use_config,$total_gene_no,$fix_ref,$confidence_ref,$fix_status_ref,$fix_know_ref) = @_;
	my @fix_value = @{$fix_ref};
	my %confidence_level = %{$confidence_ref};
	my %fix_status_generation = %{$fix_status_ref};
	my %fix_know = %{$fix_know_ref};
	my @threads;
	my $thr_count = 0;
	my $gene_name;
	my $select_top_no = 5;
	my $cutoff = 0.8;
	my $use_know = $tmp_dir."/".$know_filename."_ForStep".$gen.".txt";
	my $next_generation = $gen + 1;
	my $new_know = $tmp_dir."/".$know_filename."_ForStep".$next_generation.".txt";
	
	print STDERR "Step2: GRN decomposition\n";
	print STDERR "Step3: Parallel solving\n";
	for(my $i=0; $i<$total_gene_no; $i++){
		if($fix_value[$i]){
			$gene_name = $gene_index{$i};
			print STDERR $gene_name ." is finished\n";
		}else{
			$threads[$thr_count] = threads->new(\&run_iga, $i,$gen,$use_know,$use_config);
			sleep(1) while(scalar threads->list(threads::running) >= $threads);
			$thr_count++;
		}
	}
	foreach my $i (@threads){
		$i->join();
	}
	print STDERR "Step4: Regulation determination\n";
	open KNOW,">",$new_know;
	for(my $i=0; $i<$total_gene_no; $i++){
		$gene_name = $gene_index{$i};
		my (@regulatory_p,@regulatory_n,@regulatory_z);
		my $iga_results = $output_dir."/".$gene_name.".txt";
		if(!-e $iga_results){
			print STDERR $gene_name."'s iga result does not exist\n";die;
		}
		my @select_iga_results = select_results($select_top_no,$iga_results);#select top 5
		for(my $j=0; $j<$total_gene_no; $j++){
			$regulatory_p[$j] = 0;
			$regulatory_n[$j] = 0;
			$regulatory_z[$j] = 0;
		}
		if($model eq "HFODE"){
			#b(i),TransMax(i),n(i,j),k(i,j),deg(i)
			my $n_start = 2;
			my $n_end = 2+$total_gene_no;
			for(my $j=0; $j<$select_top_no; $j++){
				my @iga = split(/\t/,$select_iga_results[$j]);
				for(my $r=$n_start; $r<$n_end; $r++){
					my $gene_no = $r-2;
					if($iga[$r] > 0){
						if(abs($iga[$r]) > $min){
							$regulatory_p[$gene_no]++;
						}else{
							$regulatory_z[$gene_no]++;
						}
					}elsif($iga[$r] < 0){
						if(abs($iga[$r]) > $min){
							$regulatory_n[$gene_no]++;
						}else{
							$regulatory_z[$gene_no]++;
						}
			
					}else{
						$regulatory_z[$gene_no]++;
					}
				}
			}
		}else{
			#S-system for real GRN (Hij need to initial zero, only use Gij to determine the gene regulation)
			print STDERR "S-system model for real GRN under construction\n";die;
			my $n_start = 1;#Location of Gij 
			my $n_end = 1+$total_gene_no;#Location of Gij
			for(my $j=0; $j<$select_top_no; $j++){
				my @iga = split(/\t/,$select_iga_results[$j]);
				for(my $r=$n_start; $r<$n_end; $r++){
					my $gene_no = $r-1;
					if($iga[$r] > 0){
						$regulatory_p[$gene_no]++;
					}elsif($iga[$r] < 0){
						$regulatory_n[$gene_no]++;
					}else{
						$regulatory_z[$gene_no]++;
					}
				}
			}
		}#S-system model
		my $know_info;
		if(!$fix_value[$i]){#only check unfixed gene
			#print STDERR "fix[".$i."]:".$gene_name." is not fix\n";
			for(my $j=0; $j<$total_gene_no; $j++){
				my $tf_name = $gene_index{$j};
				my $each_gene_knowledge;
				my $relationship = $tf_name."_".$gene_name;
				if($regulatory_p[$j]){
					if(($regulatory_p[$j]/$select_top_no) >= $cutoff){
						$each_gene_knowledge = "+ ";
						if(!$confidence_level{$relationship}){
							$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
						}
					}
				}
				if($regulatory_n[$j]){
					if(($regulatory_n[$j]/$select_top_no) >= $cutoff){
						$each_gene_knowledge = "- ";
						if(!$confidence_level{$relationship}){
							$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
						}
					}
				}
				if($regulatory_z[$j]){
					if(($regulatory_z[$j]/$select_top_no) >= $cutoff){
						$each_gene_knowledge = "0 ";
						if(!$confidence_level{$relationship}){
							$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
						}
					}
				}
				if(!$each_gene_knowledge){
					$each_gene_knowledge = "? ";
				}
				#print STDERR "TF-gene:".$tf_name."-".$gene_name."\n";
				#ßprint STDERR "regulation:".$each_gene_knowledge."\n";
				#print STDERR "P:".$regulatory_p[$j]."\n";
				#print STDERR "N:".$regulatory_n[$j]."\n";
				#print STDERR "Z:".$regulatory_z[$j]."\n";
				#print STDERR "=============================\n";
				$know_info .= $each_gene_knowledge;
			}
			$know_info =~ s/ $//;
		}else{
			#print STDERR "fix[".$i."]:".$gene_name." is fix\n";
			#print STDERR "fix_knowledge{".$i."}:".$fix_know{$i}."\n";
			$know_info = $fix_know{$i};
		}

		
		if(!$fix_value[$i]){#only check unfixed gene
			my $gene_generation_key = $gene_name."_".$gen;
			$fix_status_generation{$gene_generation_key} = $know_info;
			my $last_generation = $gen - 1;
			my $last_gene_generation_key = $gene_name."_".$last_generation;
			print STDERR "Check the status of the determined regulations\n";
			print STDERR "fix_status_generation{".$gene_generation_key."}:".$fix_status_generation{$gene_generation_key}."\n";
			print STDERR "fix_status_generation{".$last_gene_generation_key."}:".$fix_status_generation{$last_gene_generation_key}."\n";
			if($fix_status_generation{$gene_generation_key} eq $fix_status_generation{$last_gene_generation_key}){
				print STDERR "Get the same fix regulations, so need to fix more regulations\n";
				my @select_iga_results = select_results(1,$iga_results);#select top 1
				my @iga = split(/\t/,$select_iga_results[0]);
				my $gene_knowledge;
				if($model eq "HFODE"){
					my $n_start = 2;
					my $n_end = 2+$total_gene_no;
					my $role;
					for(my $r=$n_start; $r<$n_end; $r++){
						my $tf_index = $r-2;
						my $tf_name = $gene_index{$tf_index};
						my $relationship = $tf_name."_".$gene_name;
						if($iga[$r] > 0){
							$role = "+";
							if(!$confidence_level{$relationship}){
								$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
							}
						}elsif($iga[$r] < 0){
							$role = "-";
							if(!$confidence_level{$relationship}){
								$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
							}
						}else{
							$role = 0;
							if(!$confidence_level{$relationship}){
								$confidence_level{$relationship} = sprintf("%.3f",1/$gen);
							}
						}
						$gene_knowledge .= $role." ";
					}
				}else{
					#model == "S-system"
					#S-system for real GRN (Hij need to initial zero, only use Gij to determine the gene regulation)
					print STDERR "S-system model for real GRN under construction\n";die;
					my $n_start = 1;
					my $n_end = 1+$total_gene_no;
				}
				$gene_knowledge =~ s/ $//;
				$know_info = $gene_knowledge;
				$fix_status_generation{$gene_generation_key} = $know_info;
				print STDERR "Update fix_status_generation{".$gene_generation_key."}:".$fix_status_generation{$gene_generation_key}."\n";
			}
		}
		if(!$know_info){
			print STDERR "no value, error\n";die;
		}	
		print KNOW $know_info."\n";	
	}
	return (\%confidence_level, \%fix_status_generation);
}
sub run_iga {
	my ($gene_no,$gen,$know,$conf) = @_;
	my $src_dir = dirname($program);
	if($model eq "HFODE"){
		my $ema_HFODE = $src_dir."/EMA_HFODE/EMA_HFODE";
		if(!-e $ema_HFODE){
			print STDERR $ema_HFODE." does not exist\n";die;
		}
		if(!-e $conf){
			print STDERR $conf." does not exist\n";die;
		}
		my $gene_name = $gene_index{$gene_no};
		my $command = $ema_HFODE." -i ".$gene_no." -n 30 -m 0 -G ".$generation_no." -I ".$gen." -F ".$fitness_type." -c ".$cc." ".$conf." ".$know." > ".$output_dir."/".$gene_name.".txt";
		#print STDERR $command."\n";
		`$command`;
	}

}
sub select_results {
	my ($top_no,$results) = @_;
	my (@unsort_fitness,@sort_fitness,%fitness_value);
	my @select_results;
	my ($line,@ele);
	my $no = 0;
	open IGA,"<",$results;
	while($line=<IGA>){
		chomp $line;
		@ele = split(/ /,$line);
		my $all_parameter_no = scalar @ele;
		my $iga_info;
		for(my $j=0; $j< ($all_parameter_no - 2); $j++){
			$iga_info .= $ele[$j]."\t";
		}
		$iga_info =~ s/\t$//;
		my $fitness_value = $ele[$all_parameter_no - 2];
		#print STDERR "fitness:".$fitness_value."\n";
		$unsort_fitness[$no] = $fitness_value;
		$fitness_value{$fitness_value} = $iga_info;
		$no++;
	}
	close IGA;
	@sort_fitness = sort {$a <=> $b} @unsort_fitness;
	for(my $i=0; $i<$top_no; $i++){
		my $fitness_value = $sort_fitness[$i];
		#print STDERR $i." fitness:".$fitness_value."\t".$fitness_value{$fitness_value}."\n";
		$select_results[$i] = $fitness_value{$fitness_value};
		#print STDERR "select result:".$select_results[$i]."\n";
	}
	return @select_results;
}
sub check_knowledge {
	my ($know,$gene_no, $fix_know_ref) = @_;
	my %fix_know = %{$fix_know_ref};
	my ($line,@ele);
	my $no = 0;
	open (KNOW,"<",$know) or die "open file ".$know." error\n";
	while($line=<KNOW>){
		chomp $line;
		if($no == $gene_no){
			#print STDERR $line."\n";
			my $unfix_no = 0;
			@ele = split(/ /,$line);
			foreach my $i(@ele){
				if($i eq "?"){
					$unfix_no++;
				}
			}
			if($unfix_no == 0){
				$fix_know{$no} = $line;
				return (1, \%fix_know);
			}else{
				$fix_know{$no} = "unfix";
				return (0, \%fix_know);
			}
		}
		$no++;
	}
	close KNOW;
}
sub main {
	my $command;
	my $generation = 1;#first generation;
	my $log = $output_dir."/run.log";
	my $knowledge_filename_tmp = basename($knowledge);
	my $knowledge_filename;
	my %fix_status;
	my $fix_ref;
	my %fix_knowledge;
	my $fix_knowledge_ref;
	my $config = $expression."_config";
	my @fix;
	my $all_fix = 0;
	my $total_fix_no = 0;	
	my %confidence;
	my $confidence_ref;
	if($knowledge_filename_tmp =~ /\./){
		my @filename = split(/\./,$knowledge_filename_tmp);
		$knowledge_filename = $filename[0];
	}else{
		$knowledge_filename = $knowledge_filename_tmp;
	}
	my $use_knowledge = $tmp_dir."/".$knowledge_filename."_ForStep".$generation.".txt";
	print STDERR "know_init:".$use_knowledge."\n";
	my ($total_gene,$total_repeat,$total_data_points,%profile) = read_expression($expression,$log);
	print STDERR "Number of gene:".$total_gene."\n";
	%fix_status = read_knowledge($total_gene,$use_knowledge);#generate the initial knowledge
	if(!-e $use_knowledge){
		print STDERR "generate initial knowledge failed <".$use_knowledge.">\n";die;
	}
	generate_config($config,$total_gene,$total_repeat,$total_data_points,%profile);
	
	open LOG,">",$log;
	print LOG "Total gene:".$total_gene."\nTotal repeat:".$total_repeat."\nTotal data points:".$total_data_points."\n";
	print LOG "Use know init:".$use_knowledge."\nUse config:".$config."\nUse input data:".$expression."\n";
	print LOG "Use Fitness".$fitness_type."\nUse threads:".$threads."\nUse generation:".$generation_no."\nUse CC:".$cc."\nUse min:".$min."\n";
	close LOG;

	print STDERR "Step1: Initialisation\n";	
	for(my $i=0;$i<$total_gene;$i++){
		($fix[$i], $fix_knowledge_ref) = check_knowledge($use_knowledge,$i,\%fix_knowledge);
		%fix_knowledge = %{$fix_knowledge_ref};
		$total_fix_no += $fix[$i];
	}
	print STDERR "Initial Generation ".$generation." fix: ".$total_fix_no."\n";
	if($total_fix_no == $total_gene){
		$all_fix = 1;
	}
	
	#=======finished initial=================================================
	if(!$all_fix){
		do{
			print STDERR "Generation ".$generation." start running\n";
			if(!-e $use_knowledge){
				print STDERR "use knowledge does not exist <".$use_knowledge.">\n";die;
			}
			($confidence_ref,$fix_ref) = run_EMA($generation,$knowledge_filename,$config,$total_gene,\@fix,\%confidence,\%fix_status,\%fix_knowledge);
			%confidence = %{$confidence_ref};
			%fix_status = %{$fix_ref};
			$generation++;
			my $new_knowledge_file = $tmp_dir."/".$knowledge_filename."_ForStep".$generation.".txt";
			$total_fix_no = 0;
			print STDERR "Step5: GRN combination\n";
			for(my $i=0;$i<$total_gene;$i++){
				($fix[$i],$fix_knowledge_ref) = check_knowledge($new_knowledge_file,$i,\%fix_knowledge);
				%fix_knowledge = %{$fix_knowledge_ref};
				#print STDERR "fix[".$i."]:".$fix[$i]."\n";
				#print STDERR "fix_knowledge{".$i."}:".$fix_knowledge{$i}."\n";
				$total_fix_no += $fix[$i];
			}
			my $last_gen = $generation - 1;
			print STDERR "Generation ".$last_gen." finish, fix: ".$total_fix_no."\n";
			print STDERR "Step6: Qualitative termination test\n";
			if($total_fix_no == $total_gene){
				$all_fix = 1;
				print STDERR "All regulations are determined, go to Step8\n";
			}else{
				print STDERR "Step7: Inheritance\n";
			}
		}while(!$all_fix)
	}
	my $final_knowledge = $tmp_dir."/".$knowledge_filename."_ForStep".$generation.".txt";
	if(!-e $final_knowledge){
		print STDERR "use final knowledge does not exist <".$final_knowledge.">\n";die;
	}
	print STDERR "Step8:GRN quantification\n";
	run_quantification($generation,$final_knowledge,$config,$total_gene,%confidence);
	
	
	$command = "mv calprofile* ".$tmp_dir;
	print STDERR $command."\n";
	`$command`;
	$command = "mv progress* ".$tmp_dir;
	print STDERR $command."\n";
	`$command`;
	$command = "mv printdata.txt ".$tmp_dir;
	print STDERR $command."\n";
	`$command`;
}
#=========subroutine end================================
GetOptions(
	'i=s'	=>\$expression,
	'o=s'	=>\$output_dir,
	'k=s'	=>\$knowledge,
	't=i'	=>\$threads,
	'm=s'	=>\$model,
	'g=i'	=>\$generation_no,
	'f=i'	=>\$fitness_type,
	'c=i'	=>\$cc,
	'min=s'	=>\$min,
	'h'	=>\$help,
); 
if(!$expression ||!$knowledge || !$threads || $help ||$threads < 1 || !$output_dir){
	Usage();die;
}
if($min < 0){
	print STDERR "min value must be >= 0\n";
	Usage();die;
}
if(($cc != 0)and($cc != 1)){
	print STDERR "type of cc is error ".$cc."\n";
	Usage();die;
}
if(!-e $expression){
	print STDERR $expression." does not exist\n";
	Usage();die;
}
if(!-e $knowledge){
	print STDERR $knowledge." does not exist\n";
	Usage();die;
}
if(($model ne "s-system")and($model ne "HFODE")){
	print STDERR "type of model is error\n";
	Usage();die;
}
if($generation_no < 100){
	print STDERR "generation must be > 100\n";
	Usage();die;
}
if(($fitness_type != 0)and($fitness_type != 1)and($fitness_type != 2)and($fitness_type != 3)and($fitness_type != 4)){
	print STDERR "type of fitness function is error\n";
	Usage();die;
}
if(-d $output_dir){
	$output_dir =~ s/\/$//;
	my $final = $output_dir."/final_results.txt";
	if(-e $final){
		print STDERR "This output folder already has the final_results.txt\n";die;
	}
	$tmp_dir = $output_dir."/tmp";
	if(!-d $tmp_dir){
		my $cmd = "mkdir ".$tmp_dir;
		print STDERR "create tmp folder:".$cmd."\n";
		`$cmd`;
	}
}else{
	print STDERR "Directory of output does not exist\n";
	die;
}
main();

