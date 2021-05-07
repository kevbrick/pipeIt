#!/usr/bin/perl
use strict;
use Cwd qw /abs_path/;
use List::Util qw/max min/;
use Getopt::Long;
use Time::localtime;
use Math::Round;

my $opts = $#ARGV;

GetOptions ('bam=s'        => \(my $bam),
            'f1=s'         => \(my $fqR1),
            'f2=s'         => \(my $fqR2),
            'g=s'          => \(my $genome),
            'obj=s'        => \(my $objectStoreRegex),
            'sra=s'        => \(my $SRAIDs),
            's=s'          => \(my $sample),
            'o=s'          => \(my $outName),
            'od=s'         => \(my $outDir),
            'wd=s'         => \(my $workingDir),
            'd=s'          => \(my $rundate),
            'long+'        => \(my $runLonger),
            'n=s'          => \(my $threads=16),
            'm=s'          => \(my $mem=32),
            'time=s'       => \(my $nTime),
            'lR1=s'        => \(my $read1Length="36"),
            'lR2=s'        => \(my $read2Length="40"),
            'type=s'       => \(my $PEorSR),
            'hicbins=s'    => \(my $hic_binsizes = "1000000,750000,500000,250000,100000,75000,50000,25000,10000,7500,5000,2500,1000"),
            'ssz=i'        => \(my $fastqSplitSize=10000000),
            'config=s'     => \(my $configFile),
            'dep=s'        => \(my $depJobID),
            'inode+'       => \(my $runOnCurrentNode),
            'hide+'        => \(my $hiddenScripts),
            'l+'           => \(my $runLocal),
            't+'           => \(my $outToTmp),
            'homeOK+'      => \(my $homeOK),
            'dryrun+'      => \(my $dryRun),
            'nobamyet+'    => \(my $missingBAMok),
            'debug+'       => \(my $debugMe),
            'h+'           => \(my $help),
            'H+'           => \(my $extraHelp),
            'plist+'       => \(my $extraHelpAlt),
            'gcCorr=s'     => \(my $gcCorrectionFile),
            'mod+'         => \(my $useModules),
            'test+'        => \(my $testData),
            'profile=s'    => \(my $nxfProfiles),
            'pipe=s'       => \(my $pipeList),
            'hicaligner=s' => \(my $hic_aligner = 'bwa'));

#my $todayDate = getTodayDate()sub{sprintf '%02d%02d%02d',
#    $_[5]+1900-2000, $_[4]+1, $_[3]}->(localtime);
my $todayDate = getTodayDate();

my ($restrictionEnzyme, $hic_ligation_site, $hic_restriction_site, $dnase, $phased);

## KB 7/19/18 - run each NF instance in subfolder
## Get start folder and set output folder full path if required
my $startFolder = abs_path('.').'/';

## Set extrahelp flag if --plist is used instead of --H
$extraHelp++ if ($extraHelpAlt);

## PRINT HELP OR ERRORS
if ($help || $extraHelp || $opts < 1){
	printARGS();
	if ($extraHelp){
		printPipeDetails();
	}
	exit();
}

my $errMSG;
$errMSG .= 'pipeIt cannot be run from your /home folder. Please use either /data/{user} or /data/RDCO ...'."\n" if (!$homeOK && (`pwd` =~ /(home|\~)/ && `pwd` =~ /($ENV{USER})/));
$errMSG .= 'Too many input options specified ... pick one. '."\n" if ((defined($bam) + defined($fqR1) + defined($objectStoreRegex) + defined($SRAIDs)) > 1);
$errMSG .= 'No BAM or FASTQ file(s) provided.'."\n" unless ($bam || $fqR1 || $objectStoreRegex || $SRAIDs || $missingBAMok);
$errMSG .= "BAM file does not exist [--bam $bam]"."\n"    if (not($missingBAMok) && $bam && not (-e $bam));
$errMSG .= "SRA names are incorrect (must begin with SRR / DRR / ERR) [--sra $SRAIDs]"."\n"      if ($SRAIDs && !(sraOK($SRAIDs)));
$errMSG .= "Config file specified but does not exist [--config $configFile]"."\n"    if ($configFile && not (-e $configFile));
$errMSG .= "FASTQ file does not exist [--f1 $fqR1]"."\n" if ($fqR1 && not (-e $fqR1));
$errMSG .= "FASTQ file does not exist [--f2 $fqR2]"."\n" if ($fqR2 && not (-e $fqR2));
$errMSG .= 'Env. var DSL2DIR required'."\n" unless ($ENV{DSL2DIR});
$errMSG .= 'Env. var NXF_GENOMES required (suggest: /data/RDCO/genomes)'."\n" unless ($ENV{NXF_GENOMES});
$errMSG .= 'Env. var NXF_ANNOTATION required (suggest: /data/RDCO/annotation)'."\n" unless ($ENV{NXF_ANNOTATION});
$errMSG .= 'samtools NOT found at required location (DSL2DIR = '.$ENV{DSL2DIR}.')'."\n" unless (-e $ENV{DSL2DIR}.'/samtools');
#$errMSG .= 'Restriction enzyme (--re) required for HiC pipeline'."\n" if ($pipeList =~ /hic/i && !($restrictionEnzyme));

if ($pipeList =~ /rtseq/i){
	if ($gcCorrectionFile && not (-e $gcCorrectionFile)){
		$gcCorrectionFile = $ENV{DSL2DIR}.'/accessoryFiles/rtSeq/gcCalibration/'.$gcCorrectionFile.'.tab';
		if (not (-e $gcCorrectionFile)){
			$errMSG .= 'RTseq pipe: GC correction stat file does not exist (--gcCorr) ['.$gcCorrectionFile."]\n" ;
			$errMSG .= "Options are: \n" ;
			$errMSG .= "Admera_HiSeqX \n" ;
			$errMSG .= "NIDDK_HiSeq2500 \n" ;
			$errMSG .= "Nussenzweig_HiSeq2000 \n" ;
			$errMSG .= "Yehuda_NextSeq500 \n" ;
		}
	}
	$gcCorrectionFile = abs_path($gcCorrectionFile) if ($gcCorrectionFile);
}

if ($genome eq 'guess'){
	my $inFileName;
	$inFileName = $bam  if ($bam);
	$inFileName = $fqR1 if ($fqR1);

	$genome = guessGenome($inFileName);
	$errMSG .= "Genome GUESS failed .... [$inFileName]\nGenome is required (--g)\n" unless ($genome);
}else{
	#$genome = guessGenome(".".$genome.".");
	$errMSG .= 'Genome is required (--g)'."\n" unless ($genome);
}

## KB 181128
$errMSG .= 'Genome NOT found (--g $genome)'."\n" unless (-d ($ENV{NXF_GENOMES}.'/'.$genome));

## KB 190108: Allow multiple pipes
my $initialOutDir = $outDir;
for my $pipeType(split(/,/,$pipeList)){

  ## Standardize inputs
	$pipeType = lc($pipeType);
  $pipeType = 'mm2' if ($pipeType =~ /^(minimap2|minimap|mm2)$/);

	## Assure we don't keep renaming convention from a previous pipe
	$outDir = $initialOutDir;

	if ($pipeType =~ /^(bwa|bwaaln|mm2|bowtie|ssds|ssdsorig)$/){
		$errMSG .= 'Invalid split Size (--ssz). Please provide a number ... '."\n"    unless ($fastqSplitSize =~ /^\d+$/);
		print STDERR 'Split size is very small (--ssz '.$fastqSplitSize.'). Please reconsider ... this will spawn LOTS of jobs if the BAM is large... '."\n" if ($fastqSplitSize < 5000000);
	}

	if ($SRAIDs){
		if ($pipeType !~ /(bwa|mm2|ssds|ssdsorig|minimap|ont)/){
			$errMSG .= 'SRA accessions are not yet supported for the '.$pipeType.' pipeline ... '."\n";
		}
	}

	if ($objectStoreRegex){
		if ($pipeType !~ /(ssds|ssdsorig|bwa|mm2|ont|minimap|rtseq|bowtie|hic)/){
			$errMSG .= 'Alignment from object store is not yet supported for the '.$pipeType.' pipeline ... '."\n";
		}
	}

	if ($pipeType =~ /hic/){
		if ($pipeType =~ s/(hic):([^\:]+):([^\:]+)/$1/){
			($restrictionEnzyme,$phased) = (uc($2),($3 eq "phased")?1:0);
			if ($phased){
				print STDERR "## ALERT: PHASED Hi-C pipeline is ONLY configured for B6 x CAST \n";
				if ($genome !~ /mm10_nCAST/){
					$errMSG .= "Phased HiC is only configured to be used with --genome mm10_nCAST [NOT $genome]. Please modify --genome argument.\n"
				}
			}
		}

		if (uc($restrictionEnzyme) =~ /^MBO[I1]$/){
			$hic_ligation_site = 'GATCGATC';
			$hic_restriction_site = '^GATC';
			print STDERR '## NOTE: Hi-C pipeline using '.uc($restrictionEnzyme)."\n";
		}

		if (uc($restrictionEnzyme) =~ /^DNASE$/){
			$dnase++;
			$hic_ligation_site = '';
			$hic_restriction_site = '';
			print STDERR "## NOTE: Hi-C pipeline using DNAse (micro-C)\n";
		}

		if ($pipeType =~ s/:singlecell//i){
			$nxfProfiles .= ',singlecell'
		}

		if ($pipeType =~ s/:bowtie[^:]*//i){
			$hic_aligner = 'bowtie';
		}

		$errMSG .= 'pipeit is not configured to run the hi-c pipeline with that enzyme ['.uc($restrictionEnzyme).'] ... Please talk to Kevin ...'."\n" unless ($dnase || $hic_ligation_site);
	}

	if ($errMSG){
		printARGS($errMSG);
		exit();
	}

	## SET isNFCORE FLag for nf-core pipelines
  my $isNFCORE = 0;

	if ($pipeType =~ /^(methyl)$/){
		print STDERR '## NOTE: Using an nf-core pipeline'."\n";
		$isNFCORE = 1;
	}

	## KB 7/19/18 - replace relative with absoulte paths (Just to be safe)
	$bam         = abs_path($bam)         if ($bam);
	$fqR1        = abs_path($fqR1)        if ($fqR1);
	$fqR2        = abs_path($fqR2)        if ($fqR2);

	## SET config file
	if ($configFile){
		$configFile  = abs_path($configFile)  ;
	}else{
		if ($isNFCORE){
			$configFile = abs_path('/data/RDCO/code/nf-core/config/biowulf.config.nf') ;
		}else{
			$configFile = abs_path($ENV{'DSL2DIR'}.'/config/nextflow.config.nf');
		}
	}

	## SET NEXTFLOW PROFILES
	## USE profiles to import configs for each pipe
	my $profiles;
	$profiles .= ",modules"     if ($useModules);
	$profiles .= ",singularity" if (!$useModules || $isNFCORE);
	$profiles .= ",local"       if ($runLocal);
  $profiles .= ",slurm"       if (!$runLocal);

	$profiles   .= ",ssds"     if ($pipeType =~ /ssds/);
	$profiles   .= ",align"    if ($pipeType =~ /(bwa|bwaaln|mm2|bowtie)/);
	$profiles   .= ",rtseq"    if ($pipeType =~ /^rtseq$/);
	$profiles   .= ",hic4dn"   if ($pipeType =~ /(hic)/);
	$profiles   .= ",commitFQ" if ($pipeType =~ /(commit)/);

	## Allow custom profiles (this can break things too though!)
	$profiles .= ",$nxfProfiles" if ($nxfProfiles);

	$profiles =~ s/^,//; $profiles =~ s/,,/,/g;

	## Check for dependencies (if required)
	my $dependency;
	$dependency = " --dependency=afterany:$depJobID " if ($depJobID);

	#-----------------------------------------------------------------------
	## Get output name if undefined
	unless ($outName){
		if ($SRAIDs){
			$outName = $SRAIDs;
		}else{
			if ($objectStoreRegex){
				$outName = $objectStoreRegex;
			}else{
				$outName = $bam?$bam:$fqR1;
				$outName =~ s/^.*\///;
				$outName =~ s/([\._]R1[\._]|[\._]R2[\._])/./g;
				$outName =~ s/(\.gz|\.fastq|\.fq|\.bam)/./g;
				$outName =~ s/(\.bwaMemPE|\.bwaMemSR|\.MM2|\.MM2PE|\.MM2SR|\.SSDS|\.unaligned)/\./g;
			}
		}

		for my $genomeName(`ls -d $ENV{NXF_GENOMES}\/\*\/`){
			chomp $genomeName;
			$genomeName =~ s/\/$//g;
			$genomeName =~ s/^.+\///g;
			$genomeName =~ s/[\.\/]//g;
			$outName =~ s/\.($genomeName)[\.]/\./g;
		}

		$outName =~ s/^(.+?)\.+$/$1/;
	}

	$outName =~ s/^\s+//;

	## Set other names if undefined
	my $rundate = $todayDate  unless ($rundate);

	## Get absolute path
	if (lc($pipeType) eq 'commitfq'){
		print STDERR "commitFQ uses the default output folder (/data/RDCO)\n";
		$outDir = '/data/RDCO/fastq/';
	}else{
		if ($outDir){
			if ($outDir =~ /\//){
				$outDir  = $outDir;
			}else{
				$outDir  = $startFolder.$outDir ;
			}
		}else{
			$outDir  = $startFolder.$outName ;
		}
		## Allow for pipeline params passed as colons (Hi-C pipeline req'd)
		my $outdir_pipetype = $pipeType;
		$outdir_pipetype =~ s/^(.+?)\:\S+$/$1/;
		$outDir .= "_$outdir_pipetype";
	}

	## assure path is good !!
	$outDir = abs_path($outDir);

	if (!$homeOK && ($outDir =~ /(home|\~)/ && $outDir =~ /($ENV{USER})/)){
		$errMSG .= 'pipeIt cannot use your /home folder for output. Please use either /data/{user} or /data/RDCO ...'."\n" ;
		printARGS($errMSG);
		exit();
	}

	system('mkdir -p '.$outDir);

	## Define script to run pipeline
	#my $randName         = $pipeType.'_nxf_'.int(rand()*10000000000);
	my $randName              = $pipeType.'_'.$outName.'_'.timestamp();
	my $runFolder             = $outDir.'/'.($workingDir?$workingDir:$randName).'/';
	my $outScript             = $runFolder.($hiddenScripts?'.':'').$randName.'_pipeIt.sh';
	my $outSampleSheet        = $runFolder.($hiddenScripts?'.':'').'SampleSheet_'.$randName.'.tab';
	my $outResumeScript       = $runFolder.($hiddenScripts?'.':'').'RESUME_'.$randName.'_pipeIt.sh';
	my $outResumeScriptLocal  = $runFolder.($hiddenScripts?'.':'').'local_RESUME_'.$randName.'_pipeIt.sh';

	#my $workFolder = $runFolder.'work_'.$outName;
	my $workFolder = $runFolder.'work';

	## KB Sept 21 2017: allow custom config files
	my $nextFlowArgs;
	my $myConfig = $configFile;

	if ($pipeType =~ /^hic/){
		$myConfig .= ",".abs_path($ENV{'DSL2DIR'}.'/pipelines/config/config_4dn.nf')
	}

	$nextFlowArgs .= " -c $myConfig -profile $profiles -work-dir $workFolder ";

	## KB Nov 21 2018: Set runtimes
	my ($sBatchTime,$nxfTime);
	if ($nTime){
		($sBatchTime,$nxfTime) = getTime($nTime);
	}

  ## get pipeline and ARGS
	my ($nextFlowPipe, $pipeArgs) = whatPipe($bam,$fqR1,$fqR2,$objectStoreRegex,$SRAIDs,$pipeType,$PEorSR,$runFolder);

	system('mkdir '.$runFolder);

  ## 030821 KB:
	## For nf-core pipelines, we want to have the fastqs in the run folder
	## This sidesteps issues with parsing the input string (see whatPipe function)
	if ($isNFCORE){
		system("ln -s $fqR1 $runFolder/R1.fastq.gz");
		system("ln -s $fqR2 $runFolder/R2.fastq.gz");
	}

	open SC, '>'.$outScript;

	print SC '#!/bin/bash'."\n\n";

	## Modified Sept 25 2020: Kevin Brick
	#print SC 'module load nextflow/20.01.0'."\n"; ## April 2020
	#print SC 'module load singularity/3.6.4'."\n";
  print SC 'module load nextflow/20.10.0'."\n";
	print SC 'module load singularity'."\n";
  ## END MOD Sept 25

  print SC 'module load java/1.8.0_92'."\n";
	print SC 'module load picard/2.9.2'."\n";
	print SC 'module load graphviz/2.40'."\n";
	print SC 'export TMPDIR=/lscratch/$SLURM_JOBID'."\n";

	## NOTE : this is just a placeholder so that nextflow does not die.
	##        If PYTHONPATH is empty, nextflow throws an error is the
	##        SSDS pipeline when running multiQC
	print SC 'export PYTHONPATH=/data/RDCO'."\n";

	print SC "cd $runFolder\n";

	print SC "nextflow run $nextFlowArgs $nextFlowPipe $pipeArgs ||exit -1\n";

	print SC "rm -rf $workFolder \n";
  print SC "rm $outResumeScript \n";
	print SC "rm $outResumeScriptLocal \n";
	print SC "touch pipelineComplete.OK\n";

	close SC;

	## Make resume script :)
	open RESUMESCRIPT,'>',$outResumeScript;
	open RESUMESCRIPTLOCAL,'>'.$outResumeScriptLocal;

	open IN, $outScript;
	while (<IN>){
		if ($_ =~ s/(nextflow\srun)/$1 -resume/){
			print RESUMESCRIPT      "cleanNXFwork $runFolder\n";
			print RESUMESCRIPTLOCAL "cleanNXFwork $runFolder\n";
			$_ =~ s/(\-profile\s+\S*)local/$1."slurm"/e;
			print RESUMESCRIPT      $_;

			$_ =~ s/(\-profile\s+\S*)slurm/$1."local"/e;
			print RESUMESCRIPTLOCAL $_;
		}else{
			print RESUMESCRIPT      $_;
			print RESUMESCRIPTLOCAL $_;
		}
	}
	close RESUMESCRIPT;
	close RESUMESCRIPTLOCAL;

	my $fileChk = $bam?$bam:$fqR1;
	my $inputSize = (-s $fileChk)/1000000000;
	my $pipeItJobID;

	unless ($dryRun){

		my $cmdOK; ## Assure only 1 command is run
		my $pipeExec;

		## Runs locally on an interactive node
		if ($runOnCurrentNode){
			#system("bash $outScript");
			$pipeExec = "bash $outScript";
			$cmdOK++;
		}

		if (!$cmdOK && ($inputSize > 60 || $pipeType =~ /hic/ || $runLonger || $isNFCORE)){
			#system("sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J P.".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"32":"32")."g --gres=lscratch:800 ".$dependency." --time=".($sBatchTime?$sBatchTime:"96:03:00")." --partition=norm $outScript");
			$pipeExec = "sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J nx.".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"32":"32")."g --gres=lscratch:800 ".$dependency." --time=".($sBatchTime?$sBatchTime:"96:03:00")." --partition=norm $outScript";
			$cmdOK++;
		}

		if (!$cmdOK && $inputSize > 30){
			#system("sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J P.".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"32":"32")."g --gres=lscratch:800 ".$dependency." --time=".($sBatchTime?$sBatchTime:"72:03:00")." --partition=norm $outScript");
			$pipeExec = "sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J nx".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"32":"32")."g --gres=lscratch:800 ".$dependency." --time=".($sBatchTime?$sBatchTime:"72:03:00")." --partition=norm $outScript";
			$cmdOK++;
		}

		if (!$cmdOK && $inputSize > 10){
			#system("sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J P.".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"32":"32")."g --gres=lscratch:800 ".$dependency." --time=".($sBatchTime?$sBatchTime:"48:03:00")." --partition=norm $outScript");
			$pipeExec = "sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J nx.".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"32":"32")."g --gres=lscratch:800 ".$dependency." --time=".($sBatchTime?$sBatchTime:"48:03:00")." --partition=norm $outScript";
			$cmdOK++;
		}

		if (!$cmdOK && $inputSize < 1){
			#system("sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J P.".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"16":"32")."g --gres=lscratch:800 ".$dependency." --time=".($sBatchTime?$sBatchTime:"24:03:00")." --partition=norm $outScript");
			$pipeExec = "sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J nx.".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"16":"32")."g --gres=lscratch:800 ".$dependency." --time=".($sBatchTime?$sBatchTime:"24:03:00")." --partition=norm $outScript";
			$cmdOK++;
		}

		if (!$cmdOK){
			$pipeExec = "sbatch --mail-type=BEGIN,TIME_LIMIT_90,END -J nx.".$randName." --ntasks=".($runLocal?"16":"1")." --mem=".($runLocal?"32":"32")."g --gres=lscratch:200 ".$dependency." --time=".($sBatchTime?$sBatchTime:"24:03:00")." --partition=norm $outScript";
			$cmdOK++;
		}

		print STDERR "$pipeExec \n";

		## Execute job
		$pipeItJobID = `$pipeExec`; chomp $pipeItJobID;

	}

	if ($dryRun){
		print STDERR "$nextFlowPipe pipeline NOT Started : Script created ... \n" ;
		print $outScript."\n";
	}else{
		print STDERR "$nextFlowPipe pipeline Started\n" ;
		print STDERR "pipeIt job ID : $pipeItJobID\n" ;
		print $pipeItJobID."\n";
	}

	print STDERR "-------------------------------\n";
	print STDERR "sbatch script : \n$outScript\n";
	print STDERR "-------------------------------\n";

}

########################################################################
sub whatPipe{
	my ($wBAM,$wFQ1,$wFQ2,$objRegex,$mySRA,$type,$peORsr,$runFolder) = @_;

	## Figure out what pipeline to use
	#  Returns the path to the nextflow pipeline
	#  and the ARGS to use

	my $pipeNF;
	$pipeNF = $ENV{DSL2DIR}.'/pipelines/align.nf'    if ($type =~ /^(bwa|bwaaln|mm2|bowtie)$/);
	$pipeNF = $ENV{DSL2DIR}.'/pipelines/ssds.nf'     if ($type =~ /^(ssds|ssdsorig)$/);
  $pipeNF = $ENV{DSL2DIR}.'/pipelines/commitFQ.nf' if ($type =~ /^(commitfq)$/);
	$pipeNF = $ENV{DSL2DIR}.'/pipelines/rtseq.nf'    if ($type =~ /^rtseq$/);
	$pipeNF = $ENV{DSL2DIR}.'/pipelines/hic4dn.nf'   if ($type =~ /^hic$/);

  ## Allow nf-core pipelines
	$pipeNF = 'nf-core/methylseq -r 1.5 ' if ($type =~ /^methyl$/);
	#$pipeNF = 'nf-core/hic'       if ($type =~ /^microc$/);
	#$pipeNF = 'nf-core/hic'       if ($type =~ /^hic$/);
	#$pipes{uBAM}     = $ENV{DSL2DIR}.'/pipelines/archive.nf';
	#$pipes{bismark}  = $ENV{DSL2DIR}.'/pipelines/bismark.nf';
	#$pipes{HiC}      = $ENV{DSL2DIR}.'/pipelines/hic.nf';
	#$pipes{atacseq}  = $ENV{DSL2DIR}.'/pipelines/atacseq.nf';

	my $pArgs;
	## Allow for NF-CORE PIPELINES
  if ($pipeNF =~ /nf-core/){
		## 030521 KB: nf-core pipelines don't like you to mess with this string
		##            i.e. Sample1_R{1,2}.fastq.gz can cause issues because it doesn't like text before the R{1,2} string !!!
		##            Odd behaviour, but we need to work with it. For nf-core, symlinks to fqs now go to run folder ...
		#my $readsWCString = $wFQ1; $readsWCString =~ s/R1/\R\{1\,2\}/;
		my $readsWCString = 'R{1,2}.fastq.gz';

		#KB 04-28-21: Modified for newer nf-core pipelines
		#$pArgs .= " --input \'$readsWCString\' " ;
		if ($pipeNF =~ /methylseq/){
			$pArgs .= " --reads \'$runFolder\/$readsWCString\' " ;
		}else{
			$pArgs .= " --input \'$runFolder\/$readsWCString\' " ;
		}
		$pArgs .= " --outdir $outDir " ;

		if ($type =~ /methyl/){
			$pArgs .= ' --aligner bwameth ';
			$pArgs .= ' --fasta '.$ENV{NXF_GENOMES}.'/'.$genome.'/genome.fa ' ;
			$pArgs .= ' --fasta_index '.$ENV{NXF_GENOMES}.'/'.$genome.'/genome.fa.fai ' ;
			$pArgs .= ' --bwa_meth_index '.$ENV{NXF_GENOMES}.'/'.$genome.'/BWAMethIndex/genome.fa ' if (-e $ENV{NXF_GENOMES}.'/'.$genome.'/BWAMethIndex/genome.fa.bwameth.c2t.sa');
			$pArgs .= ' --igenomes_ignore --save_reference ';
			$pArgs .= ' --single-end ' if ($wFQ1 && !$wFQ2);
		}

	}else{
		if ($type =~ /^hic$/){
			my $sampleid = $wFQ1;
			$sampleid =~ s/^.+\///;
			$sampleid =~ s/^(.+)[\._]R[12].(fastq|fastq.gz|fq|fq.gz)$/$1/;

			my $reFile = $ENV{NXF_GENOMES}.'/'.$genome.'/re/MboI_sites.txt';

			$pArgs .= ($dnase?" --dnase true ":((-e $reFile)?" --dnase false --refile $reFile":" --dnase false")).($phased?" --phased true ":" --phased false");
			$pArgs .= " --aligner ".$hic_aligner ;

			## NEEDED FOR BALANCING COOL FILE
			my $blacklist = $ENV{NXF_GENOMES}.'/'.$genome.'/blacklist/blacklist.bed';
			$pArgs .= " --blacklist $blacklist" if (- $blacklist);

			if ($fastqSplitSize < 10000000){
				print STDERR "#######################################################################\n";
				print STDERR " ALERT  : Minimum FASTQ split size for HiC pipeline is 25 Million reads\n";
				print STDERR "          Overriding specified size: $fastqSplitSize\n";
				print STDERR "          Split size set to 10000000\n";
				print STDERR "#######################################################################\n";

				$fastqSplitSize = 10000000;
			}
		}

	  $pArgs .= " --bam $wBAM     --pe ".(isBAMPE($wBAM    ,$peORsr)?"true":"false") if ($wBAM);
	  $pArgs .= " --obj $objRegex --pe ".(isOBJPE($objRegex,$peORsr)?"true":"false") if ($objRegex);
		$pArgs .= " --sra $mySRA    --pe ".(isSRAPE($mySRA   ,$peORsr)?"true":"false") if ($mySRA);
		$pArgs .= " --fq1 $wFQ1 --fq2 $wFQ2 --pe true"                                 if ($wFQ1 && $wFQ2);
		$pArgs .= " --fq1 $wFQ1 --pe false"                                            if ($wFQ1 && !$wFQ2);

		$pArgs .= " --splitSz $fastqSplitSize "                                        if ($type =~ /^(bwa|bwaaln|mm2|bowtie|ssds|ssdsorig)$/);
		$pArgs .= " --aligner $type"                                                   if ($type =~ /^(bwa|bwaaln|mm2|bowtie)$/);
		$pArgs .= ' --type map-ont '                                                   if ($type =~ /^mm2ont$/i);

		$pArgs .= " --r1Len $read1Length --r2Len $read2Length --original true "        if ($type eq 'ssdsorig');

		$pArgs .= " --name $outName --outdir $outDir --genome $genome ";
	}
	$pArgs =~ s/\s+(\S)/ $1/g;

	die("Invalid pipeline (--pipe $type)") unless ($pipeNF);
	die("Invalid args ($pArgs)")           unless ($pArgs);

  ## Should not be done for nf-core pipelines
	$pipeNF = abs_path($pipeNF) unless ($pipeNF =~ /nf-core/);

	return($pipeNF, $pArgs);
}

########################################################################
sub isBAMPE{
	my ($inBAM,$pORs) = @_;

	if ($pORs){
		return 1 if ($pORs eq 'PE');
		return 0;
	}

	## Check if BAM is PE
	my $nPE = `$ENV{DSL2DIR}/samtools view -h $inBAM |head -n 100000 |samtools view -c -f 1 -S /dev/stdin`;
	return ($nPE>0?1:0);
}

########################################################################
sub isOBJPE{
	my ($inObj,$pORs) = @_;

	if ($pORs){
		return 1 if ($pORs eq 'PE');
		return 0;
	}

	## Check if OBJ is PE
	my $nPE = `lsObj --p $inObj |grep -P 'R2.fastq' |wc -l`;
	return ($nPE>0?1:0);
}

########################################################################
sub isSRAPE{
	my ($inSRA,$pORs) = @_;

	if ($pORs){
		return 1 if ($pORs eq 'PE');
		return 0;
	}

	## Check if SRA is PE
	die ("Entrez utils not found: efetch [/usr/local/apps/edirect/10.0/esearch]") unless (-e "/usr/local/apps/edirect/10.0/esearch");
	my $nPE = `/usr/local/apps/edirect/10.0/esearch -db sra -query $inSRA|/usr/local/apps/edirect/10.0/efetch -format runinfo |grep ^$inSRA |grep ,PAIRED, |wc -l`;
	return ($nPE>0?1:0);
}

########################################################################
sub getTime{
	my ($vTime) = @_;

	my ($time4sbatch,$time4nxf);

	if ($vTime =~ /^(\d+)[hH]/){
		my $nHour = $1;
		$time4sbatch = $nHour.":00:03";
		$time4nxf    = $nHour."h";
	}

	if ($vTime =~ /^(\d+):(\d+)/){
		my ($nHour,$nMin) = ($1,$2);
		$time4sbatch = ($nHour+1).":00:03";
		$time4nxf    = ($nHour+1)."h";
	}

	return ($time4sbatch,$time4nxf);
}

########################################################################
sub guessGenome{
	my $file2guess = shift;

	## 08/02/19: KB: allow for genomes with period in name
	$file2guess =~ s/\/$//g;
	$file2guess =~ s/^.+\///g;
	$file2guess =~ s/[\.\/]//g;

	$file2guess = '.'.$file2guess.'.';
	## END 08/02/19: KB

	for my $gName(`ls -d $ENV{NXF_GENOMES}\/\*\/`){
		chomp $gName;
		$gName =~ s/\/$//g;
		$gName =~ s/^.+\///g;

		my $gFolder = $gName;

		$gName =~ s/[\.\/]//g;

		if ($gName =~ /tae/){
			my $xx = 1;
		}
		if ($file2guess =~ /\.($gName)[\.]/i){
			print STDERR "Guessed genome: $gFolder\nFrom file: $file2guess\n";
			return ($gFolder) ;
		}
	}

	return;
}

########################################################################
sub sraOK{
	my $sraList = shift;

#	for my $sraIdentifier(split(/,/,$sraList)){
#		return 0 if ($sraIdentifier !~ /^SRR\d+/)
#	}

	return 0 if ($sraList !~ /^[ESD]RR\d+/);
	return 0 if ($sraList !~ /^[ESDR\d,]+$/); ## Sept 05 19; KB: Only allow commas as delimiters

	return 1;
}

########################################################################
sub getTodayDate {
  my $t = localtime;
  return sprintf( "%04d%02d%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday);
}

########################################################################
sub timestamp {
  my $t = localtime;
  return sprintf( "%04d%02d%02d_%02d%02d%02d_rn%06d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec, rand()*100000);
}

################################################################################################################
sub printARGS{
	my $msg = shift;

	print STDERR "\n\n";
	if ($msg){
		print STDERR "***************************************************************************************************\n";
		print STDERR "INVALID ARGUMENTS: \n";
		print STDERR "---------------------------------------------------------------------------------------------------\n";
		print STDERR $msg."\n";
		print STDERR "***************************************************************************************************\n";
	}

	print STDERR "pipeIt : RDCO pipelining system: Command line arguments: \n";
	print STDERR "-----------------------------------------------------------------------------------------------------------------------\n";
	print STDERR "Input data:\n";
	print STDERR "ONE OF --bam or --obj or --sra or --f1/--f2 are required\n";
	print STDERR "\n";
	print STDERR "--bam         : bam file              SR/PE is autodetcted\n";
	print STDERR "--f1          : read 1 fastq file     SR assumed if only f1 is provided\n";
	print STDERR "--f2          : read 2 fastq file     PE assumed if f1 and f2 are provided\n";
	print STDERR "--sra         : sra ID                retreive a single SRA ID. SR/PE is autodetected\n";
	print STDERR "                                      NOTE: Only for bwa / minimap2 pipelines \n";
	print STDERR "--obj         : sample regex          This will retreive all fastq files matching the regex from \n";
	print STDERR "                                      the RDCO object storage. PE/SR is autodetected. pipeIt will\n";
	print STDERR "                                      die if regex captures > 1 sample. To allow this, use --mOK.\n";
	print STDERR "                                      To test what samples will match your regex, use :\n";
	print STDERR "                                      \$DSL2DIR/lsObj -p \$regex \n";
	print STDERR "                                      NOTE: Only for ssds / ssdsorig / bwa / minimap2 pipelines \n";
	print STDERR "\n";
	print STDERR "-----------------------------------------------------------------------------------------------------------------------\n";
	print STDERR "Pipeline details:\n";
	print STDERR "ONLY --g and --pipe are required\n";
	print STDERR "\n";
	print STDERR "--g           : genome                   from genomes folder (/data/RDCO/genomes) *\n";
	print STDERR "--n           : # of threads             default = 16 \n";
	print STDERR "--m           : memory                   default = 32 \n";
	print STDERR "--pipe        : what pipeline            default = bwa \n";
	print STDERR "                                         others  = bwaaln, minimap2, bowtie, ssds, sdssold, commitfq\n";
  #print STDERR "                                                   ubam, bismark, rtseq\n";
	print STDERR "                                         For multiple pipelines split with commas; commitfq,bwa,ssds\n";
	print STDERR "--o           : output name              default to bam name (-.bam)  \n";
	print STDERR "--od          : output folder            default to output name _ pipeline name \n";
	print STDERR "--s           : RG sample name           default = bam name \n";
	print STDERR "--d           : RG sample date           default = today [".$todayDate."] \n";
	print STDERR "--dep         : dependency job ID        add --dependency=afterany:\$jobid to sbatch \n";
	print STDERR "--ssz         : fastq split size         default = 10000000 \n";
	#print STDERR "--mOK         : allow multi-sample       default = no (see --obj above)\n";
	print STDERR "--gcCorr      : GC-correction for rtSeq  Pre-built options: \n";
	print STDERR "                                         Admera_HiSeqX / NIDDK_HiSeq2500 / Nussenzweig_HiSeq2000 / Yehuda_NextSeq500  \n";
	print STDERR "                                         If not provided, the RT-Seq pipeline will run in LONG-mode. This will generate\n";
	print STDERR "                                         a correction file that can be used for subsequent runs. The best-practice for \n";
	print STDERR "                                         new samples is to generate a correction file from a non-replicating sample \n";
	print STDERR "                                         using LONG-mode. Then use the correction file for subsequent samples.\n";
	print STDERR "\n";
	print STDERR "-----------------------------------------------------------------------------------------------------------------------\n";
	print STDERR "Pipeline logical options (no value required):\n";
	print STDERR "\n";
	print STDERR "--l           : run locally           default = FALSE \n";
	print STDERR "--long        : run for long          sets default runtime to 72 hr \n";
	print STDERR "--inode       : interactive           run on current interactive node \n";
	#print STDERR "--modules     : use modules           use modules instead of singularity \n";
	print STDERR "--hide        : make scripts hidden   default = FALSE; scripts are visible \n";
	print STDERR "--dryrun      : dry run               make scripts but dont start pipe \n";
	print STDERR "--h / --help  : show this HELP        \n";
	print STDERR "--H / --plist : show extra HELP       shows detail of available pipelines \n\n";

}

################################################################################################################
sub printPipeDetails{

	print STDERR "\n";
	print STDERR "***********************************************************************************************************************\n";
	print STDERR "Pipeline details:\n";
	print STDERR "\n";
	print STDERR "--pipe bwa         Align to reference genome using bwa mem 0.7.17. Also generates bigwigs\n";
	print STDERR "                   and sample quality metrics. QC metrics are best viewed in the multiQC\n";
	print STDERR "                   report.\n\n";
	print STDERR "--pipe bwaaln      Align to reference genome using bwa aln 0.7.17. Also generates bigwigs\n";
	print STDERR "                   and sample quality metrics. QC metrics are best viewed in the multiQC\n";
	print STDERR "                   report.\n\n";
	print STDERR "--pipe bowtie2     Align to reference genome using bowtie2. Also generates bigwigs\n";
	print STDERR "                   and sample quality metrics. QC metrics are best viewed in the multiQC\n";
	print STDERR "                   report.\n\n";
	print STDERR "--pipe minimap     Align to reference genome using minimap2/2.13.  Also generates bigwigs\n";
	print STDERR "       minimap2    and sample quality metrics. QC metrics are best viewed in the multiQC\n";
	print STDERR "       mm2         report.\n\n";
	print STDERR "--pipe ssdsorig     ORIGINAL SSDS pipeline using BWA-RA [not deprecated; use ssds pipe instead]\n";
	print STDERR "                   Align SSDS data to reference genome using SSDSpipeline/1.6. This pipeline \n";
	print STDERR "                   generates multiple output bam files, representing different types of DNA that \n";
	print STDERR "                   were identified (ssDNA, dsDNA, unclassified; see Khil et al. 2012). It \n";
	print STDERR "                   also generates BED files that represent the entire sequenced fragments. \n";
	print STDERR "                   These are important for SSDS-based peak calling. The initial bwa-aligned\n";
	print STDERR "                   bam file is also generated. This contains all reads prior to any parsing \n";
	print STDERR "                   for DNA type. The pipeline also generates bigwigs and SSDS sample quality \n";
	print STDERR "                   metrics. QC metrics are best viewed in the multiQC report.\n\n";
	print STDERR "--pipe ssds        Align SSDS data to reference genome using SSDSpipeline_long/1.6\n";
	print STDERR "                   This version of the SSDS pipeline is considerably faster because \n";
	print STDERR "                   instead of the iterative realignment of the second end, it instead\n";
	print STDERR "                   simply allows the fill-in ITR of the second end to be aligned with a\n";
	print STDERR "                   soft-clipped 5 prime end. This method is required for longer reads\n";
	print STDERR "                   like the 150bp PE reads from Admera. Early test show comparable, if\n";
	print STDERR "                   not better alignment than the original SSDS pipeline. For reads >150bp\n";
	print STDERR "                   Minimap 2 is used, otherwise BWA MEM. This was formerly the ssdslong pipe.\n\n";
  print STDERR "--pipe commitfq    This pipeline will save sequencing data to the RDCO object storage. Should ONLY \n";
  print STDERR "                   be used for original sequencing data and NOT for intermediate files. Please \n";
  print STDERR "                   discuss with Kevin if you want to add something to the object store. \n\n";
	# print STDERR "--pipe rtseq       Infer replication timing from WGS data. This pipeline currently requires an \n";
	# print STDERR "                   aligned bam file as input. The replication timing profile will be inferred \n";
	# print STDERR "                   using the method of Koren et al. 2014. Outputs are smoothed bedgraphs with \n";
	# print STDERR "                   either absolute coverage or log2 transformed replication timing.\n";
	# print STDERR "                   ** Use rtseqTest for a test rtSeq run (just 3xCS and pos 1-30000000) **\n\n";
	# print STDERR "--pipe bismark     Align bisulfite sequencing data to a reference genome using bismark. This \n";
	# print STDERR "                   pipeline has been adapted from an Uppsala University nextflow pipe and has \n";
	# print STDERR "                   not been extensively modified. CAUTION should be used and errors reported. \n\n";
	# print STDERR "--pipe ubam        Creates an unaligned BAM and quality metrics. Currently ONLY used for archiving \n";
	# print STDERR "                   Hi-C data. This is a placeholder until a Hi-C pipeline is developed. \n\n";
	print STDERR "***********************************************************************************************************************\n";

	exit;
}
