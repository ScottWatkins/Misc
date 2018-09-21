#!/usr/bin/perl -w
############################################################
# A wrapper script to run vvp/VAAST3, phevor_2, and        #
# graphite.                                                #
#                                                          #
# USAGE: runVVPpipline_v2x.pl -ids <idsfile> --vcf <vcf>   #
#                                                          #
# The ids file MUST be in the following format with the    #
# identifier in column 2. The order does not matter. Input #
# should be a file with a singleton, a trio, or a quad.    #
#                                                          #
# probandID     proband                                    #
# motherID      mother                                     #
# fatherID      father                                     #
# sibID         sibling                                    #
# affectedSibID affected_sibling                           #
#                                                          #
#                                                          #
# VVP, phevor, and bcftools must be in your path           #
# Bcftools (not vcftools) is used for better performance.  #
#                                                          #
# Set the path to the required background and reference    #
# VVP files directly below.  This is the only user         #
# input that is required to set up the script.             #
#                                                          #
# Written: 20160205 WSW  Copyright WSW and the University  #
# of Utah.                                                 #
# Revised: v2_2 for new VVP code: 20170407 WSW             #
# Revised: v2_3 for new VEP, VEP v90                       #
############################################################

use strict;
use Getopt::Long;

##------Set user VVP path and background options here------#

##___Use the full path name for each of these variables____#
#my $bkgvcf = '/archive04/swatkins/genome/vvp_background/chrALL.vep.nocall.phast.sorted.VAT2.vcf.gz';
my $bkgout = '/scratch/ucgd/lustre/u0103307/genome/vvp_v2_background/1KG.050417.vvp.db';
my $bkgnum = '2504';
my $VVP_score_variants = '/uufs/chpc.utah.edu/common/home/u0103307/bin/VVP-transcript_development/VVP';
my $VVP_burden_permutation = '/uufs/chpc.utah.edu/common/home/u0103307/bin/VVP-transcript_development/VAAST3';
my $REFSEQ = '/scratch/ucgd/lustre/u0103307/genome/human_g1k_v37_decoy_phiX.fasta';
##---------------------------------------------------------#

my $VCFTOOLS = (); my $PHEVOR_2 = (); my $bcftools = ();
my $ids = (); my $vcf = (); my $phevor=(); my $hpo = (); my $help = ();
my $genelist = (); my $ext_priors = (); my $label = 'na'; my @name = ();
my $indices = (); my $coding = (); my $snps = (); my $iht = ();
my $pval = 0.999; my $graphite = (); my $cores = 1; my $bcfcores = 1;
my $plot = (); my $pnt = 'i'; my $no_muc = ();
my $cmd = join (' ', $0, @ARGV);
my @var1info = (); my @var2info = (); my @exac = ();


MAIN ();


#----------------------------------------------------------#
#------------------------_MAIN_----------------------------#
#----------------------------------------------------------#

sub MAIN {

  INPUTS ();

  open (LOG, ">>$name[0].log");

  TIME("RUN STARTED.");
  TIME("INFO: Background file: $bkgout");
  TIME("INFO: Background file size  = $bkgnum");
  TIME("INFO: Reference genome: $REFSEQ");
  TIME("COMMAND_LINE: $cmd\n");

  $indices = EXECVVP ();  #return vvp sample indices in hash


  PICK ();

  if (defined $phevor) {

    PHEV ();

    ISECT ();

    GENELIST ();

    if ($graphite) {
      RUN_GRAPHITE ();
      PARSE_GRAPHITE ();
    }

    if ($plot) {
      PLOT ();
    }

  }


#  CLEANDIR ();

  TIME ("RUN FINISHED.");

  close (LOG);

}

#----------------------------------------------------------#
#----------------------end_MAIN----------------------------#
#----------------------------------------------------------#

#----------------------------------------------------------#
#---------------------Subroutines--------------------------#
#----------------------------------------------------------#

#---------------------sub_INPUTS---------------------------#
sub INPUTS {

  GetOptions(

	     "ids=s"          => \$ids,
	     "vcf=s"          => \$vcf,
	     "phevor"         => \$phevor,
	     "hpo=s"          => \$hpo,
	     "genelist=s"     => \$genelist,
	     "ext_priors=s"   => \$ext_priors,
             "coding"         => \$coding,
             "snps"           => \$snps,
	     "iht=s"          => \$iht,
             "label:s"        => \$label,
	     "pval:f"         => \$pval,
	     "graphite=s"     => \$graphite,
	     "cores=i"        => \$cores,
             "bcfcores=i"     => \$bcfcores,
	     "penetrance:s"   => \$pnt,
	     "plot"           => \$plot,
             "no_muc"         => \$no_muc,
	     "help"           => \$help
	    );

  if ((!defined $ids) || (!defined $vcf) || $help  )  {
    USAGE ();
  }

  unless (-e $ids) {
    ERROR ("Input ids file $ids not found.");
  }

  unless (-e $vcf) {
    ERROR ("Input VCF file $vcf not found.");
  }

  if ($phevor) {
    if (!defined $hpo) {
      ERROR ("--hpo <hpo_nodes.txt> required with --phevor option.");
      exit;
    }
    unless (-e $hpo) {
      ERROR ("File $hpo not found. Please check file $hpo.");
    }
  }

  if (defined $coding) {
    $coding = " -c ";
  } else {
    $coding = " ";
  }

  if (!defined $genelist) {
    ERROR ("Current functionality requires --genelist.  Use --genelist best_genes if you don't have a genelist.")
  }

  if ($genelist) {
    my @gl = split /\,/, $genelist;
    foreach my $g (@gl) {
      unless (-e $g) {
	unless ($g eq 'best_genes') {
	  ERROR ("File $genelist not found. Please check file $genelist");
	}
      }
      if ($g eq 'best_genes') {
	unless ($phevor) {
	  ERROR ("best_genes requires --phevor. Please add --phevor and --hpo <hpo_list> to the command line.");
	}
      }
    }
  }

  @name = split /\./, $ids;

  if ($vcf eq "$name[0].recode.vcf") {
    ERROR ("This name is not allowed for the input vcf. Please rename the vcf file.");
  }

  if (!defined $iht) {
    ERROR ("--iht required");
  }

  unless ($iht =~ /[r|d|x]/) {
    ERROR ("--iht mode $iht is not valid. Use d for de novo, r for recessive, x for x-linked recessive.");
  }

  if ($graphite) {

    my @bams = split /\,/, $graphite;

    foreach my $b (@bams) {
      unless ( -e "$b") {
	ERROR ("The BAM file $b was not found. Please check the name and path for $b!");
      }
    }

  }

  my $tbcf = `which bcftools`;
  chomp $tbcf;
  my @tbcf = split /\//,$tbcf;
  unless ($tbcf[-1] eq 'bcftools') {
    ERROR("bcftools was not found. Please put bcftools in your path to use this option.");
  }

  $PHEVOR_2 = `which phevor_2`;
  chomp $PHEVOR_2;
  my @pcheck = split /\//, $PHEVOR_2;
  unless ($pcheck[-1] =~ /phevor/) {
    ERROR("$PHEVOR_2 |phevor_2 is required but not found. Please put phevor_2 in your path.");
  }

  unless (-e "$VVP_score_variants") {
    ERROR("$VVP_score_variants was not found. You must edit this script.\n");
    exit;
  }
  unless (-e "$VVP_burden_permutation") {
    ERROR("$VVP_burden_permutation was not found. You must edit this script.\n");
    exit;
  }

  unless (-e "$bkgout.bit") {
    ERROR("The background file $bkgout  was not found. You must edit this script.\nBe sure to set the background number variable to match the new background file.\n\n");
    exit;
  }

  unless (-e "$REFSEQ") {
    ERROR("The genome reference sequence file $REFSEQ  was not found. You must edit this script.\n\n");
    exit;
  }

}

#----------------------end_sub_INPUTS----------------------#


#----------------------sub_USAGE---------------------------#
sub USAGE {
print"

This script extracts a single sample, trio, or quartet from a master joint-called
VEP-annotated vcf file and runs VVP.  Phevor2 can be run on the output and then
intersected with a list of candidate genes.  The master VCF file must be annotated
with VEP (everything). The user must set the appropriate VVP background in the script.

\n\nUSAGE:  runVVPpipeline.pl --ids <idsfile> --vcf <vep_annotated_mastervcf> --iht [r|d|x]


  --ids       A two column list of the sample, trio, or quad ids and their relations.

              Example:

                  1-02504 proband
                  1-02504-01 mother
                  1-02504-02 father
                  1-02504-03 sibling           <--- if unaffected
                  1-02504-03 affected_sibling  <--- if affected

  --vcf       A vcf file containing the samples listed in the --ids file.
              This file MUST be annotated with the ensemble VEP_90 script,
              (ensembl-tools-release-86/scripts/variant_effect_predictor \
              /variant_effect_predictor.pl), with the \"everything\" switch on.


Options:

  --phevor     Run phevor on the genes with a positive burden score.

  --hpo        <filename> A file of HPO node ids (e.g.HP:012345), one per line,
               to use with phevor.

  --ext_priors <filename> A list of user generated priors to use with phevor

  --genelist   <string> A list of genes to intersect with the burden
               and phevor results.

               Separate multiple lists using commas without spaces
               (e.g. list1,~/others/list2,/home/lists/list3).  The gene lists
               outfiles are named either xxx.intersect or xxx.final.  These
               tab and pipe delimited files can be imported into excel 
               for easy viewing.

               Use \"best_genes\" to annotate the intersect of
               VVP output and top 25 Phevor2 genes.

  --coding     Limit output to coding variants only.

  --snps       Limit output to SNPs only.

  --iht        <char> Set inheritance to recessive, dominant/denovo, or
               x-linked recessive [r|d|x].

  --pnt        Penetrance, default = i. Use c for complete penetrance.
               (not currently implemented in script).

  --label      <string> A user specified label for the output.

  --pval       <float> Maximum allowed p-value for a variant [0.999].

  --graphite   Creates a vcf file of graph adjudicated proband variants
               for the genes in the gene list. Requires indexed bamfiles for
               the proband and both parents with proband listed first
               ( e.g. --graphite I5.bam,I5mom.bam,I5dad.bam ).  Bcftools
               is also required for this function.  Use this feature to
               greatly improve the search for de novo variants.

  --cores      Number of cores to use for VVP, VAAST3, GRAPHITE [1].

  --bcfcores   Number of cores to use for bcftools extraction [1].

  --plot       Create a manhattan-style plot of the phevor output.

  --no_muc     Prevent mucins and other noisy genes from being plotted.

  --help       Print this message.


Examples:

1. Score a trio of exome data to find recessives candidate SNPs:

runVVPpipeline.pl --ids 02505.ids --vcf 02504.vcf --phevor --hpo CTDnodes.list \
 --coding --snps --label CTD_CHD -iht r --genelist best_genes

2. Score a trio of exome data to find de novo candidates snps and indels:

runVVPpipeline_graph.pl --vcf /disk_array/VeryBig.vcf --ids trio01965.ids --phevor --hpo LVOTO_phenotype_list --coding --label LVOTO --bcftools --genelist best_genes --graphite 1-01965.bam,1-01965-01.bam,1-01965-02.bam -iht d
\n\n";

exit;

}
#---------------------end_sub_USAGE------------------------#


#----------------------sub_EXECVVP-------------------------#
sub EXECVVP {

  open (IN, "<$ids") || die "$ids not found\n";
  my  %ids = (); my $cids = 0;
  while (<IN>) {
    chomp;
    my @i = split /\t|\s|,/, $_;
    $ids{$i[0]}=$i[1];
    $cids++;
  }
  close (IN);

##--generate new vcf with samples to be run.

  TIME ("Extracting $ids from $vcf ...");

    my $bcfsamples = join (",", keys(%ids));
    my $bcfcmd = ();
    $bcfcmd = "bcftools view --threads $bcfcores -Oz -s $bcfsamples -e \'ALT=\"*\"\' $vcf >$name[0].recode.vcf.gz";
    TIME ("Using bcftools to extract samples from $vcf ...");
    system ("$bcfcmd");
    system ("bcftools index $name[0].recode.vcf.gz");

#    TIME("Removing all * alleles from the extracted VCF file $name[0].recode.vcf.gz");
#    system("gzip -dc $name[0].recode.vcf.gz | awk  \'{if(\$5 \!\= \"\*\") print }\' | bgzip --threads $bcfcores >$name[0].filtered.vcf.gz");
#    sleep(3);
#    system("mv $name[0].filtered.vcf.gz $name[0].recode.vcf.gz");
#    system("bcftools index $name[0].recode.vcf.gz");


  my @p = ();
  my $sn = `bcftools query -l $name[0].recode.vcf.gz`;
  @p = split /\n/, $sn;

  ##--get index numbers for vvp config

  TIME ("Getting indices for $ids and generating the vvp command ...");
    my %idx = ();
  foreach my $k (keys(%ids)) {
    my $v = $ids{$k};
    my $c = 0;
    foreach my $x ( @p[0..$#p]) {
      if ($x eq $k) {
	$idx{$v}=$c;
      }
      $c++;
    }
  }

  my $pb=(); my $mo=(); my $fa=(); my $sib=(); my $asib=(); my $sibstat=();

  while ((my $k, my $v) = (each(%idx))) {
    if ($k =~ /pro/i) {
      $pb = $v;
    }
    if ($k =~ /mo/i) {
      $mo = $v;
    }
    if ($k =~ /fa/i) {
      $fa = $v;
    }
    if ($k eq "sibling") {
      $sib = $v;
      $sibstat = 0;
   }
    if ($k eq "affected_sibling") {
      $sib = $v;
      $sibstat = 1;
    }
  }

  my $filters = ();
  if ($cids eq 3) {
    $filters = 't';
  }  elsif ($cids eq 4) {
    $filters = 'q';
  }  elsif ($cids eq 1) {
    $filters = 'n';
  } else {
    ERROR("Invalid number of samples listed in file.  Check the ids file!");
    exit;
  }


  open (OUT, ">vvp.com");
  print OUT "$VVP_score_variants -d $bkgout -i $name[0].recode.vcf.gz -v CSQ,3,6,0,15 $coding -n $cores -o $name[0].scored_variants.out >$name[0]_parse_VVP.out\n";

  close (OUT);

  TIME ("Scoring variants with VVP ...");
  system ("chmod 700 vvp.com;  ./vvp.com");


  TIME ("Sorting the VVP output ...");

  #*****Sorting_Workaround_per_Javier|Steve_*********
  system ("sort -k2,2 $name[0].scored_variants.out >$name[0].scored_variants.tmp");
  sleep (2);
  system ("mv $name[0].scored_variants.tmp $name[0].scored_variants.out");


  open (VAAST, ">vaast.com");

  if ($cids == 1) {
    print VAAST "$VVP_burden_permutation -i $name[0].scored_variants.out -d $bkgout -t $cids -b $bkgnum -n $cores -e $iht -f $filters -r $pb  >$name[0].burden.out 2>$name[0].burden.error\n";
  } elsif ($cids == 3) {
    print VAAST "$VVP_burden_permutation -i $name[0].scored_variants.out -d $bkgout -t $cids -b $bkgnum -n $cores -e $iht -f $filters -r $pb -m $mo -w $fa >$name[0].burden.out 2>$name[0].burden.error\n";
  } elsif ($cids == 4) {
    print VAAST "$VVP_burden_permutation -i $name[0].scored_variants.out -d $bkgout -t $cids -b $bkgnum -n $cores -e $iht -f $filters -r $pb -s $sib -u $sibstat -m $mo -w $fa >$name[0].burden.out 2>$name[0].burden.error\n";
  }


  close (VAAST);


  TIME ("Running burden permutations ...");

  system ("chmod 700 vaast.com; ./vaast.com");



  my %T = ();  #Lookup hash for ensembl transcript to gene name

  open (TRLU, "<$name[0]_parse_VVP.out");
  while (<TRLU>) {
    my @l = split /\t/, $_;
    $T{$l[5]} = $l[4];
  }
  close (TRLU);

  open (BLU, "<$name[0].burden.out");
  open (BLUOUT, ">$name[0].burden.out2");
  while (<BLU>) {
    chomp;
    my @l = split /\t/, $_;
    if (exists ($T{$l[0]})) {
      my $gname = $T{$l[0]};
      my $trid = $l[0];
      $l[0] = $gname;    #rename to gene but retain transcript id.
      my $g = join ("\t", @l, $trid, "\n");
      print BLUOUT "$g";
    }
  }
  close BLUOUT;
  close BLU;


  TIME ("Cleaning burden results ...", "Reordering scores and p-values for phevor input ...");
  open (BURDIN, "<$name[0].burden.out2");
  open (BURDOUT, ">$name[0].burden.tmp");
  while (<BURDIN>) {
    chomp;
    my @b = split /\t/, $_;
    if (($b[1] > 0) && (($b[2] > -0.0000001) && ($b[2] < $pval)) ) {
      print BURDOUT "$b[0]\t$b[2]\t$b[1]\t$b[3]\t$b[4]\t$b[5]\n";
    }
  }
  close (BURDIN);
  close (BURDOUT);

  TIME ("Sorting burden results ...");
  system ("sort -k2,2n -k3,3nr $name[0].burden.tmp >$name[0].burden.sorted.tmp");
  open (FIN, "<$name[0].burden.sorted.tmp");
  open (FOUT, ">$name[0].burden.sorted");
  while (<FIN>) {
    print FOUT "$.\t$_";
  }
  close (FIN);
  close (FOUT);


  TIME("Finished VVP and VAAST runs.");

  return (\%idx);

}
#----------------------End_sub_EXEVVP----------------------#


#------------------------sub_PICK--------------------------#
# This sub takes the sorted burden file and picks the      #
# lowest p-val and highest clrt score for multi-transcript #
# burden outfile and outputs the xxx.burden.picked file.

sub PICK {

  TIME ("Selecting the best single transcripts from $name[0].burden.sorted.");

  my %P = ();

  open (IN, "<$name[0].burden.sorted");
  open (OUT, ">$name[0].burden.picked.tmp");
  while (<IN>) {
    chomp;
    my @r = split /\t/, $_;
    if (exists ($P{$r[1]})) {
      my $v = $P{$r[1]};
      my $pval = $v->[2];
      my $clrt = $v->[3];
      if ($pval >  $r[2]) {    #take lowest p-val first, the lowest clrt
	$P{$r[1]} = [ @r ];
      }
      if ( ($pval > $r[2]) && ($clrt < $r[3]) ) {
	$P{$r[1]} = [ @r ];
      }
    } else {
     $P{$r[1]} = [ @r ];
   }
  }


  while ( (my $k, my $v) = each (%P) ) {
    my $picked = join("\t", @$v);
    print OUT "$picked\n";
  }

  system ("sort -k1,1n $name[0].burden.picked.tmp >$name[0].burden.picked.tmp2");

  open (IN2, "<$name[0].burden.picked.tmp2");
  open (OUT2,">$name[0].burden.picked");
  my $c = 0;
  while (<IN2>) {
    chomp;
    $c++;
    my @r = split /\t/, $_;
    $r[0] = $c;
    my $ordered = join("\t", @r);
    print OUT2 "$ordered\n";
  }

  close (IN);
  close (OUT);
  close (IN2);
  close (OUT2);
}

#------------------------End_PICK--------------------------#


#------------------------sub_PHEV--------------------------#
# This sub runs phevor on the burden.picked output

sub PHEV {
  if ($phevor) {
    TIME ("Running phevor with terms from $hpo ...");
    system ("cut -f2,3,4 $name[0].burden.picked >$name[0].phevor.in");
    open (PIN, "<$hpo");
    my @hpo = ();
    while (<PIN>) {
      chomp;
      push (@hpo, $_);
    }
    my $nodes = join(",", @hpo);
    if ($ext_priors) {
      system("$PHEVOR_2 -c /uufs/chpc.utah.edu/common/home/u0103307/bin/Phevor2/phevor.omicia.conf -i $nodes -b GO -o $name[0].phevor.in -e $ext_priors -d $name[0].phevor2 2>phevor.error" );
    } else {
      system("$PHEVOR_2 -c /uufs/chpc.utah.edu/common/home/u0103307/bin/Phevor2/phevor.omicia.conf -i $nodes -b GO -o $name[0].phevor.in  -d $name[0].phevor2 2>phevor.error" );
      TIME("$PHEVOR_2 -c /uufs/chpc.utah.edu/common/home/u0103307/bin/Phevor2/phevor.omicia.conf -i $nodes -b GO -o $name[0].phevor.in  -d $name[0].phevor2 2>phevor.error");
    }
  }

}
#---------------------end_sub_PHEV-------------------------#

#---------------------sub_ISECT----------------------------#
# This sub intersects all genes from the VVP burden picked 
# output
# with the best phevor2 genes and writes them to best_genes
# for later annotation if "best_genes" is given as $genelist.
# Currently, only the top 100 PHEVOR gene are considered.

sub ISECT {

  my $best = 100;

  open (VIN, "<$name[0].burden.picked");
  open (PHIN, "<$name[0].phevor2");

  my %vin = (); my %phin = (); my $vinvars = 0;
  while (<VIN>) {        # grab all vaast genes in vin hash
    if ((/^[0-9]/)) {
      $vinvars++;
      my @l = split /\t|s+/, $_;
      if (defined $l[1]) {
	$vin{ $l[1]} = '1';
      }
    }
  }


  my @pscores = (); my @phin = ();
  while (<PHIN>) {       #grab top phevor genes in phin hash
    if ((/^[0-9]/) && ($. <= $best) ) {
      my @l = split /\t|\s+/, $_;
      if (defined $l[1]) {
	$phin{$l[1]} = '1';
	push (@phin, $l[1]);
	push (@pscores, $l[2]);
      }
    }
  }
  close (VIN);
  close (PHIN);


  my $pdiff = $pscores[0] - $pscores[1];
  my $ppdiff = sprintf("%4.2f", $pdiff);
  if ($pdiff > 1.5) {
    open (NOUT, ">Notice.txt");
    print NOUT "NOTICE::Phevor differential for $phin[0] is $ppdiff.  $phin[0] is likely to be the causitive gene!!\n";
    close (NOUT);
  }

my %isect = ();
  foreach (keys %vin){
   $isect{$_} = $vin{$_} if exists $phin{$_};
 }

  open (TOP, ">best_genes");

  foreach my $k (keys %isect) {
      print TOP "$k\n";
  }
  close (TOP);

}

#---------------------end_sub_ISECT------------------------#


#---------------------sub_GENELIST-------------------------#
sub GENELIST {

  my @genelist = split /\,/, $genelist;  #pull in file $genelist

  foreach my $glist (@genelist) {

    TIME ("Intersecting outputs with genes from $glist ...");

    my @listname = split /\//, $glist;
    my $listname = $listname[-1];

    my %list=();
    open (IN, "<$glist");
    while (<IN>) {
      chomp;
      $list{$_} = 1;
    }
    close (IN);

    my @vvplist = (); my @phevlist = ();

    open(IN2, "<$name[0].burden.picked");

    while (<IN2>) {
      my @i = split /\t/, $_;
      if (exists $list{$i[1]}) {
	push (@vvplist, $i[1]);
      }
    }
    close (IN2);

    if (-e "$name[0].phevor2") {
      open(IN3, "<$name[0].phevor2");
      while (<IN3>) {
	my @j = split /\s+\t/, $_;
	if (exists $list{$j[1]}) {
	  push (@phevlist, $j[1]);
	}
      }
      close (IN3);
    }

    my @vvponly = ();

    my %lup = ();
    @lup{@phevlist} = ();
    foreach my $x (@vvplist) {
      push(@vvponly, $x) unless exists $lup{$x};
    }

    open (OUT, ">$name[0].$listname.tmp");
    foreach my $i (@phevlist) {
      print OUT  "$i\n";
    }
    foreach my $i (@vvponly) {
      print OUT "$i\n";
    }
    close (OUT);

    open (OUT2, ">$name[0].$listname.intersect"); #reprint indices hash in output
    print OUT2 "#Sample\tUserLabel\tGene\tVVP_rank\tPhevor_rank\tVaast_pval\tCLRT\tPhevor_score\tInheritance\tAllele1\tMaxAF1\tAllele2\tMaxAF2\tExpGTfreq\ttype1|impact1|HGVSc1|HGVSp1|rsid1|sift1|polyphen1|ClinSig1|GenePheno\ttype2|impact2|HGVSc2|HGVSp2|rsid2|sift2|polyphen2|ClinSig2|GenePheno\tgt1|gt2\n";
    print LOG  "INFO::";
    while ((my $k, my $v) = (each(%$indices))) {
      print LOG "INDEX:$k=$v ";
    }
    print LOG "\n";
    print "Processing $name[0].$listname.tmp...\n";

    open (IN4, "<$name[0].$listname.tmp");
    open (GRAPHVARS, ">>graphite_vars.in");

    while (<IN4>) {
      chomp;

      my $bgenes = `grep -P '$_\t' $name[0].burden.picked`; #trailing tab needed
      chomp ($bgenes);

      my @bgenes = split /\t/, $bgenes;

      my $transcript = $bgenes[6];

      my @pinfo = ();
      my $phvrank = 'na';

      if (-e "$name[0].phevor2") {
	my $pinfo =  `grep -P '$_' $name[0].phevor2`;
	@pinfo = split /\s+\t/, $pinfo;
	$phvrank = $pinfo[0] + 1;   # convert to 1-based rank
      }

      #Processing of the variants

      # Note: for the denovo case, sometimes the single variant
      # is in the VVP 2nd position. Since there is only one variant,
      # the 2nd is swapped into the first position for processing.

      if ($iht eq 'd') {
	if ($bgenes[4] !~ /[0-9]/) {
	  my $swap = $bgenes[4];
	  $bgenes[4] = $bgenes[5];
	  $bgenes[5] = $swap;
	}
      }

      my @vars1 = split /\,/, $bgenes[4];
      my @vars2 = split /\,/, $bgenes[5];

      my $tvar1 = (); my $tvar2 = (); my $sv1 = (); my $sv2 = ();
      my @sv1 = (); my @sv2 = (); my $allele1 = (); my $allele2 = ();
      my $csq1 = (); my $csq2 = ();



      foreach my $i (@vars1) {
	if ($i =~ /\d/) {
	  $tvar1 = $i;

	  # decomposing the vcf can yield multiple variants with
	  # the same position. Split on \n and the take only
	  # the one matching the burden.picked. Filter the * alt alleles,
	  # and missing data, then make sure there's just one allele by
          # grepping the actual site from the burden.picked file.

          # Now add transcript level parsing of the VEP annotation
          # by selecting the correct transcript with the damaged
          # variant(s).

	  $sv1 = `bcftools view -H  -r $tvar1 $name[0].recode.vcf.gz | grep $tvar1`;


	  my @sv1 = split /\n/, $sv1;

	  if (scalar(@sv1) > 1 ) {
	    @sv1 = ();
	    $allele1 = 'multiallelic';;
	    last;
	  } elsif (scalar(@sv1) == 1) {
	    my @r = split /\t/, $sv1[0];


	    my @tranno1 = split /,/, $sv1;
	    foreach (@tranno1) {
	      if (/$transcript/) {
		$csq1 = $_;        #transcript-based csq1
	      }
	    }

	    @sv1 = ();
	    @sv1 = @r;
	    $allele1 = $sv1[4];
	    last;
	  } else {
	  $allele1 = 'missing_data'
	  }
	}
      }

      foreach my $i (@vars2) {   # eg. .,12:123456,.

	if ($i =~ /\d/) {
	  $tvar2 = $i;

	  $sv2 = `bcftools view -H  -r $tvar2 $name[0].recode.vcf.gz | grep $tvar2`;

	  my @sv2 = split /\n/, $sv2;

	  if (scalar(@sv2) > 1 ) {
	    @sv2 = ();
	    $allele2 = 'multiallelic';;
	    last;
	  } elsif (scalar(@sv2) == 1 ) {
	    my @r = split /\t/, $sv2[0];


	    my @tranno2 = split /,/, $sv2;
	    foreach (@tranno2) {
	      if (/$transcript/) {   #one transcript for both alleles
		$csq2 = $_;          #transcript-based csq2
	      }
	    }

	    @sv2 = ();
	    @sv2 = @r;
	    $allele2 = $sv2[4];
	    last;
	  } else {
	    $allele2 = 'missing_data';
	  }
	} else {
	  $tvar2 = 'na';
	  $allele2 = 'na';
	  $csq2 = 'na';
	}
      }

#     print "\n\n****$tvar1,$tvar2,$allele1,$allele2\n";

#     feed positions and alleles to PARSE_VARS and return data

      (my $exac_ref, my $var1info_ref, my $var2info_ref, my $gts)= PARSE_VARS ($tvar1,$tvar2,$allele1,$allele2,$csq1,$csq2);


      my $expgeno = ();
      my $type = ();

      if ($iht eq 'd') {
	$type = 'de_novo';
	$expgeno = @$exac_ref[0];
      } elsif ($iht eq 'r') {
	if ($tvar2 eq 'na') {
	  $type = 'recessive';
	} else {
	  $type = 'compound_het';
	}
      }

      if ($type eq 'recessive') {
	if (@$exac_ref[0] ne 'na') {

	  $expgeno = @$exac_ref[0]**2;
	} else {
	  $expgeno = 0.000001**2;
	}
      }

      if ($type eq 'compound_het') {
	if ((@$exac_ref[0] ne 'na') && (@$exac_ref[1] ne 'na')) {
	  $expgeno = @$exac_ref[0] * @$exac_ref[1];
	} elsif (@$exac_ref[0] ne 'na') {
	  $expgeno = @$exac_ref[0] * 0.000001;
	} elsif (@$exac_ref[1] ne 'na') {
	  $expgeno = 0.000001 * @$exac_ref[1];
	} else {
	  $expgeno = 0.000001**2;
	}
      }

      my $phevor2 = `grep $bgenes[1] $name[0].phevor2`;
      my @phevor2 = split /\s+|\t/, $phevor2;

      my $var1info = join("|", @$var1info_ref);
      my $var2info = join("|", @$var2info_ref);

      print OUT2 "$name[0]\t$label\t$bgenes[1]\t$bgenes[0]\t$phvrank\t$bgenes[2]\t$bgenes[3]\t$phevor2[2]\t$type\t$tvar1\t@$exac_ref[0]\t$tvar2\t@$exac_ref[1]\t$expgeno\t$var1info\t$var2info\t$gts\n";


      (my $v1 = $tvar1) =~ s/\:/\t/; # create a list of vars to extract for graphite
      print GRAPHVARS "$v1\n";        # list will include all vars from all gene lists
      if ($tvar2 ne 'na') {
	(my $v2 = $tvar2) =~ s/\:/\t/;
	print GRAPHVARS "$v2\n";
      }


    }

    close (IN4);
    close (OUT2);
    close (GRAPHVARS);
  }
}

#-----------------------end_sub_GENELIST-------------------#

#-----------------------sub_PARSE_VARS---------------------#

# This subroutine scans the recode.vcf and parses the CSQ
# annotations for 2 variant positions and their alleles.
# The results are returned as references
# in the 3 global variables @exac, @var1info and @var2info
# and printed in subroutine GENELIST.
# The returned info may be changed to accomondate the desired
# CSQ fields depending on what infomation you want. Alter
# the @csq array below. See vep90 documentation for
# field explainations.

sub PARSE_VARS {


  my $a1 = (); my $a2 = ();
  my $csq1 = $_[4]; my $csq2 = $_[5];
  my @gts1 = (); my @gts2 = ();

  $a1 = $_[0];

  unless ($_[1] =~ /[n|m]/) {
    $a2 = $_[1];
  }

  @exac = ();  @var1info = ();  @var2info = (); #reinitialize globals


  # CSQ fields used from vep90 --everything 1=Consequence, 2=IMPACT,
  # 10=HGVSc, 11=HGVSp, 17=Existing_variation, 32=GENE_PHENO, 33=SIFT, 
  # 34=PolyPhen
  # 56=CLIN_SIG, 58=PHENO, 59=PUBMED; USE 54=MAX_AF gnomAD below for screening
  my @csq = qw(1 2 10 11 17 33 34 56 58);

  my @info1 = (); my @info2 = (); my $exac1 = (); my $exac2 = ();

  my $mvar1 = `bcftools view -H -r $a1 $name[0].recode.vcf.gz | grep $a1`;

  my @VARfields1 = split /\t|;/, $mvar1;

  foreach my $a (@VARfields1) {   #vt annotation adds new fields to vcf so
    if ($a =~ /^CSQ/) {           #get the CSQ fields only
      $mvar1 = $a;
    }
    if ( $a =~ /^([0|1]\/[0|1])\:/ )  { #get genotypes in order
      push (@gts1, $1);
    }
  }

  my @p1  = split /CSQ\=/, $mvar1;

  if (!defined $p1[1]) {
    @info1 = ('No data', 'Variant_fails_filters');
    $exac1 = 'na';
    @gts1 = ('na');
  } else {
    if (defined $csq1) {
      @info1 = split /\|/, $csq1;    #input from transcript parsed $csq1 input
    }
  }

  my $gts1 = join (",", @gts1);

  my $gts2 = 'na';                #default recessive

  if (!defined $info1[54]) {      #exac total freq
    $exac1 = 'na';
  } elsif (length($info1[54]) == 0) {
    $exac1 = 'na';
  } elsif (length($info1[54]) >= 1) {

    $exac1 = $info1[54];

  } else {
    $exac1 = 'na';
  }

  foreach my $i (@csq) {
    if (!defined $info1[$i]) {
      push (@var1info, 'na');
    } elsif (length($info1[$i]) >= 1) {   #!!! some defined blank values
      push (@var1info, $info1[$i]);
    } else {
      push (@var1info, 'na');
    }
  }

  my @p2 = ();
  my $mvar2 = ();

  if ($_[1] eq 'na') {
    @var2info = ('na|na|na|na|na|na|na|na|na');
  }

  if ($_[1] ne 'na') {

    $mvar2 = `bcftools view -H -r $a2 $name[0].recode.vcf.gz | grep $a2`;

    my @VARfields2 = split /\t|;/, $mvar2;

    foreach my $a (@VARfields2) {   #vt annotation adds new fields to vcf,
      if ($a =~ /^CSQ/) {           #get the CSQ fields only
	$mvar2 = $a;
      }

      if ($a =~ /^([0|1]\/[0|1])\:/)  { #get genotypes in order
	push (@gts2, $1);
      }
    }

    @p2  = split /CSQ\=/, $mvar2;

    if (!defined $p2[1]) {
      @info2 = ('No data', 'Variant_fails_filters');
      $exac2 = 'na';
      @gts2 = ('na');
    } else {
      if (defined $csq2) {
	@info2 = split /\|/, $csq2;    #transcript $csq2 from above
      }
    }


    $gts2 = join (",", @gts2);

    if (!defined $info2[54]) {
      $exac2 = 'na';
    } elsif (length($info2[54]) == 0) {
      $exac2 = 'na';
    } elsif (length($info2[54]) >= 1 ) {

      $exac2=$info2[54];

    } else {
      $exac2 = 'na';
    }

    foreach my $i (@csq) {
      if (!defined $info2[$i]) {
	push (@var2info, 'na');
      } elsif (length($info2[$i]) == 0) {
	push (@var2info, 'na');	
      } elsif (length($info2[$i]) >= 1) {
	push (@var2info, $info2[$i]);
      } else {
	push (@var2info, 'na');
      }
    }
  } else {
    $exac2 = 'na';
  }

  my $gts = join ("|", $gts1, $gts2);

  @exac = ($exac1, $exac2);



  return (\@exac, \@var1info, \@var2info), $gts;

}
#---------------------end_PARSE_VARS-----------------------#


#-------------------------sub_RUN_GRAPHITE---------------------#
# This subroutine examines the variants in the xxx.intersect
# output (makes a vcf file) and runs graphite on the vcf file.


sub RUN_GRAPHITE {

  TIME("Running graphite on candidate variants ...");

  my $BAMS = ();


  ($BAMS = $graphite) =~ s/\,/ /g;
  my @bams = split /\,/, $graphite;


  my @check = split /\//, $bams[0];

  open (CHECK, "<$ids");   #make sure that the first bam is the proband
  
  while (<CHECK>) {
    chomp;
    my @c = split /\t|\s+/, $_;
    if ($c[1] eq 'proband') {
      my $pbbam = join (".", $c[0], 'bam');

#      print "--->$pbbam and $check[-1]\n";

      unless ($pbbam eq $check[-1]) {
	print "The sample name of the first listed command line bam \"$bams[0]\" is not the same as the name expected for the proband bam  \"$pbbam\".\n  The first bam file must be the proband, and it should have the same name as the bam (e.g. patient1.bam and patient1). Please correct the names.\n";
	exit;
      }
    }
  }

# Make the list of valid variants to test for either de novo or recessive
# Exclude the * alleles, non-polymorphic sites, non-genotyped sites.
# Remake the graphite_vars.in.  Allow missing geno ./. in extraction
# because transhets don't require parent to valid. Clean the file
# to get rid of extraneous lines with invalid formats.

  system ("sort -k1,1n -k2,2n  graphite_vars.in |uniq > sorted_vars.tmp");  #remove duplicates for graphite runs

  open (CLEAN, "<sorted_vars.tmp");
  open (CLEANED, ">graphite_vars.in");

  while (<CLEAN>) {
    if (/^[0-9|X|Y][0-9]?\t[0-9]+$/) {
      print CLEANED "$_";
    }
  }

  close (CLEAN);
  close (CLEANED);

  TIME("Preparing all candidate sites ...");
  system ("bcftools view  -R graphite_vars.in $name[0].recode.vcf.gz  >tmp1 ");

  open (GRAF, "<tmp1");
  open (GRAFOUT, ">graphite_in.vcf");

  my %seen = ();
  while (<GRAF>) {
    my @l = split /\t/, $_;

    if (/^#/) {
      print GRAFOUT "$_";
    } else {
      my $site = join (":", $l[0], $l[1]);
      $seen{$site}++;

      if ($seen{$site} > 1) {   #omit additional  decomposed sites that 
                                #have same position - take first shortcut!?!
	next;
      }

      unless (/\t\*\t/) {
	print GRAFOUT "$_";
      }
    }
  }
  close (GRAF);
  close (GRAFOUT);

  unshift (@bams, " ");
  my $GBAMS = join (" -b ", @bams);

  system ("mkdir graphout");


  TIME("Running graphite:  COMMAND_LINE: graphite -t $cores -d -f $REFSEQ  $GBAMS -v graphite_in.vcf -o graphout");
  system ("graphite -t $cores -d -f $REFSEQ $GBAMS -v graphite_in.vcf -o graphout");

  system ("mv graphout/graphite_in.vcf graphite_out.vcf");
  system ("rmdir graphout");
  TIME("Graphite run completed.")

}

#-------------------------end_RUN_GRAPHITE---------------------#

#-------------------------sub_PARSE_GRAPHITE-------------------#
# This sub take the graphite vcf and an applies depth and strand
# coverage thresholds to determine positive adjudication status.
# These positive sites are then pulled from the list of potential
# variants, and a final tables is created.

sub PARSE_GRAPHITE {

  TIME("Determining the adjudication status of each variant.");

  my %ids = (); my $pbidx = (); my $moidx = (); my $faidx =();
  my %ajudvars = (); my $pb_dpnfp = (); my $mo_dpnfp = ();
  my $fa_dpnfp =(); my @pb_dp4nfp = (); my @mo_dp4nfp = ();
  my @fa_dp4nfp = (); my $sibidx = (); my @sib_db4nfp = ();

  open (ADJUD, ">Adjudication.out");

  open (IDS, "<$name[0].ids");  #read in id file
  while (<IDS>) {
    chomp;
    my @i = split /\t|\s+/, $_;
    $ids{$i[0]} = $i[1];
  }
  close (IDS);

  open (GRIN, "<graphite_out.vcf");  #get indices of proband,parents
  while (<GRIN>) {
    chomp;
    if (/^#C/) {
      my @r = split /\t/, $_;
      my @g = @r[9 .. $#r];
      for (my $n = 0; $n <= $#g; $n++) {
	my $v = $ids{$g[$n]};

	if ($v eq 'proband') {
	  $pbidx = $n;
	}
	if ($v eq 'mother') {
	  $moidx = $n;
	}
	if ($v eq 'father') {
	  $faidx = $n;
	}
	if ($v eq 'sibling') {
	  $sibidx = $n;
	}
	if ($v eq 'affected_sibling') {
	  $sibidx = $n;
	}
      }
    }

    unless (/#/) {
      my @r = split /\t/, $_;
      my @g = @r[9 .. $#r];

      my @id = split /:/, $r[8];  #get the indices of the fields of interest
      my $dpnfpidx = (); my $dp4nfpidx = ();
      for (my $c = -1; $c <= $#id; $c++) {
	if ($id[$c] eq 'DP_NFP') {
	  $dpnfpidx = $c;
	}
	if ($id[$c] eq 'DP4_NFP') {
	  $dp4nfpidx = $c;
	}
      }

      #----------------Trios_section---------------#

      if (scalar(@g) == 3) {     #trios
	my @pb_rdat = split /:/, $g[$pbidx];
	$pb_dpnfp = $pb_rdat[$dpnfpidx];
	@pb_dp4nfp = split /,/, $pb_rdat[$dp4nfpidx];

	my @mo_rdat = split /:/, $g[$moidx];
	$mo_dpnfp = $mo_rdat[$dpnfpidx];
	@mo_dp4nfp = split /,/, $mo_rdat[$dp4nfpidx];

	my @fa_rdat = split /:/, $g[$faidx];
	$fa_dpnfp = $fa_rdat[$dpnfpidx];
	@fa_dp4nfp = split /,/, $fa_rdat[$dp4nfpidx];

	if ($iht eq 'r') {

	  if ( ($pb_dpnfp >= 6) && (($mo_dpnfp >= 6) || ($fa_dpnfp >= 6)) ) {
	    if (  (($pb_dp4nfp[2] >= 1) || ($pb_dp4nfp[3] >= 1 )) && (   (($mo_dp4nfp[2] >= 2) && ($mo_dp4nfp[3] >= 2 )) ||  (($fa_dp4nfp[2] >= 2) && ($fa_dp4nfp[3] >= 2 ))    )   ) {
	      my $pass = join (":", $r[0],$r[1]);

	      print ADJUD "Positive graphite ajudication (95%) for $r[0]:$r[1] => Sufficient depth    => proband:$pb_dpnfp:@pb_dp4nfp,mother:$mo_dpnfp:@mo_dp4nfp,father:$fa_dpnfp:@fa_dp4nfp\n";


	      if (  (($pb_dp4nfp[0] <= 2) || ($pb_dp4nfp[1] <= 2 )) && (   (($mo_dp4nfp[0] <= 2) || ($mo_dp4nfp[1] <= 2 )) ||  (($fa_dp4nfp[0] <= 2) || ($fa_dp4nfp[1] <= 2 ))    )   ) {
		print ADJUD "Negative graphite ajudication (95%) for $r[0]:$r[1] => Bad Strandedness  => proband:$pb_dpnfp:@pb_dp4nfp,mother:$mo_dpnfp:@mo_dp4nfp,father:$fa_dpnfp:@fa_dp4nfp\n";
		next;
	      }

	      $ajudvars{$pass} = '1';  # add pos vars to ajud hash for later


	    } else {
	      print ADJUD "Negative ajudication (95%) for $r[0]:$r[1] => Sufficient depth   => proband:$pb_dpnfp:@pb_dp4nfp,mother:$mo_dpnfp:@mo_dp4nfp,father:$fa_dpnfp:@fa_dp4nfp\n"
	    }

	  } else {
	      print ADJUD "Negative ajudication (95%) for $r[0]:$r[1] => Insufficient depth => proband:$pb_dpnfp:@pb_dp4nfp,mother:$mo_dpnfp:@mo_dp4nfp,father:$fa_dpnfp:@fa_dp4nfp\n"
	  }

	}

	if ($iht eq 'd') {
	  #print "pb = $pb_dpnfp mo = $mo_dpnfp fa = $fa_dpnfp\n";
	  my $falim = 0;
	  my $molim = 0;
	  if ($fa_dpnfp > 50) {
	    $falim = $fa_dpnfp * 0.03;
	  }
	  if ($mo_dpnfp > 50) {
	    $molim = $mo_dpnfp * 0.03;
	  }

	  my $pblim = 1;
	  if ($pb_dpnfp >= 33) {
	    $pblim = $pb_dpnfp * 0.075;
	  }

	  if ( ($pb_dpnfp >= 8) && ($mo_dpnfp >= 8) && ($fa_dpnfp >= 8) ) {
	    if (  (($pb_dp4nfp[2] >= $pblim) && ($pb_dp4nfp[3] >= $pblim )) && (   (($mo_dp4nfp[2] <= $molim) && ($mo_dp4nfp[3] <= $molim )) &&  (($fa_dp4nfp[2] <= $falim) && ($fa_dp4nfp[3] <= $falim ))    )   ) {

	      my $pass = join (":", $r[0],$r[1]);
	      $ajudvars{$pass} = '1';

	      print ADJUD "Positive ajudication (95%) for $r[0]:$r[1] => Sufficient depth   => proband:$pb_dpnfp:@pb_dp4nfp,mother:$mo_dpnfp:@mo_dp4nfp,father:$fa_dpnfp:@fa_dp4nfp\n"

	    } else {
	      print ADJUD "Negative ajudication (95%) for $r[0]:$r[1] => Sufficient depth   => proband:$pb_dpnfp:@pb_dp4nfp,mother:$mo_dpnfp:@mo_dp4nfp,father:$fa_dpnfp:@fa_dp4nfp\n"
	    }

	  } else {
	      print ADJUD "Negative ajudication (95%) for $r[0]:$r[1] => Insufficient depth => proband:$pb_dpnfp:@pb_dp4nfp,mother:$mo_dpnfp:@mo_dp4nfp,father:$fa_dpnfp:@fa_dp4nfp\n"
	  }

	}
      }

      #----------------End_Trio_section----------------#


      #----------------Singleton_section---------------#

      if (scalar (@g) == 1) {
#	TIME("Adjudicating the single proband's variants.");

        my @pb_rdat = split /:/, $g[$pbidx];
        $pb_dpnfp = $pb_rdat[$dpnfpidx];
        my @pb_dp4nfp = split /,/, $pb_rdat[$dp4nfpidx];

	print "@pb_dp4nfp\n";


	if ($iht eq 'r') {

	  if ( $pb_dpnfp >= 10 ) {

	    if ( ($pb_dp4nfp[2] >= 3) || ($pb_dp4nfp[3] >= 3) )  {    #alt allele has 3 reads on either strand
	      my $pass = join (":", $r[0],$r[1]);

	      print "$pass\n";

	      print ADJUD "Positive graphite ajudication (95%) for $r[0]:$r[1] => Sufficient depth and strandedness   => proband:$pb_dpnfp:@pb_dp4nfp\n";


	      if ( ($pb_dp4nfp[2] < 3) || ($pb_dp4nfp[3] < 3 ) ) {
		print ADJUD "Negative graphite ajudication (95%) for $r[0]:$r[1] => Bad Strandedness  => proband:$pb_dpnfp:@pb_dp4nfp\n";
		next;
	      }

	      $ajudvars{$pass} = '1';  # add pos vars to ajud hash for later


	    } else {
	      print ADJUD "Negative ajudication (95%) for $r[0]:$r[1] => Sufficient depth   => proband:$pb_dpnfp:@pb_dp4nfp\n";
	    }

	  } else {
	      print ADJUD "Negative ajudication (95%) for $r[0]:$r[1] => Insufficient depth => proband:$pb_dpnfp:@pb_dp4nfp\n";
	  }

	}

	if ($iht eq 'd') {

	  my $pblim = 3;

	  if ( ($pb_dpnfp >= 10) ) {
	    if ( ($pb_dp4nfp[2] >= $pblim) && ($pb_dp4nfp[3] >= $pblim) ) {

	      my $pass = join (":", $r[0],$r[1]);
	      $ajudvars{$pass} = '1';

	      print ADJUD "Positive ajudication (95%) for $r[0]:$r[1] => Sufficient depth   => proband:$pb_dpnfp:@pb_dp4nfp\n";

	    } else {
	      print ADJUD "Negative ajudication (95%) for $r[0]:$r[1] => Sufficient depth, bad strandedness  => proband:$pb_dpnfp:@pb_dp4nfp\n";
	    }

	  } else {
	      print ADJUD "Negative ajudication (95%) for $r[0]:$r[1] => Insufficient total depth => proband:$pb_dpnfp:@pb_dp4nfp\n";
	  }

	}

      }

      #----------------End_Singleton_section-----------#


    }
  }

# Now test each of the alleles in the interesect output list to make sure
# both allele are pos. adjud if compound het. Then, write a new output list
# with the adjudicated variants only.

  my @genelist = split /\,/, $genelist;

  foreach my $g (@genelist) {
    my @lname = split /\//, $g;  # gene lists may be full path

    open (LISTIN, "<$name[0].best_genes.intersect");
    open (LISTOUT, ">$name[0].best_genes.final");

    my $grank = 0;

    while (<LISTIN>) {
      if (/^#/) {
	chomp $_;
	print LISTOUT "$_\tFinal_rank\tGraphite1|Graphite2\n";
      }
      my @l = split /\t/, $_;

     # must test both alleles,  2nd allele is na in recessive or de novo

      if ( (exists ( $ajudvars{$l[9]} )) && ( $l[11] eq 'na' ) )  {
	$grank++;
	chomp $_;
	my $r1 = `grep $l[9] Adjudication.out`;
	chomp $r1;
	my @adj1 = split /\s=>\s/, $r1;
	my $adj2 = 'na';
	print LISTOUT "$_\t$grank\t$adj1[2]|$adj2\n";
      }	elsif ( (exists ( $ajudvars{$l[9]} ))  &&  (exists ( $ajudvars{$l[11]} ))  ) {
	$grank++;
	chomp $_;
	my $r1 = `grep $l[9] Adjudication.out`;
	chomp $r1;
	my @adj1 = split /\s=>\s/, $r1;
	my $r2 = `grep $l[9] Adjudication.out`;
	chomp $r2;
	my @adj2 = split /\s=>\s/, $r2;
	print LISTOUT "$_\t$grank\t$adj1[2]|$adj2[2]\n";

     } else {
	next;
      }
    }
  }
  close (LISTIN);
  close (LISTOUT);
  close (ADJUD);
}

#-----------------------end_PARSE_GRAPHITE-------------------#


#------------------------sub_TIME--------------------------#

sub TIME {

  my $time = `date`;
  chomp $time;
  foreach my $t (@_) {
    print "$time\: $t\n";
    print LOG "$time\: $t\n";
  }

}

#-----------------------end_sub_TIME-----------------------#

#------------------------sub_ERROR-------------------------#
sub ERROR {
  print "\nERROR: $_[0]\n\n";
  exit;
}
#----------------------end_sub_ERROR-----------------------#

#------------------------sub_CLEANDIR----------------------#
sub CLEANDIR {

#  unlink ("$name[0].burden.out");
#  if (-e "best_genes") {
#    unlink ("best_genes");
#  }

#  unlink ("$name[0].scored_variants.out");
#  unlink ("$name[0].formatted.vvp");
#  if (-e "graphite_vars.in") {
#    unlink ("graphite_vars.in");
#  }
#  unlink ("Adjudicated.vcf");
#  unlink ("$name[0].recode.vcf");
  system ("mv $name[0].recode.vcf.gz $name[0].vcf.gz;");

#  unlink ("vvp.config");
#  unlink ("$name[0].log");

#  unlink ("graphite.recode.vcf graphite.recode_1.vcf; graphite_vars.in");
  unlink ("$name[0].burden.sorted.tmp");
  unlink ("$name[0].burden.tmp");
  unlink ("$name[0].burden.picked.tmp", "$name[0].burden.picked.tmp2");
  unlink ("$name[0].phevor.in");
  unlink ("phevor.error", "tmp1");
#  unlink ("$name[0].$listname.tmp");
 # system("rm $name[0].phevor2.tmp");


}
#----------------------end_sub_ERROR-----------------------#


#------------------------sub_PLOT--------------------------#
# This subroutine checks for the customized phevor plotting
# script and makes a small plot of the phevor2
# results.  A png file is created.
# Recommend Anaconda2-4.3.0-xxx-x86_64  python2.7
# If graphite was run, then the failing genes are removed
# before plotting.

sub PLOT {

  if ($graphite) {
    TIME("Removing genes that failed graphite from the phevor plot.");
    my @best = `cut -f3 $name[0].best_genes.final`;
    my @inter = `cut -f3 $name[0].best_genes.intersect`;
    chomp @best;
    chomp @inter;

    my %pass; my @failed;

    #Bad genes and high par genes from PCGC runs
    my @muc = qw(ADH7 AGAP6 AHNAK2 AKAP3 ARHGEF5 ATN1 ATXN1 AXDND1 C11orf40 CCDC66 CD177 CDIB CDON CLEC18B CLEC18C CPS1 DSPP FAM8A1 FOXD4L6 GGT2 GOLGA6L18 GOLGA8K HLA HLA-A HLA-B HLA-C HLA-DQB1 HLA-DRB1 HLA-DRB5 HLA-DPB1 HRC IGFN1 KIAA2018 KRT10 KRT2 KRTAP10 KRTAP10-1 KRTAP5 KRTAP5-5 L1TD1 MAFA MAGEF1 MAML2 MAP3K1 MAP3K4 METRN MST1 MUC1 MUC12 MUC13 MUC15 MUC16 MUC17 MUC2 MUC20 MUC22 MUC3A MUC4 MUC5B MUC6 MUC7 NBPF1 NBPF15 NBPF9 NCOA3 NOP9 OR11H1 OVGP1 PABPC3 PDE4DIP PLEKHG5 PLIN4 PRKCSH QRICH2 RIMBP3B RNF212 RPL1 OR51A7 SERPINb3 SERPINB SH3BGR SIRPA SPATA31A3 TBC1D3G TPTE TRAK1 URI1 VEZF1 WNK1 ZNF384 ZNF527);

    if ($no_muc) {
      push(@failed, @muc);
    }


    @pass{@best} = ();
    foreach my $bad (@inter) {
      push (@failed, $bad) unless exists $pass{$bad};
    }

    my %del = ();
    @del{@failed} = ();


    open (IN, "<$name[0].phevor2");
    open (OUT, ">$name[0].phevor2_cleaned");
    while (<IN>) {
      my @r = split /\t|\s+/, $_;
      print OUT unless exists ($del{$r[1]});
    }
    close (IN);
    close (OUT);

  }

  TIME("Creating the manhattan-styled plot of the best genes.");
  my $PLOTTER = `which vaast_manhattan_custom.py`;
  chomp $PLOTTER;
  my @PLOTTER = split /\//,$PLOTTER;
  unless ($PLOTTER[-1] eq 'vaast_manhattan_custom.py') {
    ERROR("\n\nvaast_manhattan_custom.py was not found!\nPlease put vaast_manhattan_custom.py in your path to use this option.\n This option also requires python2.7 and matplotlib in your path\n\n)");
    exit;
  }

  if (-e "$name[0].phevor2_cleaned") {
    system("tail -n +2 $name[0].phevor2_cleaned >$name[0].phevor2.tmp");
  } else {
    system("tail -n +2 $name[0].phevor2 >$name[0].phevor2.tmp");
  }
  my $best = `head -n1 $name[0].phevor2.tmp`;
  my @best = split /\s+|\t/, $best;

  my $newpos =`grep $best[1] $name[0].best_genes.final`;


  my @newpos = split /\t/, $newpos;
  my @mypos = split /:/, $newpos[9];
  my @info1=split /\|/, $newpos[14];
  my @info2=split /\|/, $newpos[15];


  if ($mypos[0] =~ /X/i) {
    $mypos[0] = 23;
  }
  if ($mypos[0] =~ /Y/i) {
    $mypos[0] = 24;
  }
  if ($mypos[0] =~ /MT/i) {
    $mypos[0] = 25;
  }

  my $lloc = 'center';

  if ($mypos[0] <= 8) {
    $lloc = 'left';
  } elsif (($mypos[0] > 8 ) && ($mypos[0] < 12)) {
    $lloc = 'center';
  } elsif ($mypos[0] >= 12) {
    $lloc = 'right';
  }


  sub AACONV {
    my $conv = '';
    ($conv = $_[0]) =~ s/Ala/A/ig;
    ($conv = $conv) =~ s/Cys/C/ig;
    ($conv = $conv) =~ s/Asp/D/ig;
    ($conv = $conv) =~ s/Glu/E/ig;
    ($conv = $conv) =~ s/Phe/F/ig;
    ($conv = $conv) =~ s/Gly/G/ig;
    ($conv = $conv) =~ s/His/H/ig;
    ($conv = $conv) =~ s/Ile/I/ig;
    ($conv = $conv) =~ s/Lys/K/ig;
    ($conv = $conv) =~ s/Leu/L/ig;
    ($conv = $conv) =~ s/Met/M/ig;
    ($conv = $conv) =~ s/Asn/N/ig;
    ($conv = $conv) =~ s/Pro/P/ig;
    ($conv = $conv) =~ s/Gln/Q/ig;
    ($conv = $conv) =~ s/Arg/R/ig;
    ($conv = $conv) =~ s/Ser/A/ig;
    ($conv = $conv) =~ s/Thr/T/ig;
    ($conv = $conv) =~ s/Val/V/ig;
    ($conv = $conv) =~ s/Trp/W/ig;
    ($conv = $conv) =~ s/Tyr/Y/ig;
    ($conv = $conv) =~ s/Ter/X/ig;
    return $conv;
  }

  my $predlabel = "";
  my $redlabel = "";
  my $aa1 = "";
  my $aa2 = "";
  my @var1 = ();
  my @var2 = ();

  if (defined $info2[3]) {        #Variants for plots

    if ($info2[3] =~ /na/i) {     #check if var2 is aa else take dna
      @var2 = split /:/, $info2[2];
      $var2[1] = 'na';
    } else {
      @var2 = split /:/, $info2[3];
    }

    if ($info1[3] =~ /na/i) {     # check var1 ...
      @var1 = split /:/, $info1[2];
    } else {
      @var1 = split /:/, $info1[3];
    }
  }


  if ( (!defined $var1[1]) ) {    #when no annotation
    $var1[1] = 'na';
  }
  if ( (!defined $var2[1]) ) {
    $var1[1] = 'na';
  }

  if ($iht =~ /r/) {
    if ($var1[1] =~ /^p/) {        #use 1 letter aa codes
      my $conv = AACONV($var1[1]);
      $aa1 = $conv;
    } else {
      $aa1 = $var1[1];
    }

    if ($var2[1] =~ /^p/) {

      my $conv = AACONV($var2[1]);
      $aa2 = $conv;
    } elsif (defined $var2[1]) {
      $aa2 = $var2[1];
    }

    if (($var2[1] =~ /na/ig) && ($iht eq 'r')) {
      $aa2 = $aa1;
      $redlabel = join("/", $aa1, $aa2);
    }  else {
      $redlabel = join("/", $aa1, $aa2);
    }
  }

  if ($iht =~ /d/) {
 #   if (($var2[1] =~ /na/) && ($var1[1] !~ /na/ig) && ($iht eq 'd')) {
      $redlabel = $var1[1];
#    }
  }

($predlabel = $redlabel) =~ s/\>/\\>/g;

  TIME("Running: $PLOTTER --phevor2 $name[0].phevor2.tmp $name[0]_phevor2.png  ~/bin/Phevor2/ontos/genes.gff3 --title $name[0] --mm_width 50 --mm_height 50 --red_genes $best[1] --rpos $lloc  --red_label $predlabel");

  system("$PLOTTER --phevor2 $name[0].phevor2.tmp $name[0]_phevor2.png  ~/bin/Phevor2/ontos/genes.gff3 --title $name[0] --mm_width 50 --mm_height 50 --red_genes $best[1]  --rpos $lloc --red_label $predlabel");

}

#----------------------end_sub_PLOT-------------------------#
