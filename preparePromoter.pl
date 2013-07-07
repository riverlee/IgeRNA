#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Sat Jul  6 16:14:13 2013
###################################
use strict;
use warnings;
use Getopt::Long;
use LWP::Simple;
use Cwd;

my $currentdir=getcwd;
##############################
# Get the script folder
$0 =~ /^(.*)\//;
my $homeDir = $1;
if (!defined($homeDir) || $homeDir eq '') {
    $homeDir = ".";
}
$homeDir .= "/";
my $pwd = '';
unless ($homeDir =~ /^\//) {
    #`pwd > .ls`;
    open IN, "pwd|";
    while (<IN>) {
        chomp;
        $pwd = $_;
        last;
    }
    close IN;
    $homeDir = $pwd . "/" . $homeDir;
}
chdir($homeDir) or die $!;
##################################

##################################
# Prepare promoter sequences based on the refGene.txt 
# which could be downloaded by downloadData.pl script
##################################
my ($build,$help,$upstreamExtend,$downstreamExtend) = ("hg19",0,2000,500);

# Read parameters from the command line
 unless(GetOptions("build=s"=>\$build,"h|help"=>\$help,"up=i"=>\$upstreamExtend,"down=i"=>\$downstreamExtend)){
    help();
    die $!;
}

# Check parameters
help(1) if ($help);
if($build ne "hg19" && $build ne "hg18" && $build ne "hg17"){
    help(1);
}

# Main Program to download the data
my $dirname="data";
mkdir $dirname unless ( -d $dirname);
chdir $dirname;

# Get the refGene.txt file
my $refgene="${build}_refGene.txt";
my $refgenome="${build}_genome.fa";

if( ! -e $refgene){
    die  "$refgene doesn't exists, please download it by downloadData.pl -build $build \n";
}
if( ! -e $refgenome){
    die "$refgenome doesn't exists, please download it by downloadData.pl -build $build \n";
}


##################################
# Main program
info("Loading refGene to get TSS ...",1);
my %refgene;
loadRefGene($refgene,\%refgene);
info("There are ".scalar(keys %refgene)." Genes",1);

#refineRefGene(\%refgene);

info("Loading genome sequeneces ...",1);
my %genome; 
loadGenome($refgenome,\%genome);
info("There are ".scalar(keys %genome)," Chromosomes",1);
#
info("Output promoter sequences ...",1);
my $output="${build}_promoter_upstream${upstreamExtend}_downstream${downstreamExtend}.fa";
chdir($currentdir) or die $!;
getPromoter(\%refgene,\%genome,$upstreamExtend,$downstreamExtend,$output);
info("Finished",1);

sub getPromoter{
    my ($ref,$genomeref,$upstreamExtend,$downstreamExtend,$output) = @_;
    open(OUT,">$output") or die $!;

    # Some genes may have multiple transcripts,
    # However, they may share the same TSS or similar TSS, thus 
    # the promoter may not changes too much
    foreach my $gene (sort keys %{$ref}){
        my @NMIDs = sort keys %{$ref->{$gene}};
        if (scalar(@NMIDs) ==1 ){
                my $nmid= $NMIDs[0];
                my $strand = $ref->{$gene}->{$nmid}->{strand};
                my $tss = $ref->{$gene}->{$nmid}->{TSS};
                my $chr = $ref->{$gene}->{$nmid}->{chr};

                my $start;
                my $seq="";
                if($strand eq "+"){
                    $start = $tss - $upstreamExtend;
                    $start = 0 if ($start <1);
                    $seq = substr($genomeref->{$chr},$start,$upstreamExtend+$downstreamExtend);
                }else{
                    $start = $tss - $downstreamExtend;
                    $start = 0 if ($start <1);
                    $seq = substr($genomeref->{$chr},$start,$upstreamExtend+$downstreamExtend);
                    # DO complemetary
                    $seq = reverse_complement($seq); 
                }

                # for the sequence name
                my $end = $start+$upstreamExtend+$downstreamExtend;
                $seq = formatSeq($seq,50);
                my $header = "$gene|$nmid";
                $header.="|$strand|TSS=$chr:$tss|Range=$chr:$start-$end";
                print OUT ">$header\n$seq\n";
                next;
        }
        my %tmp;
        #my %keep;
        foreach my $nm (@NMIDs){
            my $strand = $ref->{$gene}->{$nm}->{strand};
            my $chr = $ref->{$gene}->{$nm}->{chr};
            my $tss = $ref->{$gene}->{$nm}->{TSS};
            push @{$tmp{$chr}->{$strand}->{$tss}},$nm;
        }

        # Get the 5 most left TSS
        foreach my $chr (sort keys %tmp){
            foreach my $strand ( sort keys %{$tmp{$chr}}){
                my @tss;
                if ($strand eq "+"){
                    @tss = sort {$a <=> $b} keys %{$tmp{$chr}->{$strand}};
                }else{
                    @tss = sort {$b <=> $a} keys %{$tmp{$chr}->{$strand}};
                }
                my $tss = shift @tss;
                
                my @nms = @{$tmp{$chr}->{$strand}->{$tss}};

                # Get sequences
                my $seq="";
                my $start;
               # next if (!exists($genomeref->{$chr}));
                if($strand eq "+"){
                    $start = $tss - $upstreamExtend;
                    $start = 0 if ($start <1);
                    $seq = substr($genomeref->{$chr},$start,$upstreamExtend+$downstreamExtend);
                }else{
                    $start = $tss - $downstreamExtend;
                    $start = 0 if ($start <1);
                    $seq = substr($genomeref->{$chr},$start,$upstreamExtend+$downstreamExtend);
                    # DO complemetary
                    $seq = reverse_complement($seq); 
                }

                # for the sequence name
                my $end = $start+$upstreamExtend+$downstreamExtend;
                $seq = formatSeq($seq,50);
                my $header = "$gene|";
                $header.= join ";",@nms;
                $header.="|$strand|TSS=$chr:$tss|Range=$chr:$start-$end";
                print OUT ">$header\n$seq\n";
            }
        }
    }
    close OUT;
}

sub formatSeq{
    my ($seq,$len) = @_;
    
    # default to 50 characters of the sequences per line
    $len = 50 unless $len;
    my $formatted_seq="";
    while( my $chunk = substr($seq,0,$len,"")){
        $formatted_seq.="$chunk\n";
    }
    chomp($formatted_seq);
    return $formatted_seq;
}

sub reverse_complement{
    my $dna = shift;

    # Reverse the dna sequence
    my $revcomp =reverse($dna);

    # Complement the reversed DNA sequence
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

sub refineRefGene{
    # Some genes may have multiple transcripts,
    # However, they may share the same TSS or similar TSS, thus 
    # the promoter may not changes too much
    my ($ref) = @_;
    foreach my $gene (sort keys %{$ref}){
        my @NMIDs = sort keys %{$ref->{$gene}};
        next if (scalar(@NMIDs) ==1 );
        my %tmp;
        my %keep;
        foreach my $nm (@NMIDs){
            my $strand = $ref->{$gene}->{$nm}->{strand};
            my $chr = $ref->{$gene}->{$nm}->{chr};
            my $tss = $ref->{$gene}->{$nm}->{TSS};
            push @{$tmp{$chr}->{$strand}->{$tss}},$nm;
        }

        # Get the 5 most left TSS
        foreach my $chr (sort keys %tmp){
            foreach my $strand ( sort keys %{$tmp{$chr}}){
                my @tss;
                if ($strand eq "+"){
                    @tss = sort {$a <=> $b} keys %{$tmp{$chr}->{$strand}};
                }else{
                    @tss = sort {$b <=> $a} keys %{$tmp{$chr}->{$strand}};
                }
                my $tss = shift @tss;
                
                foreach my $nm (@{$tmp{$chr}->{$strand}->{$tss}}){
                    $keep{$nm}=1;
                }
            }
        }
        foreach my $nm(@NMIDs){
            unless (exists($keep{$nm})){
                delete $ref->{$gene}->{$nm};
            }
        }
    }
}

sub loadRefGene{
    my($in,$ref) = @_;
    open(IN,$in) or die $!;
    while(<IN>){
        s/\r|\n//g;
        # $name is transcript id (NM_) while $name2 the gene symbol
        my($bin,$name,$chr,$strand,$txStart,$txEnd,$cdsStart,$cdsEnd,$exoncount,$exonStarts,$exonEnds,$score,$name2,$cdsStartstat,$cdsendstat,$Exonframes) = split "\t";
       if($chr=~/chr[0-9XY]+$/){
           if($strand eq "+"){
               $ref->{$name2}->{$name}->{TSS}=$txStart;
           }else{
               $ref->{$name2}->{$name}->{TSS}=$txEnd;
           }
           $ref->{$name2}->{$name}->{strand}=$strand;
           $ref->{$name2}->{$name}->{chr}=$chr;
       }
    }
    close IN;
}

sub loadGenome{
    my ($in,$ref) = @_;
    open(IN,$in) or die $!;
    my $chr="";
    my $seq="";
    while(<IN>){
        s/\r|\n//g;
        if(/^>/){
            s/>//g;
            # if this is the first sequence name, just get the name of $chr
            if($chr eq ""){
                $chr=$_;
                next;
            }else{
                if($chr=~/chr[0-9XY]+$/){
                    $ref->{$chr}=$seq;
                }
                $chr=$_;
                $seq="";
                next;
            }
        }else{
            $seq.=$_;
        }
    }
    # Deal with the last sequence
    $ref->{$chr}=$seq if($chr=~/chr[0-9XY]+$/);
    close IN;
}





sub info{
    my ($str,$flag) = @_;
    print "[",scalar(localtime),"] $str";
    print "\n" if ($flag);
}
sub help{
    my ($flag) = @_;
    print <<USAGE;
Usage: preparePromoter.pl -build [hg19] -up [2000] -down [500] -h
      -build either hg19,hg18 or hg17, default is hg19
      -up    upstream of TSS, default is 2000bp
      -down  downstream of TSS default is 500bp
      -h     display this help
########################################################
# Prepare the promoter sequences in fasta format
########################################################
USAGE
    exit(1) if($flag);
}
