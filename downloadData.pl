#!/usr/bin/env perl
###################################
# Author: Jiang Li
# Email: riverlee2008@gmail.com
# Date: Sat Jul  6 14:13:00 2013
###################################
use strict;
use warnings;
use Getopt::Long;
use LWP::Simple;


# Download the refGene table from UCSC
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
my ($build,$help) = ("hg19",0);

# Read parameters from the command line
 unless(GetOptions("build=s"=>\$build,"h|help"=>\$help)){
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

my $outfile="${build}_refGene.txt.gz";
my $url="http://hgdownload.cse.ucsc.edu/goldenPath/$build/database/refGene.txt.gz";

# Download the file
my $rc = getstore($url,$outfile);

# unzip the file
`gunzip  $outfile`;




sub help{
    my ($flag) = @_;
    print <<USAGE;
Usage: downloadData.pl -build <hg19> -h
      -build either hg19,hg18 or hg17, default is hg19
      -h     display this help
########################################################
# This program will download refGene.txt and genome refernce
# sequences from UCSC and store them in the  data folder
########################################################
USAGE
    exit(1) if($flag);
}
