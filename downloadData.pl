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

# refGene.txt file name and chromosome name
my $outfile="${build}_refGene.txt.gz";
my $url="http://hgdownload.cse.ucsc.edu/goldenPath/$build/database/refGene.txt.gz";

my $outfile2="${build}_genome.fa";
my $url2="http://hgdownload.cse.ucsc.edu/goldenPath/$build/bigZips/chromFa.tar.gz";

# Download the file
info("Download refGene.txt.gz file...",0);
my $rc = getstore($url,$outfile);
if(is_success($rc)){
    print "OK\n";
}else{
    print "Failed \n";
    exit(1);
}

# unzip the file
`gunzip  $outfile`;

#
info("Download chromFa.tar.gz file...",0);
my $rc2 = getstore($url2,"chromFa.tar.gz");
if(is_success($rc2)){
    print "OK\n";
}else{
    print "Failed \n";
    exit(1);
}

# unzip the file
`tar -zxvf chromFa.tar.gz`;
`cat chr*.fa >$outfile2`;
`rm -rf chromFa.tar.gz`;



sub info{
    my ($str,$flag) = @_;
    print "[",scalar(localtime),"] $str";
    print "\n" if ($flag);
}
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
