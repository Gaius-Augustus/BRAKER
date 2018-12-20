#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# gatech_pmp2hints.pl                                                                              #
# Convert Georgia Tech (Atlanta) intron information from protein mapping pipeline to braker.pl     #
# hints format.                                                                                    #
# Authors: Katharina Hoff                                                                          #
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
#                                                                                                  #
# This script is under the Artistic Licence                                                        #
# (http://www.opensource.org/licenses/artistic-license.php)                                        #
#                                                                                                  #                                  #
# Usage: gatech_pmp2hints.pl < in > out
####################################################################################################

use strict;
use warnings;

while(<STDIN>){
	my @t = split(/\t/);
	print "$t[0]\t$t[1]\tintron\t$t[3]\t$t[4]\t$t[5]\t$t[6]\t$t[7]\tsrc=P;mult=$t[5];\n";
}