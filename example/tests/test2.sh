wd=test2

if [ -d $wd ]; then
    rm -r $wd
fi

# Attention:
#
# We are in the process of changing protein handling within BRAKER.
# The files ../prothint_augustus_reg.gff, ../evidence_augustus.gff and ../top_chains.gff are currently not
# standard output files of ProtHint.
# The following manipulations have been performed:
#
# ../evidence_augustus.gff: src=P has been replaced by src=C
# ../prothint_augustus_reg.gff: original prothint_augustus.gff was generated with OrthoDB Arthropoda excluding order
#                               of Drosophila; multiplicity of hints was normalized by median and subsequently, logistic
#                               regression was performed to classify the hints into two classes: "less reliable" = 0, and
#                               "more reliable" = 2.
#                               my $y = -4.00529 + 4.73909 * $normalized_mults[$i] + 9.09026 * $al_score;
#                               my $class_label = 0;
#                               if($y >= 0.85) {
#                                   $class_label = 2;
#                               }
#
# Be aware that this test uses a custom exstrinsic.cfg file at the moment.

( time braker.pl --genome=../genome.fa --hints=../prothint_augustus_reg.gff,../evidence_augustus.gff,../top_chains.gff --evidence=../evidence.gff --prothints=../prothint.gff --epmode --softmasking --workingdir=$wd --extrinsicCfgFiles=../../scripts/cfg/top_chain_order.cfg ) &> test2.log
