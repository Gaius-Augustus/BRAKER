#!/usr/bin/env bash
# ==============================================================
# Tomas Bruna
#
# Compute accuracies of a prediction for a list of features
# ==============================================================


compare () {
    label="$1"
    flags="$2"
    acc="$($(dirname $0)/compare_intervals_exact.pl --pseudo $pseudo \
        --f1 "$annot" --f2 "$sortedPrediction" --$flags | cut -f4 | tail -3)"
    paste <(echo -e "${label}_Sn\n${label}_Sp") <(echo "$acc") -d"\t"
}

if  [ "$#" -lt 4 ]; then
    echo "Usage: $0 annot.gtf pseudogenes prediction.gtf feature_1 .. feature_n"
    echo "Example: $0 annot.gtf pseudogenes.gff prediction.gtf gene cds intron"
    exit
fi

annot=$1; shift
pseudo="$1"; shift
prediction="$1"; shift

sortedPrediction=$(mktemp)

sort -k1,1 -k4,4n -k5,5n $prediction > $sortedPrediction

for type in "$@"; do
    if [ $type == "start" ] || [ $type == "stop" ] || [ $type == "intron" ]; then
        flags="$type --no"
    else
        flags="$type"
    fi
    compare $type "$flags"
done

rm $sortedPrediction
