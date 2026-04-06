#!/bin/sh
#################################
# RM.out to gff3 with Colors!!! #
#################################
# Author: Clément Goubert
# Modified: 2025-08-14 — treat attributes as a single string, filter lines starting with number
#################################

if [ -z "$1" ]; then
    echo "No input provided (RepeatMasker .out file)"
    echo "USAGE: ./rm2gff3.sh RM.out"
    exit 1
fi

awk -v OFS="\t" \
    -v LINEcol="#3399ff" \
    -v SINEcol="#800080" \
    -v DNAcol="#ff6666" \
    -v LTRcol="#00cc44" \
    -v RCcol="#ff6600" \
    -v Low_complexitycol="#d1d1e0" \
    -v Satellitecol="#ff99ff" \
    -v Simple_repeatcol="#8686ac" \
    -v Unknowncol="#f2f2f2" '
BEGIN {
    print "##gff-version 3"
}
$1 ~ /^[0-9]/ {
    sw  = $1;  div = $2;  del = $3;  ins = $4;
    seq = $5;  qs  = $6;  qe  = $7;  ql  = $8;
    str = ($9 == "C" ? "-" : $9);
    rep = $10; class = $11; rb = $12; re = $13; rl = $14; id = $15;

    if (str == "+") {
        target = "Target=" rep " " class " " rb " " re " " rl;
    } else {
        target = "Target=" rep " " class " " re " " rb " " rl;
    }

    attrs = target ";Div=" div ";Del=" del ";Ins=" ins ";SWscore=" sw;

    color = Unknowncol;
    if (class ~ /^LINE\//)             color = LINEcol;
    else if (class ~ /^SINE\//)        color = SINEcol;
    else if (class ~ /^DNA\//)         color = DNAcol;
    else if (class ~ /^LTR\//)         color = LTRcol;
    else if (class ~ /^RC\//)          color = RCcol;
    else if (class ~ /Low_complexity/) color = Low_complexitycol;
    else if (class ~ /Satellite/)      color = Satellitecol;
    else if (class ~ /Simple_repeat/)  color = Simple_repeatcol;

    attrs = attrs ";color=" color;

    print seq, "RepeatMasker-4.0.1", "similarity", qs, qe, sw, str, ".", attrs;
}
' "$1"
