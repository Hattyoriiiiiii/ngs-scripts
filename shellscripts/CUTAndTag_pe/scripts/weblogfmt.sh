#!/bin/bash
# usage: Rscript scripts/log_summary.R && sh scripts/weblogfmt.sh summary/log_summary.txt > weblog.html

function tagit()
{
    printf '<%s>%s</%s>\n' "${1}" "${2}" "${1}"
}

# basic header tags
echo "<!DOCTYPE html>"
echo "<html>"
echo "<body>"
# echo "<h1>$1</h1>"  # タイトル
echo "<h1>web summary</h1>"  # タイトル

echo "<table border=1>"
echo "<tr>"  # new table row
echo "<th>Cell</th>"  # column header
echo "<th>Condition</th>"
echo "<th>Replicate</th>"
echo "<th>Experiment</th>"
echo "<th>Sequencing Depth</th>"
echo "<th>Mapped Fragment数</th>"
echo "<th>Alignment Rate (%)</th>"
echo "<th>Alignment Uniquely (%)</th>"

echo "</tr>"


while read f1 f2 f3 f4 f5 f6 f7 f8
do
    echo "<tr>"

    tagit "td" "${f1}"
    tagit "td" "${f2}"
    tagit "td" "${f3}"
    tagit "td" "${f4}"
    tagit "td" "${f5}"
    tagit "td" "${f6}"
    tagit "td" "${f7}"
    tagit "td" "${f8}"

    echo "</tr>"
done < $1


# close tags
echo "</table>"
echo "</body>"
echo "</html>"