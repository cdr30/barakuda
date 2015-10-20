#!/bin/sh


list=`\ls *.html`

for ff in ${list}; do
    cat ${ff} | sed -e s/"MEDUZA"/"BARAKUDA"/g > out
    mv -f out ${ff}
done

