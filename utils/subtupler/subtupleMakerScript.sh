#!/bin/sh

#Runs Takes in an array of files and runs the subtupler across them.  Outputs
#The files to the defined output directory
while test $# -gt 0; do
    case "$1" in
        --help)
            echo "-datadir, specify data directory to subtuple"
            echo "-outdir, specifiy output directory for subtuples"
            exit 0
            ;;
        -d)
            shift
            if test $# -gt 0; then
                export DATADIRECTORY=$1
            else
                echo "define your data directory"
                exit 1 
            fi
            shift
            ;;
        -t)
            shift
            if test $# -gt 0; then
                export DATTYPE=$1
            else
                echo "Define the rootfile type (DATA or MC)"
                exit 1
                #mkdir output
                #OUTDIRECTORY=./output/
            fi
            shift
            ;;
        -o)
            shift
            if test $# -gt 0; then
                export OUTDIRECTORY=$1
            else
                echo "define your output directory"
                exit 1
                #mkdir output
                #OUTDIRECTORY=./output/
            fi
            shift
            ;;
    esac
done 

FILES=( ${DATADIRECTORY}* )

for f in ${FILES[@]};
do
    echo $f
    fstrip=${f##$DATADIRECTORY}
    echo $fstrip
    subtup_strip=${fstrip//.ntuple.root/}_subtup.ntuple.root
    echo $subtup_strip
    echo "DATTYPE: "
    echo ${DATTYPE}
    #if ((${DATTYPE}=="DATA"))
    #then
    #  ./subtupleMaker_AmBe --data $f --out ${OUTDIRECTORY}$subtup_strip
    #fi

    if ((${DATTYPE}=="MC"))
    then
      ./subtupleMaker_AmBe --MC $f --out ${OUTDIRECTORY}$subtup_strip
    fi
done;
