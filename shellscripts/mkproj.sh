#!/bin/bash

########## <<< functions >>> ##########

function mkproj()
{
    if [ ! -d "$@"} ]; then
        printf "\n Making a directory named \"${@}\" \n"
        mkdir -p "$@"
    else
        printf "\n There is already a directory named \"${@}\" \n"
    fi
}


########## <<< config >>> ##########

DIR=(data log others summary scripts 01_qc/{before/multiqc_before,after/multiqc_after} 02_trim 03_align 04_filtering 05_peaks notebooks/{data,figures})
LINK=(01_qc/before/multiqc_before 01_qc/after/multiqc_after others scripts summary log notebooks)

########## Making directories ##########

printf "\n <<< Type project name (e.g. ChIP_Stra8) >>> \n\n"
read projname

mkproj ~/project; cd ~/project

if [ ! -d ${projname} ]
then
    printf "\n ---> Making a directory named \"$projname\" \n"
    printf "\n <<< ----- Initialization ----- >>> \n\n"
    mkdir -p ${projname}
else
    printf "\n ---> There is already a directory named \"$projname\"\n\n"
fi


cd ${projname}
for directory in "${DIR[@]}"; do
    mkproj ${directory}
done

# For loop proccesing
touch others/{sample.txt,rename.sh,rename.csv}

# For scripts
# git clone

########## change directory to make symbolic link ##########
mkproj REPORTS_${projname}; cd REPORTS_${projname}

for path in "${LINK[@]}"; do
    basename=`echo ${path} | xargs basename`
    ln -s ../${path} ${basename}
done

mkproj ~/PROJECT_REPORTS; cd ~/PROJECT_REPORTS; ln -s ~/project/${projname}/REPORTS_${projname} REPORTS_${projname}

if type tree
then
    tree ~/project/${projname}
else
    printf \
        """
        ${projname}
        ├── 01_qc
        │   ├── after
        │   │   └── multiqc_after
        │   └── before
        │       └── multiqc_before
        ├── 02_trim
        ├── 03_align
        ├── 04_filtering
        ├── 05_peaks
        ├── REPORTS_${projname}
        │   ├── log -> ../log
        │   ├── multiqc_after -> ../01_qc/after/multiqc_after
        │   ├── multiqc_before -> ../01_qc/before/multiqc_before
        │   ├── notebooks -> ../notebooks
        │   ├── others -> ../others
        │   ├── scripts -> ../scripts
        │   └── summary -> ../summary
        ├── data
        ├── log
        ├── notebooks
        │   ├── data
        │   └── figures
        ├── others
        │   ├── rename.csv
        │   ├── rename.sh
        │   └── sample.txt
        ├── scripts
        └── summary
        """
fi

echo