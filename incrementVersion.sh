#!/bin/bash

## Short script for incrementing the version number in all of PopDel's files.
if [ "$#" -ne "1" ]
then
   echo "Usage: incrementVersion.sh major|minor|patch"
   exit 1
fi

##Path of the Makefile
makefile="Makefile"
##Path to README.md
readme="README.md"

#get todays date
today=$(date +%Y-%m-%d)
#Get current version from Makefile:
cVer=$(grep "^VERSION=" ${makefile} | sed -e 's;VERSION=;;')
#Get current version from Readme:
cVerRm=$(grep "\sPopDel version: " ${readme} | sed -e 's;\s*PopDel version: ;;')

## Check for inconsistencies and abort if any were found.
if [ "${cVer}" != "${cVerRm}" ]
then
   echo "Inconsistent versions between files!"
   echo -e "${makefile}:\t${cVer}"
   echo -e "${readme}:\t${cVerRm}"
   echo "Please fix versions before continuing."
   exit 1
fi

## Calculate the new version
if [ "$1" == "major" ]
then
   nVer=$(echo $cVer | awk -F '.' 'OFS="." {print $1+1, 0, 0}')
elif [ "$1" == "minor" ]
then
   nVer=$(echo $cVer | awk -F '.' 'OFS="." {print $1, $2+1, 0}')
elif [ "$1" == "patch" ]
then
   nVer=$(echo $cVer | awk -F '.' 'OFS="." {print $1, $2, $3+1}')
else
   echo "Please specify which position of the version string to update. Use one of \"major\", \"minor\", \"patch\" as first and only argument."
   exit 1
fi

## Ask for confirmation and apply the changes.
read -p "Update the version strings from ${cVer} to ${nVer} and the date to ${today}? [y/n]" -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
   sed -i -e "s;^VERSION=${cVer};VERSION=${nVer};" ${makefile}
   sed -i -e "s;^DATE=on .*;DATE=on ${today};" ${makefile}
   sed -i -e "s;PopDel version: ${cVer};PopDel version: ${nVer};" ${readme}
   sed -i -e "s;Last update: .*;Last update: ${today};" ${readme}
   echo "All done. Release the Kraken!"
   exit 0
else
   echo "Operation aborted by user."
   exit 1
fi
