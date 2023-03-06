#!/bin/bash

CUTOFF=10
GFF_PATH=false
BAM_PATH=false

show_help() {
cat << EOF

Usage: ${0##*/} [-g -b -c... -h ] [-g path/to/file ] [-c int]...
This tool uses Samtools coverage, to determine which gff features were present in any given dataset
based on BAM alignment. 

	-g			Gff-file input [path/to/file]
	-b			Bam file input [path/to/file]
	-c			Coverage threshold, regions below this number will be ignored
	--genesonly		Only look at gene features, ignore all other features
	-h | --help		Display help

EOF
}                

##########################################################################################################
# Command line argument parsing is done here mostly using binary switches
##########################################################################################################
while :; do
	case "$1" in
	-g)
		# If there is an argument given by the user, this will show up as $2
		# We double check that this argument is an actual argument and not just another
		# option by checking if there is no dash at the start of the argument
		# '-option' is ignored while 'option' is considered a super provided argument
		if [ "$2" ] && [[ ! $2 == -* ]]; then
			if [ -f "$2" ]; then
			  ### Take action if file  exists
				echo "Gff-file provided by user: ${2}";
				GFF_PATH=$(readlink -e ${2});
			else
			  ###  Control will jump here if file does NOT exists ###
				echo "Gff-file not supplied or does not exist";
				exit 1;
			fi
		# If an argument is passed, an additional shift is needed so the script doesn't then
		# read the argument again as an option
		shift
		fi
		;;
	-b)
		if [ "$2" ] && [[ ! $2 == -* ]]; then
			if [ -f "$2" ]; then
				echo "Bam file provided by user: ${2}";
				BAM_PATH=$(readlink -e ${2})
			else
				echo "Bam file not supplied or does not exist";
				exit 1;
			fi
		shift
		fi
		;;
	-c)
		if [ "$2" ] && [[ ! $2 == -* ]]; then
			if [[ $2 =~ ^-?[0-9]+$ ]] ; then
				CUTOFF="$2"
				echo "Coverage cutoff value is set to: $CUTOFF";
			else
				echo "overage cutoff provided by user is non numerical, default is used [10]";
			fi
			shift
		else
			continue
		fi
		;;
	--genesonly)
		GENESONLY_SWITCH=true
		;;
	-h|--help) 
		show_help >&2
		;;
	*)
		break
	esac
	shift
done

### Check if user actually provided files
if [[ "$GFF_PATH" = false ]] || [[ "$BAM_PATH" = false ]] ;  then
	echo "Please provide a Gff file of the reference, and a BAM file to query using the options -g and -b";
	exit 1;
fi

#Check if the bam file has an index
if [ ! -f "${BAM_PATH}.bai" ]; then
	samtools index ${BAM_PATH}
fi

#Create an empty output file
echo -n '' > Features_in_dataset.tsv

### Now go through each line of the GFF, if option -genesonly is used, it will ignore everything except lines containing 'gene'
while read FEATURE; do
	if [[ "$FEATURE" != \#* ]];then
		IFS=$'\t' read -r -a FEATURE_SPLIT <<< "$FEATURE"
		if [ "$GENESONLY_SWITCH" = true ]; then
			if [ ! "${FEATURE_SPLIT[2]}" = "gene" ]; then
				continue
			fi
		fi 
		RANGE="${FEATURE_SPLIT[0]}:${FEATURE_SPLIT[3]}-${FEATURE_SPLIT[4]}"
		SAM_OUTPUT="$(samtools coverage -r ${RANGE} $BAM_PATH)"
		IFS=$'\t'
		COVERAGE_VALUES=($SAM_OUTPUT)
		if [ "$(printf '%s\n' "${COVERAGE_VALUES[13]}" "${CUTOFF}" | sort -V | head -n1)" = "${CUTOFF}" ]; then 
			FEATURE="${FEATURE%%[[:space:]]}	;coverage=${COVERAGE_VALUES[13]}" 
			echo "${FEATURE}" >> Features_in_dataset.tsv
		fi
	fi
done < ${GFF_PATH}

