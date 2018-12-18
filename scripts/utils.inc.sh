## utils.inc.sh
##
## Copyright (c) 2017 Institut Curie                               
## Author(s): Eric Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

###########################
## trap handler
###########################
function trap_error()
{   
    echo "Error: $1 - line $2 - exit status of last command: $?. Exit" >&2
    exit 1
}

function trap_exit()
{
    ##Since bash-4.0 $LINENO is reset to 1 when the trap is triggered
    if [ "$?" != "0" ]; then
	echo "Error: exit status detected. Exit." >&2
	touch ${ODIR}/exit-err
    fi

    if [ -e ${ODIR}/mapping/tmp ]; then 
	echo -e "Cleaning temporary folders ..." >&2
	/bin/rm -rf ${ODIR}/mapping/tmp; 
    fi
}

trap 'trap_error "$0" "$LINENO"' ERR
trap 'trap_exit' 0 1 2 3

set -E ## export trap to functions
set -o pipefail  ## trace ERR through pipes         
#set -o errexit   ## set -e : exit the script if any statement returns a non-true return value


###########################
## Subroutine for pipelines
###########################

die() 
{ 
    echo "Exit: $@" >&2
    exit 1
}

exec_cmd()
{
    echo $*
    if [ -z "${DRY_RUN+x}" ]; then
	eval "$@" || die 'Error'
    fi
}

exec_ret()
{
    if [ -z "${DRY_RUN+x}" ]; then
	eval "$@" || die 'Error'
    fi
}

abspath() 
{
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

read_config()
{
    eval "$(sed -e '/^$/d' -e '/^#/d' -e 's/ =/=/' -e 's/= /=/' -e 's/ *$//' $1 | \
awk -F"=" '{printf("%s=\"%s\"; export %s;\n", $1, $2, $1)} $1~"PATH" && $2!=""{printf("export PATH=%s:$PATH;\n", $2)}')"
}

is_in_path()
{
    type -P $1 > /dev/null && echo 1 || echo 0
}


file_exists()
{

    if [ ! -e $1 ]; then
 	echo -e "Error: The file ${1} was not found. Exit" >&2
 	echo
 	exit 1
    elif [ ! -s $1 ]; then
 	echo -e "Error: The file ${1} was found but is empty. Exit" >&2
 	echo
 	exit 1
    fi
}

get_fastq_prefix()
{
    basename $1 | sed -e 's/[\._]*R*[12]*.fastq\(.gz\)*//'
}

check_output_files()
{
    local odir=$1
    local ofile=$2

    ## Must remove last folder of ${odir}
    odir=$(dirname $odir)
    
    if [ -s ${ofile} ]; then
	while read f; do
	    if [ -e ${odir}/${f} ]; then
		continue
	    else
		die "Check failed on file ${odir}/${f}"
	    fi
	done < ${ofile}
    else
	die "Input file '$ofile' not found !"
    fi
}
