#!/bin/bash
#===============================================================================
#
#          FILE:  test_extrapol_field.sh
# 
#         USAGE:  ./test_extrapol_field.sh 
# 
#   DESCRIPTION:  test that the outputs are identical from the current extrapol_field
# 		and the recorded test
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   (), 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  02/27/2012 11:27:36 AM EST
#      REVISION:  ---
#===============================================================================

out=$(R CMD BATCH test_extrapol_field.r output.lst)
echo $out

diff=$(diff usamples.txt usamplestest.txt)
if [[ $diff != "" ]] ; then 
	echo "diff in u"
fi
diff=$(diff wsamples.txt wsamplestest.txt)
if [[ $diff != "" ]] ; then 
	echo "diff in w"
fi

