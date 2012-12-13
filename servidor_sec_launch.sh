#!/bin/bash
#===============================================================================
#
#          FILE:  sec_launch_servidor.sh
# 
#         USAGE:  spatcontrol/sec_launch_servidor.sh myRscript.R
# 
#   DESCRIPTION:  send a command to computation server
# 
#       OPTIONS:  ---
#  REQUIREMENTS:  ---
#          BUGS:  ---
#         NOTES:  ---
#        AUTHOR:   (), 
#       COMPANY:  
#       VERSION:  1.0
#       CREATED:  12/08/2012 09:47:39 AM PET
#      REVISION:  ---
#===============================================================================

ipServidor="192.168.1.2"
folderExec="$PWD"
shortFolderExec="$(basename $PWD)"
fileExec="$1"

echo "Copiando \n $folderExec \n en $ipServidor:$shortFolderExec"
puja_a_servidor.sh $folderExec
echo "Terminado de copiar \n"

echo "Haciendo:
$ipServidor:$shortFolderExec/spatcontrol/sec_launch.sh $fileExec"

ssh -n -f $ipServidor "sh -c 'cd $shortFolderExec ; nohup spatcontrol/sec_launch.sh $fileExec 2>&1 &'"

# echo "Contenido de \"$USER@192.168.1.2:$shortFolderExec/sec_launch.log:"

sleep 2

echo "Puede checkar como va el calculo directamente en servidor o con:"
agara_output_servidor.sh
"

