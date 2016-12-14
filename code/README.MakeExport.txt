globalement
il faut
ADAPTER le fichier Makefile.gen pour l'environnement
puis
CONSTRUCTION des lib

cd Algebres 
make 
cd Image
make
cd Tomo 
make 

ATTENTION il faut peut être placer d'abord les fichiers .h dans le
répertoire d'include (chez moi ~/include) .... (les binaires sont dans
~/bin)

puis utilisation des lib dans des executables
je fournis deux exemples d'utilisation de la classe GFBsinogramme
(dans Tomo) mais attention un des deux exemples GFBLineSimulData
donnent des résultats (étranges) qu'il me faudrait analyser

make GFBLineSimulData  
# usage GFBLineSimulData< entryGFBLineSimulData mais je recommande
# de lire d'abord entryGFBLineSimulData qui sont des réponses typiques
# aux questions et d'executer "entryGFBLineSimulData" en répondand à
# la main les premières utilisations
RQ : les codes ici utilises le phantom de Shepp et Logan  (10 ellipses
...(les deux zéro qui suivent sont 0 disque et 0 triangle.... les
triangles, c'est historique.... ce n'est pas supporté
actuellement.... mais ça pourrait)

make GFBCircularSimulData 
#usage : GFBCircularSimulData < entryGFBCircularSimulData


FIN du README

(ci dessous un mémo pour moi)
################################################################ 
################################################################ 
################################################################ 
ce que j'ai fait
#########################
# dans ~/Datas/Modules/Tomo
cd ~/Datas/Modules/Tomo

cp projection.C sinogramme.C FBDetLine.C GFBsinogramme.C FBsinogramme.C FBconvol.C pseudofct.C PsConvol.C integ.C  ../ExportPourThomas/Tomo/           

cp projection.h    sinogramme.h FBDetLine.h GFBsinogramme.h FBsinogramme.h FBconvol.h      PsConvol.h      pseudofct.h   integ.h  tomo.h tomotypes.h ../ExportPourThomas/Tomo/

########################
# dans ~/Datas/Modules/Image (on oublie deformation.o  deformation3D.o
# en adaptant le Makefile)
cd ~/Datas/Modules/Image

cp Makefile ../ExportPourThomas/Image/
cp realimage.h realimages3d.h realimage.C realimages3d.C

###########################
# algèbre qui contient entre autres des codes des numerical recipes
# (j'ai la licence)
cd ~/Datas/Modules
cp -a Algebres  ExportPourThomas/

# il y a des warning avec le make qui cesse(raie)nt avec l'option
#-w à  de g++  (mais je n'aime pas trop.... je l'ai juste testé dans
#le makefile des executables....)


###########################################
# executables pour test
# ça compile et ça semble s'executer a peu près mais il faut
# regarder...


# dans /Users/desbat/Datas/src/C++/FBI copie dans ExportPourThomas/src/  de deux exemple d'utilisation
cd /Users/desbat/Datas/src/C++/FBI
cp GFBCircularSimulData.C GFBLineSimulData.C ~/Datas/Modules/ExportPourThomas/src/

# et adaptation du Makefile
cp Makefile.thomas ~/Datas/Modules/ExportPourThomas/src/Makefile
# et adaptation du Makefile pour Thomas

# environnement de test
mkdir ~/Datas/Modules/ExportPourThomas/src/TEST

# cp du fichier de définition du phantom de Shepp et Logan
cp defSL ~/Datas/Modules/ExportPourThomas/src/TEST/
cp entryGFBLineSimulData ~/Datas/Modules/ExportPourThomas/src/TEST/
# usage : GFBLineSimulData < entryGFBLineSimulData 


# WARNING GFBCircularSimulData semble avoir des pb... ??? NaN ??? à voir....
# WARNING need numerical recipes for random noise
cp TEST/entryGFBCircularSimulData ~/Datas/Modules/ExportPourThomas/src/TEST/
  
usage : GFBCircularSimulData < entryGFBCircularSimulData


