#! /bin/sh

#==========#
#  INPUTS  #
#==========#
# tag dans cvs de la version a distribuer
TAGORIG=$1
# path of the place where the created distribution will be put (full path)
PATHBUILD=$2"/openfem"
#=====================#
#  PATH TO CUSTOMIZE  #
#=====================#  
# make_distrib path [CUSTOM]
PATHMAKE="/home/domi/tmp/matlab_distrib"


PATHDEPOT=$2"/openfem_depot"

if [ $# != 2 ]; then
    echo "$0 fromdir todir"
    echo "$0 : bad args"
    exit 2
fi

if [ ! -d $2 ]; then
    mkdir $2
fi

if [ -d $PATHBUILD ]; then
    echo "openfem directory exists. Specify an other path to build the distribution."
    exit
fi

if [ ! -d $PATHBUILD ]; then
    mkdir $PATHBUILD
fi

if [ -d $PATHDEPOT ]; then
    echo "openfem_depot directory exists. Specify an other path to build the distribution."
   exit
fi

if [ ! -d $PATHDEPOT ]; then
    mkdir $PATHDEPOT
fi

#==================================================#
#  retrait de depot cvs de la version a distribuer #
#==================================================#
cd $PATHDEPOT
CVSOPEN="export -r "$1" openfem"
CVSSRC="export -r"$1" src"
CVSTEX="export -r"$1" tex"
CVSHTML="export -r"$1" html"
cvs $CVSOPEN
cvs $CVSSRC
cvs $CVSTEX
cvs $CVSHTML
cd ../

PATHORIG=$PATHDEPOT"/openfem"

#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#
# DON'T MODIFY WHAT IS FOLLOWING !                                                            #
#---------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------#                 
PATHCUR=`pwd`

#====================================================#
#             directories creation                   #
#====================================================#
DEMOS=$PATHBUILD"/demos"
HTML=$PATHBUILD"/html"
MEX=$PATHBUILD"/mex"
OFACT=$PATHBUILD"/@ofact"
SDT3=$PATHBUILD"/sdt3"
TEST=$PATHBUILD"/test"
TEX=$PATHBUILD"/tex"
SRC=$PATHBUILD"/src"
VISU1=$DEMOS"/visu"
VISU2=$TEST"/visu"
PLOTS=$TEX"/plots"

mkdir $DEMOS
mkdir $HTML
mkdir $MEX
mkdir $OFACT
mkdir $SDT3
mkdir $TEST
mkdir $TEX
mkdir $SRC
mkdir $VISU1
mkdir $VISU2
mkdir $PLOTS

#====================================================#
#              principal directory                   #
#====================================================#
CURDIR=$PATHORIG"/*.m"
LGPL=$PATHORIG"/LGPL.txt"
INST=$PATHORIG"/INSTALL.txt"
MARQ=$PATHORIG"/MARQUE.html"
TRAD=$PATHORIG"/TRADEMARK.html"
#OFUT=$PATHMAKE"/util/ofutil.m"

CURCOP=$PATHBUILD"/."

cp $CURDIR $CURCOP
cp $LGPL $CURCOP
cp $INST $CURCOP
cp $MARQ $CURCOP
cp $TRAD $CURCOP
#cp $OFUT $CURCOP

#=========================================#
#            demos directory              #
#=========================================#
DEMOSDIR=$PATHORIG"/demos/*"
DEMOSCOP=$DEMOS"/."
cp $DEMOSDIR $DEMOSCOP

#=========================================#
#             html directory              #
#=========================================#
HTMLDIR1=$PATHORIG"/../html/*.html"
HTMLDIR2=$PATHORIG"/../html/*.gif"
HTMLCOP=$HTML"/."
cp $HTMLDIR1 $HTMLCOP
cp $HTMLDIR2 $HTMLCOP

#=========================================#
#              mex directory              #
#=========================================#
MEXDIR=$PATHORIG"/mex/*.c"
MEXCOP=$MEX"/."
cp $MEXDIR $MEXCOP

#=========================================#
#           ofact directory               #
#=========================================#
OFACTDIR=$PATHORIG"/@ofact/*.m"
OFACTCOP=$OFACT"/."
cp $OFACTDIR $OFACTCOP

#=========================================#
#            sdt3 directory               #
#=========================================#
SDT3DIR=$PATHORIG"/sdt3/*.m"
SDT3COP=$SDT3"/."
cp $SDT3DIR $SDT3COP

#=========================================#
#            test directory               #
#=========================================#
TESTDIR=$PATHORIG"/test/*"
TESTCOP=$TEST"/."
cp $TESTDIR $TESTCOP

#=========================================#
#             tex directory               #
#=========================================#
TEXDIR=$PATHORIG"/../tex/*"
PLOTSDIR=$PATHORIG"/../tex/plots/*"
TEXCOP=$TEX"/."
PLOTSCOP=$PLOTS"/."
cp $TEXDIR $TEXCOP
cp $PLOTSDIR $PLOTSCOP

#=========================================#
#             src directory               #
#=========================================#
SRCDIR1=$PATHORIG"/../src/*.f"
SRCDIR2=$PATHORIG"/../src/*.c"
SRCDIR3=$PATHORIG"/../src/*.ins"
SRCDIR4=$PATHORIG"/../src/f2c.h"
SRCCOP=$SRC"/."

cp $SRCDIR1 $SRCCOP
cp $SRCDIR2 $SRCCOP
cp $SRCDIR3 $SRCCOP
cp $SRCDIR4 $SRCCOP