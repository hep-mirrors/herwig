#
export DIR=/theory/Herwig++/
export REPOSITORY=$DIR/Repository
export PROJECT=/var/pcce/usera/gieseke/work/mc/Herwig++
#
# --- First: create Repository
#
cd $DIR
mkdir Repository/
cvs -d $REPOSITORY init
#
# --- Second: import the project in the Repository
#
cd $PROJECT
cvs -d $REPOSITORY import -m "Imported Herwig++" Herwig++ yoyo start
#
# --- Third: rename project as "old", and use instead the
#            one from the Repository.
#
cd ../
mv Herwig++ old.Herwig++
#
cvs -d $REPOSITORY checkout Herwig++
#
