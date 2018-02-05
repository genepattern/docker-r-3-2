# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.

: ${R_LIBS_S3=/genepattern-server/Rlibraries/R303/rlibs}
export R_LIBS=/usr/local/lib/R/site-library

##################################################
# MODIFICATION FOR R PACKAGE INSTALLATION
##################################################
# mount pre-compiled libs from S3
echo "FOR R3. 13 CUSTOMIZING: aws s3 sync $S3_ROOT$R_LIBS_S3 $R_LIBS --quiet"

ls $R_LIBS


if [ -f "r.package.info" ]
then
        echo "$TASKLIB/r.package.info found."
        Rscript /build/source/installPackages.R r.package.info
else
        echo "$TASKLIB/r.package.info not found."
fi



