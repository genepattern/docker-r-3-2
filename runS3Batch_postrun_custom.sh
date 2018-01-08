# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.

echo "R3.11 POST RUN CUSTOM: PERFORMING aws s3 sync $R_LIBS $S3_ROOT$R_LIBS_S3"
aws s3 sync $R_LIBS $S3_ROOT$R_LIBS_S3 --quiet

