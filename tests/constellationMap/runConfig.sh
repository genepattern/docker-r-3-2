#!/bin/sh

TASKLIB=$PWD/src
INPUT_FILE_DIRECTORIES=$PWD/data
S3_ROOT=s3://moduleiotest
WORKING_DIR=$PWD/job_3113
RLIB=$PWD/rlib

COMMAND_LINE="Rscript $TASKLIB/run_ConstellationMap.R  $TASKLIB/ $TASKLIB/ $RLIB java   --input.gct.file=$INPUT_FILE_DIRECTORIES/GSE27473_series_matrix_collapsed_to_symbols.PROJ.gct --input.cls.file=$INPUT_FILE_DIRECTORIES/WT_vs_Knockdown.cls --gene.sets.database=c2.all.v5.1.symbols.gmt  --top.n=100  --direction=positive  --image.format=PDF  --jaccard.threshold=0.1 "

#COMMAND_LINE=" ls /usr/local/lib/R/site-library/rlib/3.0/site-library"


DOCKER_CONTAINER=genepattern/docker-r-3-2

JOB_DEFINITION_NAME="R32_Generic"
JOB_ID=gp_job_r32_Revealer_$1
JOB_QUEUE=TedTest

EXTRA_LOCAL=" -v $RLIB:$RLIB"


echo "extra mount is " $EXTRA_LOCAL
