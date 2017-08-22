#!/bin/sh

TASKLIB=$PWD/src
INPUT_FILE_DIRECTORIES=$PWD/data
S3_ROOT=s3://moduleiotest
WORKING_DIR=$PWD/job_1111
RLIB=$PWD/rlib
RHOME=/packages/R-3.2.5/

COMMAND_LINE="Rscript --no-save --quiet --slave --no-restore $TASKLIB/correlated-features.r --libdir $TASKLIB --min.correlation 0.4"

DOCKER_CONTAINER=genepattern/docker-r-3-2

JOB_DEFINITION_NAME="R32_Generic"
JOB_ID=gp_job_r32_CoClustering_$1
JOB_QUEUE=TedTest


