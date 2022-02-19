#! /bin/bash

# Run all available stability tests and generate a json-report for further analysis

NUM_TEST_CASES=1000 # number of etst cases per each class

REPORT_DIR="./reports"
TS=$(date +%Y%m%d%H%M%S)
HASH=$(git rev-parse --short HEAD)
BRANCH=$(git rev-parse --abbrev-ref HEAD)
REPORT_NAME="${REPORT_DIR}/report_${BRANCH}_${HASH}_${TS}.json"

cd stability
mkdir -p $REPORT_DIR
pytest --num_test_cases=$NUM_TEST_CASES --workers auto --json-report --json-report-file=$REPORT_NAME
