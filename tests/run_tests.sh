#! /bin/bash

# Run all available tests and generate a json-report for further analysis

NUM_TEST_CASES=1000

TS=$(date +%Y%m%d%H%M%S)
HASH=$(git rev-parse --short HEAD)
BRANCH=$(git rev-parse --abbrev-ref HEAD)
REPORT_NAME="./report_${BRANCH}_${HASH}_${TS}.json"

pytest --num_test_cases=$NUM_TEST_CASES --json-report --json-report-file=$REPORT_NAME
