#! /bin/bash

# Run regression tests

cd regression
pytest --workers auto
