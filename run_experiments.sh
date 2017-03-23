#!/bin/bash

rm config.py

for f in config*.py
do
  cp $f config.py
  echo $f
  rm config.py
  ./run.sh -fpgcde
done