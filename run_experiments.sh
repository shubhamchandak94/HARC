#!/bin/bash

rm config.py

for f in config*.py
do
  cp $f config.py
  echo $f
  ./run.sh -fpgcde
  rm config.py
done
