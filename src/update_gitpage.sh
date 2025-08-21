#!/bin/bash 
cp *_pipeline.html index.html
git add .
git commit -m "$1"
git push