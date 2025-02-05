#! /bin/bash

# This script installs the workflow as a submodule of a dataset

name=imputation

echo "Preparing config templates"
mkdir -p config/$name
cp modules/$name/config/template*.yaml config/$name
echo "Edit config templates to match your dataset"

echo "Preparing resources"
mkdir -p resources
cp --no-clobber -a modules/$name/resources/. resources/
echo "Resource files copied to dataset"
echo "Commit .dvc and .gitignore files to the dataset repo"
echo "Then 'dvc pull -R resources/' to download the resources to the dataset"
