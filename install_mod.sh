#! /bin/bash

# This script installs the workflow as a submodule of a dataset

repo=git@github.com:swarbricklab/imputation.git
name=imputation

echo "Installing $name as git submodule"
git submodule add $repo modules/$name 
cd modules/$name
git submodule init
git submodule update
cd -

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
