#!/bin/bash

echo "------ script started ------"

echo "-- creating flask virtual environment --"

cd flask

python3 -m venv venv

source venv/bin/activate

pip install -r requirements.txt

deactivate

cd ..

echo "-- creating jupyter-notebook virtual environment --"

cd jupyter-notebook

python3 -m venv venv

source venv/bin/activate

pip install -r requirements.txt

deactivate

cd ..

echo "------ script finished ------"
