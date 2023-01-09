#!/bin/bash

echo "------ 'virtual_environment.sh' started ------"

echo -n " creating virtual environment... "
python3 -m venv venv
echo "done"

echo -n " activating virtual environment... "
source venv/bin/activate
echo "done"

echo -n " installing dependencies from 'requirements.txt'... "
pip install -q -r ./scripts/requirements.txt
echo "done"

echo -n " deactivating virtual environment... "
deactivate
echo "done"

echo "------ script finished ------"
