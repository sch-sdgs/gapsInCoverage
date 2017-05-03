#!/usr/bin/env bash

cd /results/plugins/gapsInCoverage

virtualenv gapsInCoverage_venv

source /results/plugins/gapsInCoverage/gapsInCoverage_venv/bin/activate

pip install -r requirements.txt

version=`git --git-dir="/results/plugins/gapsInCoverage/.git" describe --tags|sed ':a;N;$!ba;s/\n/ /g'`

mv gapsInCoverage.py gapsInCoverage.py.old

echo ${version}

sed "s|version = '1.0'|version = '$version'|g" gapsInCoverage.py.old > gapsInCoverage.py


