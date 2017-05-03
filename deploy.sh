#!/usr/bin/env bash

chmod +x bin/*

export PATH="$PATH:/home/ionadmin/.local/bin:/results/plugins/gapsInCoverage/bin"


cd /results/plugins/gapsInCoverage

virtualenv gapsInCoverage_venv

source /results/plugins/gapsInCoverage/gapsInCoverage_venv/bin/activate

pip install --upgrade pip

pip install -r requirements.txt

git clone https://github.com/sch-sdgs/SDGSOpenCGA.git

cd SDGSOpenCGA

git checkout develop

python setup.py install

pip install -r requirements.txt

cd /results/plugins/gapsInCoverage

git clone https://github.com/sch-sdgs/SDGSCommonLibs.git

cd SDGSCommonLibs

git checkout master

python setup.py install

pip install -r requirements.txt

cd /results/plugins/gapsInCoverage

git clone https://github.com/sch-sdgs/SDGSDataModels.git

cd SDGSDataModels

git checkout master

python build.py

python setup.py install

pip install -r requirements.txt

cd /results/plugins/gapsInCoverage

pip install suds

version=`git --git-dir="/results/plugins/gapsInCoverage/.git" describe --tags|sed ':a;N;$!ba;s/\n/ /g'`

mv gapsInCoverage.py gapsInCoverage.py.old

echo ${version}

sed "s|version = '1.0'|version = '$version'|g" gapsInCoverage.py.old > gapsInCoverage.py


