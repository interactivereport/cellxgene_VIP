#!/usr/bin/env bash
source $(dirname $0)/.env
eval $VIPenv
python -u $(dirname $0)/plotH5ad.py "$@"
#echo $(dirname $0)/plotH5ad.py