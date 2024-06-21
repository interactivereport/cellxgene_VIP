#!/usr/bin/env bash
set -e

if [[ "$#" -lt 2 ]];then
  echo -e "\n\tslimH5ad.sh <path/to/a/h5ad/file> <path/to/save/the/slimmed/h5ad/file>\n"
else
  source $(dirname $0)/.env
  eval $VIPenv >/dev/null 2>&1
  python -u $(dirname $0)/slimH5ad.py "$@"
fi
