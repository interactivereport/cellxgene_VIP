#!/usr/bin/env bash
set -e

if [[ "$#" -lt 2 ]];then
  echo -e "\n\tslimH5ad.sh <path/to/a/h5ad/file> <path/to/save/the/slimmed/h5ad/file>\n"
else
  if [[ ! -f "$(dirname $0)/.env" ]]; then
    echo "Missing .env, please let admin know the setup is not completed"
    exit 1
  fi
  source $(dirname $0)/.env
  if [[ -n "$VIPenv" ]]; then
    eval $VIPenv >/dev/null 2>&1
    python -u $(dirname $0)/slimH5ad.py "$@"
  elif [[ -n "$scRNAview_dockerName" ]]; then
    docker exec $scRNAview_dockerName python scRNAview/bin/slimH5ad.py $@
  fi
fi





  
