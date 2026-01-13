#!/usr/bin/env bash
set -e

if [[ ! -f "$(dirname $0)/.env" ]]; then
  echo "Missing .env, please let admin know the setup is not completed"
  exit 1
fi
source $(dirname $0)/.env
if [[ -n "$VIPenv" ]]; then
  eval $VIPenv >/dev/null 2>&1
  python -u $(dirname $0)/getH5adInfo.py "$@"
elif [[ -n "$scRNAview_dockerName" ]]; then
  docker exec $scRNAview_dockerName python scRNAview/bin/getH5adInfo.py $@
fi
