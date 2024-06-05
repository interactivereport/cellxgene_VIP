#!/usr/bin/env bash
set -e
source $(dirname $0)/.env
eval $VIPenv >/dev/null 2>&1
python -u $(dirname $0)/getH5adInfo.py "$@"
