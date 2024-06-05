#!/usr/bin/env bash

if [[ "$#" -ge 4 ]]; then
  set -e
  source $(dirname $0)/.env
  eval $VIPenv >/dev/null 2>&1
  python -W ignore -u $(dirname $0)/plotH5adPlotly.py "$@"
else
  echo "==============="
  echo "plotH5adPlotly path/to/h5ad plot/type [other arguments]"
  echo ""
  echo "    plotH5adPlotly path/to/h5ad violin a/gene/name an/annotation/name -n cell/number -g gene/cutoff -v <gene/column/header>"
  echo "    plotH5adPlotly path/to/h5ad dot gene/names an/annotation/name -n cell/number -g gene/cutoff -e min,max/exp/scale -p max/percentage/scale -l max/value/log -v <gene/column/header>"
  echo "Gene names are NOT case-sensitive, and separated by (,)"
  echo "Return a html string which can be inserted into a DIV tag"
  echo "==============="
fi
